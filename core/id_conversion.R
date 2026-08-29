# ==============================================================================
# core/id_conversion.R  —  Gene-ID type detection, conversion, and panel
# unification.  Controlled by ID_TYPE in config.R ("symbol" or "ensembl").
#
# All conversions use org.Hs.eg.db (no internet required) and are cached on
# disk via cache_rds() so they only run once per project.
# ==============================================================================

# --- Heuristic: what ID type does a character vector contain? -----------------
# Returns "ensembl" if >80 % of non-NA entries start with ENSG, else "symbol".
detect_id_type <- function(genes) {
  g <- genes[!is.na(genes) & nzchar(genes)]
  if (length(g) == 0) return("symbol")
  if (mean(grepl("^ENSG", g)) > 0.8) "ensembl" else "symbol"
}

# --- Strip Ensembl version suffix (ENSG00000123456.3 -> ENSG00000123456) ------
strip_ensembl_version <- function(ids) sub("\\.\\d+$", "", ids)

# --- One-way conversion via org.Hs.eg.db (cached) -----------------------------
# We cache the COMPLETE from→to mapping for the entire organism once, then do
# per-call subset lookups in memory.  This avoids the previous bug where a
# fixed cache key caused different gene-set calls to all return the same result.
.get_full_id_map <- function(from, to) {
  if (!requireNamespace("org.Hs.eg.db", quietly = TRUE))
    stop("org.Hs.eg.db is required for ID conversion. ",
         "Install with: BiocManager::install('org.Hs.eg.db')")
  key <- paste0("full_id_map_", tolower(from), "_to_", tolower(to))
  cache_rds(key, function() {
    all_keys <- AnnotationDbi::keys(org.Hs.eg.db::org.Hs.eg.db, keytype = from)
    map <- suppressMessages(
      AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db,
                            keys      = all_keys,
                            column    = to,
                            keytype   = from,
                            multiVals = "first")
    )
    map[!is.na(map)]   # named vector: from-ID → to-ID, NAs dropped
  })
}

# Public conversion function. `from` / `to` are AnnotationDbi keytype strings
# ("SYMBOL" or "ENSEMBL"). Uses the cached full map for O(1) per-gene lookups.
convert_gene_ids <- function(genes, from = "SYMBOL", to = "ENSEMBL") {
  map         <- .get_full_id_map(from, to)
  genes_clean <- unique(strip_ensembl_version(genes))
  result      <- map[genes_clean]          # named, may contain NAs for misses
  n_lost      <- sum(is.na(result))
  if (n_lost > 0)
    message(sprintf("  ID conversion %s->%s: %d/%d genes unmapped (dropped).",
                    from, to, n_lost, length(genes_clean)))
  unique(unname(result[!is.na(result)]))
}

# --- Unify every signature in the panel to ID_TYPE ---------------------------
# Signatures already in the target type are returned unchanged.
# Mixed panels (some ensembl, some symbol) are handled gene-set by gene-set.
unify_panel_ids <- function(panel, target = ID_TYPE) {
  if (!target %in% c("symbol", "ensembl"))
    stop("ID_TYPE must be 'symbol' or 'ensembl', got: ", target)

  from_col <- if (target == "ensembl") "SYMBOL"  else "ENSEMBL"
  to_col   <- if (target == "ensembl") "ENSEMBL" else "SYMBOL"

  converted <- lapply(names(panel), function(nm) {
    genes    <- panel[[nm]]
    detected <- detect_id_type(genes)
    if (detected == target) {
      # Already correct type — just strip version suffixes if ensembl.
      if (target == "ensembl") strip_ensembl_version(genes) else genes
    } else {
      message(sprintf("  [id_conversion] '%s': %s -> %s", nm, detected, target))
      convert_gene_ids(genes, from = from_col, to = to_col)
    }
  })
  names(converted) <- names(panel)
  converted <- converted[lengths(converted) > 0]
  message(sprintf("Panel unified to '%s': %d signatures retained.", target, length(converted)))
  converted
}

# --- Expression-matrix row-name harmonisation ---------------------------------
# Call this after building any expression matrix to ensure its row IDs match
# ID_TYPE.  `current_type` is auto-detected if not supplied.
harmonise_matrix_ids <- function(mat, current_type = NULL) {
  if (is.null(current_type)) current_type <- detect_id_type(rownames(mat))
  target <- ID_TYPE
  if (current_type == target) {
    if (target == "ensembl") rownames(mat) <- strip_ensembl_version(rownames(mat))
    return(mat)
  }
  from_col <- if (target == "ensembl") "SYMBOL"  else "ENSEMBL"
  to_col   <- if (target == "ensembl") "ENSEMBL" else "SYMBOL"

  orig_ids <- strip_ensembl_version(rownames(mat))
  new_ids  <- suppressMessages(
    AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db,
                          keys = orig_ids, column = to_col, keytype = from_col,
                          multiVals = "first")
  )

  # Alias fallback (symbol-keyed matrices only). Older datasets label rows with
  # symbols HGNC has since renamed (e.g. GSE50760, 2014: SDPR->CAVIN2,
  # HIST1H4H->H4C8, PVRL4->NECTIN4). keytype "SYMBOL" holds only the CURRENT
  # symbol, so those rows drop out even though the gene is measured. For the
  # unmapped residue we retry via keytype "ALIAS", accepting only aliases that
  # resolve to exactly ONE target id — ambiguous aliases (an old symbol shared
  # by several genes) stay dropped rather than risk a wrong assignment.
  if (from_col == "SYMBOL" && anyNA(new_ids)) {
    miss  <- which(is.na(new_ids))
    resid <- orig_ids[miss]
    al    <- suppressMessages(
      AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db,
                            keys = unique(resid), column = to_col,
                            keytype = "ALIAS", multiVals = "list"))
    al1   <- vapply(al, function(x) { x <- x[!is.na(x)]
                    if (length(x) == 1L) x else NA_character_ }, character(1))
    filled <- unname(al1[resid])
    new_ids[miss] <- filled
    n_rec <- sum(!is.na(filled))
    if (n_rec > 0)
      message(sprintf("  harmonise_matrix_ids: +%d rows recovered via ALIAS (renamed symbols).",
                      n_rec))
  }

  keep          <- !is.na(new_ids)
  n_lost        <- sum(!keep)
  if (n_lost > 0)
    message(sprintf("  harmonise_matrix_ids %s->%s: %d/%d rows unmapped (dropped).",
                    current_type, target, n_lost, nrow(mat)))
  mat           <- mat[keep, , drop = FALSE]
  rownames(mat) <- unname(new_ids[keep])
  limma::avereps(mat)   # collapse any duplicates introduced by many-to-one mapping
}
