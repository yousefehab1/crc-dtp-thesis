# ==============================================================================
# core/expression.R  —  Per-platform normalisation, all converging on a
# log2-scale, gene-symbol-indexed matrix ready for ssGSEA (#2, #3).
#
# Convention now shared by every module: ssGSEA is fed a LOG2-scale matrix.
# (ssGSEA is a within-sample rank statistic, so log2 vs linear is rank-
# invariant — but we standardise on log2 everywhere so the inputs are
# literally one convention, and we drop the pan-cancer back-transform (#2b).)
# ==============================================================================

# --- Affymetrix microarray (GSE39582): probe -> symbol or ensembl, log2 ------
# Target ID type is governed by ID_TYPE (config.R). When ID_TYPE == "ensembl"
# the function maps probes directly to ENSEMBL IDs; when "symbol" it behaves
# identically to the original implementation.
prep_microarray_symbols <- function(exp_matrix, anno_db = NULL, id_type = ID_TYPE) {
  if (is.null(anno_db)) anno_db <- hgu133plus2.db::hgu133plus2.db
  target_col <- if (id_type == "ensembl") "ENSEMBL" else "SYMBOL"
  id_label   <- if (id_type == "ensembl") "ensembl" else "symbol"

  probe_ids <- AnnotationDbi::mapIds(anno_db, keys = rownames(exp_matrix),
                                     column = target_col, keytype = "PROBEID",
                                     multiVals = "first")
  probe_all <- AnnotationDbi::mapIds(anno_db, keys = rownames(exp_matrix),
                                     column = target_col, keytype = "PROBEID",
                                     multiVals = "CharacterList")
  n_multi <- sum(lengths(probe_all) > 1)
  message(sprintf("  probe annotation: %d/%d (%.1f%%) multi-map; multiVals='first'.",
                  n_multi, length(probe_all), 100 * n_multi / length(probe_all)))

  mat <- exp_matrix
  rownames(mat) <- unname(probe_ids)
  n_unmapped <- sum(is.na(rownames(mat)))
  message(sprintf("  probe->%s: %d/%d (%.1f%%) unmapped, dropped.",
                  id_label, n_unmapped, nrow(mat), 100 * n_unmapped / nrow(mat)))
  mat <- mat[!is.na(rownames(mat)), , drop = FALSE]
  if (id_type == "ensembl")
    rownames(mat) <- strip_ensembl_version(rownames(mat))
  message(sprintf("  collapsing duplicates: %d probes -> %d %ss.",
                  nrow(mat), length(unique(rownames(mat))), id_label))
  limma::avereps(mat)
}

# --- TCGA STAR TPM (CRC-COAD and pan-cancer): log2(TPM+1), log-space collapse,
#     primary-only filter (#2a, #3-TCGA, #17). Used identically by both arms.
# Row names follow ID_TYPE: gene_name (symbols) or gene_id (Ensembl, version
# suffix stripped) from rowData.
prep_tcga_tpm <- function(se, tumour_codes = TUMOUR_CODES, id_type = ID_TYPE) {
  if (!"tpm_unstrand" %in% SummarizedExperiment::assayNames(se))
    stop("Assay 'tpm_unstrand' not found. Available: ",
         paste(SummarizedExperiment::assayNames(se), collapse = ", "))

  rd  <- SummarizedExperiment::rowData(se)
  mat <- SummarizedExperiment::assay(se, "tpm_unstrand")

  if (id_type == "ensembl") {
    if (!"gene_id" %in% colnames(rd))
      stop("rowData has no 'gene_id' column; cannot use id_type='ensembl' for TCGA.")
    rownames(mat) <- strip_ensembl_version(rd$gene_id)
  } else {
    rownames(mat) <- rd$gene_name
  }

  mat <- mat[!is.na(rownames(mat)), , drop = FALSE]
  mat <- mat[rowSums(is.na(mat)) == 0, , drop = FALSE]

  mat <- log2(mat + 1)                      # standardise on log2 (#2b)
  mat <- limma::avereps(mat)                # geometric-mean collapse, matches pan-cancer (#3)

  # primary-only sample filter (#17)
  cd   <- as.data.frame(SummarizedExperiment::colData(se))
  keep <- cd$sample_type %in% tumour_codes
  if (sum(keep) == 0) stop("No samples of type(s): ", paste(tumour_codes, collapse = ", "))
  mat[, keep, drop = FALSE]
}

# --- Gene-set coverage check (pan-cancer gate) --------------------------------
gene_set_coverage <- function(mat, genes) mean(genes %in% rownames(mat)) * 100
