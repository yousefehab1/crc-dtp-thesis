# ==============================================================================
# core/signatures.R  —  ONE gene-set loader for all modules (#1).
# sig.csv holds every signature (one column each). Each module selects which
# columns it needs by name; pan-cancer just selects "Up"/"Down" and derives
# Composite after scoring (see core/scoring.R).
# ==============================================================================

# Load sig.csv into a named list of character vectors (one per column).
load_signature_panel <- function(path = SIG_FILE) {
  stopifnot(file.exists(path))
  # Read with check.names=FALSE and strip a leading UTF-8 BOM from the header
  # BEFORE make.names(): otherwise a BOM + "Down" header sanitises to "X...Down",
  # silently renaming that signature. (fileEncoding="UTF-8-BOM" is NOT usable here
  # — it truncates this file at its first non-UTF-8 byte, dropping most rows.)
  df <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  names(df) <- make.names(sub("^﻿", "", names(df)))
  sets <- lapply(colnames(df), function(cn) {
    g <- unique(as.character(df[[cn]]))
    g <- g[!is.na(g) & trimws(g) != "" & !tolower(g) %in% c("na", "none")]
    g
  })
  names(sets) <- colnames(df)
  sets <- sets[lengths(sets) > 0]
  message(sprintf("Signature panel: %d sets loaded (%s).",
                  length(sets), paste(names(sets), collapse = ", ")))
  sets
}

# Build a GeneSetCollection from selected panel entries (default: all).
build_gene_sets <- function(panel, which = NULL) {
  if (!is.null(which)) {
    miss <- setdiff(which, names(panel))
    if (length(miss)) stop("Signatures not found in panel: ", paste(miss, collapse = ", "))
    panel <- panel[which]
  }
  GSEABase::GeneSetCollection(
    lapply(names(panel), function(nm) GSEABase::GeneSet(panel[[nm]], setName = nm))
  )
}

# Up/Down overlap warning (a gene in both contributes +/- simultaneously).
check_updown_overlap <- function(panel, up = "Up", down = "Down") {
  if (all(c(up, down) %in% names(panel))) {
    ov <- intersect(panel[[up]], panel[[down]])
    if (length(ov))
      warning(length(ov), " gene(s) in both ", up, " and ", down, ": ",
              paste(head(ov, 10), collapse = ", "), call. = FALSE)
  }
}
