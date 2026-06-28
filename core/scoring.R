# ==============================================================================
# core/scoring.R  —  ONE ssGSEA wrapper + ONE score-naming convention (#12).
# Every score column is "<Signature>_ssGSEA". Composite is derived AFTER
# scoring (it's a score operation, not a gene-set), keeping #1 clean.
# ==============================================================================

# Score a log2-scale matrix; returns a data frame (samples x signatures),
# rownames = sample/barcode, columns "<Signature>_ssGSEA".
run_ssgsea <- function(mat, gene_sets, suffix = SCORE_SUFFIX) {
  param <- GSVA::ssgseaParam(exprData = as.matrix(mat), geneSets = gene_sets)
  res   <- GSVA::gsva(param)
  df    <- as.data.frame(t(res))
  colnames(df) <- paste0(colnames(df), suffix)
  df
}

# Composite = Up - Down, on the scored columns.
add_composite <- function(scores_df, up = "Up", down = "Down", suffix = SCORE_SUFFIX) {
  uc <- paste0(up, suffix); dc <- paste0(down, suffix)
  if (all(c(uc, dc) %in% colnames(scores_df)))
    scores_df[[paste0("Composite", suffix)]] <- scores_df[[uc]] - scores_df[[dc]]
  scores_df
}

# TCGA barcode -> 12-char patient ID.
tcga_patient_id <- function(barcodes)
  stringr::str_extract(gsub("\\.", "-", toupper(barcodes)), "TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}")
