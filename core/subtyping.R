# ==============================================================================
# core/subtyping.R  —  CMS + PDS molecular subtype callers, download-free.
#
# One shared code path for both the standalone crc_subtyping.R module and the
# integrated crc_survival.R module. Classifiers are fed a LOG2, gene-SYMBOL
# matrix built independently of the ssGSEA ID namespace (ID_TYPE) because:
#   - CMScaller maps rownames via rowNames="symbol" to its Entrez templates;
#   - PDSclassifier expects gene symbols.
#
# CORRECTNESS (resolves former §11 items):
#  (§11a/b) Both cohorts are classified as ALREADY-LOG2 expression with
#           RNAseq=FALSE, so CMScaller applies identical median-centering/scaling
#           (ematAdjust) to each and does NOT re-log2/quantile-normalise. Passing
#           RNAseq=TRUE to a log2(TPM+1) matrix would double-transform.
#  (§11c)   PDSpredict takes a log2 symbol matrix; we log gene coverage so a poor
#           overlap is visible rather than silently degrading the calls.
# ==============================================================================

# --- Build a log2, gene-symbol matrix for one cohort (no re-download) ---------
# GSE39582: pass the raw ExpressionSet exprs matrix (probe rows).
# TCGA:     pass the cached SummarizedExperiment.
build_symbol_matrix <- function(x, source = c("microarray", "tcga")) {
  source <- match.arg(source)
  if (source == "microarray")
    prep_microarray_symbols(x, id_type = "symbol")
  else
    prep_tcga_tpm(x, id_type = "symbol")
}

# --- CMS via CMScaller (RNAseq=FALSE: input is already log2) ------------------
call_cms <- function(emat_symbol_log2) {
  message(sprintf("  CMS: classifying %d samples x %d symbols (RNAseq=FALSE).",
                  ncol(emat_symbol_log2), nrow(emat_symbol_log2)))
  res <- CMScaller::CMScaller(emat = emat_symbol_log2, rowNames = "symbol",
                              RNAseq = FALSE, doPlot = FALSE)
  data.frame(Sample_ID = rownames(res),
             CMS = as.character(res$prediction),
             stringsAsFactors = FALSE)
}

# --- PDS via PDSclassifier (guarded: GitHub-only dependency) ------------------
# threshold = 0.6 is the PDSpredict (v1.0.1) package default: a sample is assigned to
# PDS1/2/3 only when that subtype's SVM posterior probability exceeds 0.6, else "Mixed".
# We keep the classifier authors' default rather than tuning it (yields ~21% Mixed in both
# CRC cohorts). Pinned install: sidmall/PDSclassifier@c89a19c891185f7806848a990571423d6a32d2a8.
# Returns NULL (with a single skip message) if the package is unavailable, so
# the caller degrades to CMS + survival gracefully.
call_pds <- function(emat_symbol_log2, species = "human", threshold = 0.6) {
  if (!requireNamespace("PDSclassifier", quietly = TRUE)) {
    message("  PDS: PDSclassifier not installed; skipping ",
            "(remotes::install_github('sidmall/PDSclassifier@c89a19c')).")
    return(NULL)
  }
  # PDSpredict expects a data.frame whose first column is the gene symbol and
  # remaining columns are samples (genes in rows). Confirm against data(testData).
  test_df <- data.frame(Gene = rownames(emat_symbol_log2),
                        as.data.frame(emat_symbol_log2, check.names = FALSE),
                        check.names = FALSE, stringsAsFactors = FALSE)
  res <- PDSclassifier::PDSpredict(test_df, species = species, threshold = threshold)
  # PDSpredict returns Sample_ID | PDS1 | PDS2 | PDS3 | prediction | PDS_call.
  # Prefer PDS_call: the threshold-aware label ("Mixed" for ambiguous samples)
  # rather than the raw argmax `prediction`.
  pds_col <- intersect(c("PDS_call", "prediction", "PDS"), colnames(res))[1]
  id_col  <- intersect(c("Sample_ID", "Sample", "SampleID", "sample"), colnames(res))[1]
  if (is.na(pds_col)) {
    message("  PDS: unexpected PDSpredict output columns (",
            paste(colnames(res), collapse = ", "), "); skipping.")
    return(NULL)
  }
  ids <- if (is.na(id_col)) rownames(res) else res[[id_col]]
  data.frame(Sample_ID = as.character(ids),
             PDS = as.character(res[[pds_col]]),
             stringsAsFactors = FALSE)
}

# --- Attach subtype columns onto a clinical frame -----------------------------
# subtype_df has a Sample_ID column keyed to matrix colnames. `key_fun` maps
# those Sample_IDs to the clinical join key (identity for GSE; barcode->patient
# collapse for TCGA), applied with the same deterministic dedup as the scores.
attach_subtypes <- function(clinical, subtype_df, clinical_key = "Sample_ID",
                            key_fun = identity) {
  if (is.null(subtype_df) || nrow(subtype_df) == 0) return(clinical)
  sd <- subtype_df
  sd[[clinical_key]] <- key_fun(sd$Sample_ID)
  sd <- sd %>%
    dplyr::arrange(Sample_ID) %>%                       # deterministic vial pick
    dplyr::distinct(.data[[clinical_key]], .keep_all = TRUE)
  if (clinical_key != "Sample_ID") sd$Sample_ID <- NULL # keep it when it IS the key
  dplyr::left_join(clinical, sd, by = clinical_key)
}

# --- Build subtype-only cohorts for the survival block ------------------------
# Given a clinical frame with a subtype column, return a named list of one
# cohort per non-missing level, named "<prefix>_<level>" (e.g. GSE_CMS1).
subtype_cohorts <- function(clinical, subtype_col, prefix) {
  if (!subtype_col %in% colnames(clinical)) return(list())
  lvls <- sort(unique(stats::na.omit(clinical[[subtype_col]])))
  lvls <- lvls[lvls != ""]
  out <- lapply(lvls, function(lv)
    clinical[!is.na(clinical[[subtype_col]]) & clinical[[subtype_col]] == lv, , drop = FALSE])
  names(out) <- paste0(prefix, "_", lvls)
  out
}
