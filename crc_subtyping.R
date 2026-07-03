# ==============================================================================
# crc_subtyping.R  —  CMS + PDS molecular subtyping (standalone runner).
#
# Shares ONE implementation with the integrated CRC survival module: the actual
# classifiers live in core/subtyping.R (call_cms / call_pds / build_symbol_matrix).
# This script only downloads the two CRC cohorts and delegates. It can be run on
# its own for subtype calls without the survival analysis.
#
#   source("core/config.R"); source("core/io.R"); source("core/expression.R")
#   source("core/scoring.R"); source("core/subtyping.R")
#   root <- init_run("subtyping")
#   run_crc_subtyping(root)
#
# PDS requires the GitHub-only PDSclassifier package. Pinned to the exact commit that
# produced the reported PDS calls (v1.0.1) for reproducibility:
#   remotes::install_github("sidmall/PDSclassifier@c89a19c891185f7806848a990571423d6a32d2a8")
# If absent, PDS is skipped and CMS still runs.
#
# NOTE: confirm PDSpredict input orientation against data(testData) in the
# package before relying on PDS calls (README is thin on matrix format).
# ------------------------------------------------------------------------------

run_crc_subtyping <- function(out_root) {
  module <- "crc_subtyping"
  out    <- file.path(out_root, module)
  dir.create(out, recursive = TRUE, showWarnings = FALSE)

  # --- GSE39582 (microarray) --------------------------------------------------
  eset    <- cache_rds("GSE39582_eset", function() { g <- GEOquery::getGEO("GSE39582", GSEMatrix = TRUE); g[[1]] })
  mat_gse <- build_symbol_matrix(Biobase::exprs(eset), "microarray")

  # --- TCGA-COAD (TPM) --------------------------------------------------------
  query <- TCGAbiolinks::GDCquery(project = "TCGA-COAD",
                                  data.category = "Transcriptome Profiling",
                                  data.type     = "Gene Expression Quantification",
                                  workflow.type = "STAR - Counts",
                                  sample.type   = "Primary Tumor")
  se       <- cache_rds("TCGA_COAD_se", function() { TCGAbiolinks::GDCdownload(query); TCGAbiolinks::GDCprepare(query) })
  mat_tcga <- build_symbol_matrix(se, "tcga")

  # --- CMS + PDS via shared callers (RNAseq=FALSE, log2 symbol matrices) ------
  tag <- function(df, cohort) { if (is.null(df)) return(NULL); df$Cohort <- cohort; df }

  cms <- dplyr::bind_rows(
    tag(tryCatch(call_cms(mat_gse),  error = function(e) { message("CMS GSE: ",  e$message); NULL }), "GSE39582"),
    tag(tryCatch(call_cms(mat_tcga), error = function(e) { message("CMS TCGA: ", e$message); NULL }), "TCGA-COAD"))
  if (!is.null(cms) && nrow(cms) > 0)
    write.csv(cms, file.path(out, "CMS_calls.csv"), row.names = FALSE)

  pds <- dplyr::bind_rows(
    tag(call_pds(mat_gse),  "GSE39582"),
    tag(call_pds(mat_tcga), "TCGA-COAD"))
  if (!is.null(pds) && nrow(pds) > 0)
    write.csv(pds, file.path(out, "PDS_calls.csv"), row.names = FALSE)

  message("Subtyping complete (standalone): CMS + PDS.")
}
