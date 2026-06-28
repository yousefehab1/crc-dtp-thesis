# ==============================================================================
# subtyping/crc_subtyping.R  —  CMS / CRIS / PDS molecular subtyping.
#
# PULLED OUT of the main CRC survival module per decision #11. It is NOT sourced
# or run by main.R. It is kept here, standalone, because its fate is undecided:
# it may be re-integrated, kept as a side analysis, or removed entirely. Run it
# only on its own when you want subtype calls for the two CRC cohorts.
#
#   source("core/config.R"); source("core/io.R"); source("core/expression.R")
#   root <- init_run("subtyping")
#   run_crc_subtyping(root)
#
# ------------------------------------------------------------------------------
# !!! OPEN §11 ITEMS TO RESOLVE BEFORE TRUSTING THESE CALLS (deferred, to discuss)
# !!!  (a) ematAdjust asymmetry: the original applied median-centering for the
# !!!      CMS classifier on ONE cohort but not the other. Decide whether both
# !!!      cohorts should be centered identically.
# !!!  (b) CMScaller RNA-seq flag: passing a log2(TPM+1) matrix while also
# !!!      letting CMScaller apply its own RNA-seq transform risks a double
# !!!      transform. Confirm the intended input scale for each classifier.
# !!!  (c) PDS normalization is undocumented in the original; document the exact
# !!!      expected input (raw/log/z-scored) before relying on PDS calls.
# ------------------------------------------------------------------------------

run_crc_subtyping <- function(out_root) {
  module <- "crc_subtyping"
  dir.create(file.path(out_root, module), recursive = TRUE, showWarnings = FALSE)

  # --- GSE39582 (microarray) --------------------------------------------------
  eset    <- cache_rds("GSE39582_eset", function() { g <- GEOquery::getGEO("GSE39582", GSEMatrix = TRUE); g[[1]] })
  mat_gse <- prep_microarray_symbols(Biobase::exprs(eset))

  # --- TCGA-COAD (TPM) --------------------------------------------------------
  query <- TCGAbiolinks::GDCquery(project = "TCGA-COAD",
                                  data.category = "Transcriptome Profiling",
                                  data.type     = "Gene Expression Quantification",
                                  workflow.type = "STAR - Counts",
                                  sample.type   = "Primary Tumor")
  se       <- cache_rds("TCGA_COAD_se", function() { TCGAbiolinks::GDCdownload(query); TCGAbiolinks::GDCprepare(query) })
  mat_tcga <- prep_tcga_tpm(se)

  # --- CMS (CMScaller) --------------------------------------------------------
  # NOTE (§11b): verify RNA-seq flag vs. input scale before trusting calls.
  cms_call <- function(mat, label, is_rnaseq) {
    res <- CMScaller::CMScaller(emat = mat, RNAseq = is_rnaseq, doPlot = FALSE)
    data.frame(Sample_ID = rownames(res), CMS = res$prediction, Cohort = label, stringsAsFactors = FALSE)
  }
  cms_gse  <- tryCatch(cms_call(mat_gse,  "GSE39582",  is_rnaseq = FALSE), error = function(e) { message("CMS GSE: ", e$message); NULL })
  cms_tcga <- tryCatch(cms_call(mat_tcga, "TCGA-COAD", is_rnaseq = TRUE),  error = function(e) { message("CMS TCGA: ", e$message); NULL })

  # --- CRIS (CMScaller::CRIS or templates) ------------------------------------
  cris_call <- function(mat, label) {
    res <- CMScaller::CMScaller(emat = mat, templates = CMScaller::geneSets.CRIS, RNAseq = TRUE, doPlot = FALSE)
    data.frame(Sample_ID = rownames(res), CRIS = res$prediction, Cohort = label, stringsAsFactors = FALSE)
  }
  cris_gse  <- tryCatch(cris_call(mat_gse,  "GSE39582"),  error = function(e) { message("CRIS GSE: ", e$message); NULL })
  cris_tcga <- tryCatch(cris_call(mat_tcga, "TCGA-COAD"), error = function(e) { message("CRIS TCGA: ", e$message); NULL })

  # --- PDS --------------------------------------------------------------------
  # NOTE (§11c): document the exact expected input scale before relying on this.
  # Placeholder hook left intentionally — wire to the original PDS routine once
  # the normalization question is resolved.
  message("PDS: deferred pending normalization decision (§11c).")

  out <- file.path(out_root, module)
  if (!is.null(cms_gse) || !is.null(cms_tcga))
    write.csv(dplyr::bind_rows(cms_gse, cms_tcga),   file.path(out, "CMS_calls.csv"),  row.names = FALSE)
  if (!is.null(cris_gse) || !is.null(cris_tcga))
    write.csv(dplyr::bind_rows(cris_gse, cris_tcga), file.path(out, "CRIS_calls.csv"), row.names = FALSE)

  message("Subtyping complete (standalone). Review §11 open items before use.")
}
