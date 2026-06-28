# ==============================================================================
# main.R  —  Orchestrator for the CRC DTP pipeline.
#
# Sources the shared core, then runs the three analysis modules against one
# sig.csv panel and one timestamped output root. Subtyping is intentionally
# NOT run here (decision #11) — see subtyping/crc_subtyping.R.
#
# Required inputs in the working directory:
#   sig.csv       master signature panel (one column per signature; must include
#                 the columns named below). Pan-cancer needs "Up" and "Down".
#   TCGA-CDR.csv  TCGA Clinical Data Resource (Liu et al. 2018).
#
# Run:  Rscript main.R          (or source("main.R") in an interactive session)
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(tibble); library(stringr); library(readr)
  library(ggplot2); library(limma); library(survival)
  library(GSVA); library(GSEABase); library(patchwork)
})

# Pan-cancer downloads large SummarizedExperiment objects; TCGA-BRCA alone
# requires ~20 GB of peak RAM.  R_MAX_VSIZE must be set before the session
# starts (e.g. add R_MAX_VSIZE=32Gb to ~/.Renviron).  We emit a clear warning
# here so the cause of any OOM failure is obvious in the log.
local({
  lim_gb <- mem.maxVSize() / 1e9
  if (!is.infinite(lim_gb) && lim_gb < 24)
    warning(sprintf(
      "R memory limit is %.0f GB — large TCGA cohorts (e.g. BRCA) may fail.\n",
      lim_gb),
      "Add  R_MAX_VSIZE=32Gb  to ~/.Renviron and restart R to raise the limit.",
      call. = FALSE)
})

# ---- Source core (order matters: config -> io -> the rest) -------------------
source("core/config.R")
source("core/io.R")
source("core/id_conversion.R")   # depends on config (ID_TYPE) and io (cache_rds)
source("core/signatures.R")
source("core/expression.R")
source("core/scoring.R")
source("core/clinical.R")
source("core/stats.R")
source("core/plotting.R")
source("core/gsea_plots.R")      # depends on pub_theme (plotting.R)

# ---- Source modules ---------------------------------------------------------
source("modules/crc_survival.R")
source("modules/pancan_survival.R")
source("modules/mets_de.R")

# ---- Which signatures each analysis uses (#1) -------------------------------
# Edit these to point at the relevant sig.csv columns.
CRC_SIGNATURES  <- NULL          # NULL = use every column in sig.csv
PANCAN_UP       <- "Up"
PANCAN_DOWN     <- "Down"
METS_CORE_SIG   <- "Up"          # signature highlighted in the Mets PCA overlay

# ---- Run --------------------------------------------------------------------
out_root <- init_run("CRC_DTP")
panel    <- unify_panel_ids(load_signature_panel(SIG_FILE))   # convert to ID_TYPE

run_crc_survival(out_root, panel, crc_signatures = CRC_SIGNATURES)
run_pancan_survival(out_root, panel, up_name = PANCAN_UP, down_name = PANCAN_DOWN)
run_mets_de(out_root, panel, core_signature = METS_CORE_SIG)

finalize_run(out_root)
message("\nAll modules finished. Outputs in: ", out_root)
