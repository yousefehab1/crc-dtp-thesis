# ==============================================================================
# core/config.R  —  Single source of truth for every shared constant.
# Changing a threshold here changes it for ALL modules (decisions #9, #14, #15).
# ==============================================================================

GLOBAL_SEED    <- 42          # one global seed (#15). Only GSEA tie-breaking is stochastic.
DAYS_PER_MONTH <- 30.44       # CDR stores days; convert once at load.

# --- Landmark / survival ------------------------------------------------------
OS_CUTPOINT    <- 36          # months (3-yr landmark)
RFS_CUTPOINT   <- 36

# --- Statistical gating thresholds (#9) — applied identically everywhere ------
MIN_GROUP_N    <- 10          # Wilcoxon: minimum per group
MIN_KM_GROUP_N <- 5           # KM: minimum per High/Low group
MIN_EVENTS     <- 5           # KM/Cox: minimum events
MIN_COX_N      <- 10          # Cox: minimum n

# --- Cohort / coverage gates --------------------------------------------------
MIN_COHORT_N   <- 100         # pan-cancer: min matched patients AFTER RNA-seq join
MIN_GENE_COV   <- 70          # min % of a gene set present in the expression matrix

# --- TCGA conventions ---------------------------------------------------------
# Primary-only in BOTH CRC and pan-cancer arms (#17). The metastatic fallback
# that pan-cancer used previously is removed so the COAD overlap is defined
# identically across arms.
TUMOUR_CODES     <- c("Primary Tumor")
EXCLUDE_PROJECTS <- c("TCGA-LAML", "TCGA-DLBC")   # non-solid; DTP biology N/A

# --- Clinical source (#4) -----------------------------------------------------
# One shared CDR recurrence endpoint applied to COAD in both arms. PFI = broad
# coverage, CDR-recommended pan-cancer. Switch to "DFI" for a CRC sensitivity run.
CDR_FILE     <- "Data/TCGA-CDR.csv"
CDR_RFS_TYPE <- "PFI"

# --- Inputs / cache -----------------------------------------------------------
SIG_FILE  <- "Data/Fatemeh.csv"        # master panel; each module selects which columns it uses (#1)
CACHE_DIR <- "cache"          # shared GEO/GDC object cache (#16)

# --- Gene ID type -------------------------------------------------------------
# "symbol"  → HGNC gene symbols throughout (default, matches legacy behaviour).
# "ensembl" → Ensembl gene IDs (ENSG…) throughout.
# Changing this converts sig.csv gene sets AND expression matrix row names so
# both always use the same namespace before ssGSEA is run.
ID_TYPE <- "ensembl"

# --- Score naming (#12) -------------------------------------------------------
SCORE_SUFFIX <- "_ssGSEA"     # every ssGSEA column is "<Signature>_ssGSEA"
