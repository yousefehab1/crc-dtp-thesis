# ==============================================================================
# modules/pancan_treated.R  —  Pan-cancer DTP survival, TREATED patients only.
#
# Standalone add-on (mirrors the run_composites decoupling). Reuses the already-
# scored MASTER_Pan_Cancer_Clinical.csv from a completed pan-cancer run — there is
# NO RNA re-download or re-scoring here. "Treated" = TCGA pharmaceutical therapy
# OR radiation therapy ("yes" in either field) — the same definition the CRC
# module uses (crc_survival.R:169), so one treated definition holds across the
# whole thesis. Treatment status is pulled once per cohort via GDCquery_clinic
# (small clinical tables) and cached as the RAW yes/no flags, so the definition
# can be changed later without another network pull.
#
# Outputs (under <run>/pancan_survival/treated/):
#   treatment_status.csv                    per-cohort treated / total counts
#   treatment_status_perpatient.csv         per-patient status
#   MASTER_Pan_Cancer_Clinical_Treated.csv  treated subset of the scored master
#   FDR_Stats_Summary_Treated.csv           per-cohort + batch-corrected aggregate
#   Fig6{A,B,C}_pancan_forest_*_Treated.(png|pdf)  per-cohort Up / Composite / Down
#   Fig6_pancan_survival_composite_Treated.(png|pdf)  the three panels assembled
#   Fig7_pancan_violins_Treated.(png|pdf)      pooled score by 3-year outcome
#   Fig_Pancan_aggregate_forest.(png|pdf)      batch-corrected pooled aggregate
#   publication/                            title-stripped copies of the above
#
# Depends on modules/composites.R (build_pancan_forest / build_pancan_aggregate),
# core/stats.R (run_all_stats / apply_fdr), core/io.R (cache_rds), core/config.R.
# ==============================================================================

# Per-cohort TCGA treatment flags. Cached per project (same filename scheme as
# cache_rds) so the GDC clinical pull happens once and is incremental. We cache the
# RAW pharmaceutical / radiation yes-no flags rather than a derived label, so the
# "treated" definition can be changed downstream without re-pulling from GDC. A
# network FAILURE is NOT cached (it would otherwise pin the cohort at zero treated
# forever) — only a real response is, including a cohort that genuinely lacks one
# or both treatment fields.
.pancan_treatment_flags <- function(projects) {
  has_value <- function(col, val) vapply(col, function(x) val %in% x, logical(1))
  # A cohort can be missing either field independently; absent => FALSE both ways,
  # which lands the patient in "Unknown/Not Reported" downstream.
  flag <- function(cl, field, val)
    if (field %in% names(cl)) has_value(cl[[field]], val) else rep(FALSE, nrow(cl))
  empty <- data.frame(Patient_ID = character(), Project_ID = character(),
                      pharma_yes = logical(), pharma_no = logical(),
                      rad_yes = logical(), rad_no = logical(), stringsAsFactors = FALSE)
  PHARMA <- "treatments_pharmaceutical_treatment_or_therapy"
  RAD    <- "treatments_radiation_treatment_or_therapy"
  cfile <- function(proj)
    file.path(CACHE_DIR, paste0(gsub("[^A-Za-z0-9_]+", "_", paste0("treat_flags_", proj)), ".rds"))
  dplyr::bind_rows(lapply(projects, function(proj) {
    f <- cfile(proj)
    if (file.exists(f)) { message("[cache] treat_flags_", proj, " (load)"); return(readRDS(f)) }
    cl <- tryCatch(TCGAbiolinks::GDCquery_clinic(project = proj, type = "clinical"),
                   error = function(e) e)
    if (inherits(cl, "error")) {
      message("  [!] clinical pull failed for ", proj, ": ", conditionMessage(cl),
              " — NOT cached, will retry next run")
      return(empty)
    }
    if (!any(c(PHARMA, RAD) %in% names(cl)))
      message("  [!] ", proj, ": no pharmaceutical or radiation field — all Unknown")
    df <- data.frame(
      Patient_ID = trimws(toupper(gsub("\\.", "-", cl$submitter_id))),
      Project_ID = proj,
      pharma_yes = flag(cl, PHARMA, "yes"), pharma_no = flag(cl, PHARMA, "no"),
      rad_yes    = flag(cl, RAD,    "yes"), rad_no    = flag(cl, RAD,    "no"),
      stringsAsFactors = FALSE) %>%
      dplyr::distinct(Patient_ID, .keep_all = TRUE)
    dir.create(CACHE_DIR, showWarnings = FALSE, recursive = TRUE)
    saveRDS(df, f)          # cache successful pulls only (incl. legit no-field cohorts)
    df
  }))
}

run_pancan_treated <- function(out_root, up_name = "Up", down_name = "Down") {
  module     <- "pancan_survival"
  in_dir     <- file.path(out_root, module)
  out_dir    <- file.path(in_dir, "treated")
  master_csv <- file.path(in_dir, "MASTER_Pan_Cancer_Clinical.csv")
  if (!file.exists(master_csv))
    stop("run_pancan_treated(): master not found: ", master_csv,
         " (run the pan-cancer survival module first).")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  message("\n======== PAN-CANCER SURVIVAL — TREATED (pharmaceutical or radiation) ========")
  message("  run: ", out_root)

  up_col   <- paste0(up_name,   SCORE_SUFFIX)
  down_col <- paste0(down_name, SCORE_SUFFIX)
  comp_col <- paste0("Composite", SCORE_SUFFIX)
  raw_cols <- c(up_col, down_col, comp_col)
  num <- function(x) suppressWarnings(as.numeric(x))

  # ---- reuse the existing scored master (NO RNA re-run) ----------------------
  master <- read.csv(master_csv, check.names = FALSE, stringsAsFactors = FALSE)
  master <- master %>%
    dplyr::mutate(dplyr::across(dplyr::all_of(c(raw_cols, "OS3Y_event", "OS3Y_delay",
                                                "RFS3Y_event", "RFS3Y_delay")), num))
  master$Patient_ID <- trimws(toupper(master$Patient_ID))
  # Restore the endpoint factor orientation derive_endpoints() set (read.csv drops
  # it), so the Wilcoxon/KM rows keep the module's event-first sign convention.
  master$Surv_3yr       <- factor(master$Surv_3yr,       levels = c("Dead_3yr", "Alive_3yr"))
  master$Recurrence_3yr <- factor(master$Recurrence_3yr, levels = c("Recurred", "Recurrence-Free"))

  # ---- treatment status (pharmaceutical OR radiation), cached per cohort -----
  # Same case_when ladder as the CRC module (crc_survival.R:169) so both sections
  # of the thesis use one treated definition.
  treat <- .pancan_treatment_flags(sort(unique(master$Project_ID)))
  master <- master %>%
    dplyr::left_join(dplyr::select(treat, Patient_ID, pharma_yes, pharma_no, rad_yes, rad_no),
                     by = "Patient_ID") %>%
    dplyr::mutate(Treatment_Status = dplyr::case_when(
      pharma_yes | rad_yes ~ "Treated",
      pharma_no  & rad_no  ~ "Not Treated",
      TRUE                 ~ "Unknown/Not Reported")) %>%
    dplyr::select(-pharma_yes, -pharma_no, -rad_yes, -rad_no)

  # ---- attrition report (auditable treated/total per cohort) -----------------
  counts <- master %>%
    dplyr::group_by(Project_ID) %>%
    dplyr::summarise(N_total     = dplyr::n(),
                     N_treated   = sum(Treatment_Status == "Treated"),
                     N_untreated = sum(Treatment_Status == "Not Treated"),
                     N_unknown   = sum(Treatment_Status == "Unknown/Not Reported"),
                     .groups = "drop") %>%
    dplyr::arrange(dplyr::desc(N_treated))
  write.csv(counts, file.path(out_dir, "treatment_status.csv"), row.names = FALSE)
  write.csv(master %>% dplyr::select(Patient_ID, Project_ID, Treatment_Status),
            file.path(out_dir, "treatment_status_perpatient.csv"), row.names = FALSE)
  message(sprintf("  Treated: %d / %d (%.0f%% pharma or radiation)  |  Not treated: %d  |  Unknown: %d",
                  sum(master$Treatment_Status == "Treated"), nrow(master),
                  100 * mean(master$Treatment_Status == "Treated"),
                  sum(master$Treatment_Status == "Not Treated"),
                  sum(master$Treatment_Status == "Unknown/Not Reported")))
  for (i in seq_len(nrow(counts)))
    message(sprintf("    %-14s treated %4d / %4d", counts$Project_ID[i],
                    counts$N_treated[i], counts$N_total[i]))

  mt_all <- master %>% dplyr::filter(Treatment_Status == "Treated")
  if (!nrow(mt_all)) stop("run_pancan_treated(): no treated patients found.")
  write.csv(mt_all, file.path(out_dir, "MASTER_Pan_Cancer_Clinical_Treated.csv"), row.names = FALSE)

  # ---- Phase 4a: per-cohort stats (uncorrected) on EVERY treated cohort ------
  # Feeds the per-cohort forest, the FDR_Stats CSV, and the aggregate forest's
  # cross-cohort direction tally. Treated subsets are small, so many cohorts will
  # fail the Cox gate and render as n/t — run_all_stats gates via MIN_COX_N /
  # MIN_EVENTS (Is_Testable); that is expected, not evidence of absence.
  rows <- list()
  for (proj in sort(unique(mt_all$Project_ID))) {
    cd <- mt_all %>% dplyr::filter(Project_ID == proj)
    for (sc in raw_cols) {
      rows[[length(rows) + 1]] <- run_all_stats(cd, sc, "OS",  "PanCancer", proj, sc, "Surv_3yr")
      rows[[length(rows) + 1]] <- run_all_stats(cd, sc, "RFS", "PanCancer", proj, sc, "Recurrence_3yr")
    }
  }

  # ---- Phase 4b: batch-corrected aggregate over well-populated treated cohorts -
  # removeBatchEffect needs populated batches; require >= MIN_COX_N treated per
  # cohort so a singleton batch cannot distort the pooled corrected scores.
  AGG_MIN <- MIN_COX_N
  keep    <- counts$Project_ID[counts$N_treated >= AGG_MIN]
  mt_agg  <- mt_all %>% dplyr::filter(Project_ID %in% keep)
  message(sprintf("  aggregate: %d of %d treated cohorts have >= %d treated patients",
                  length(keep), dplyr::n_distinct(mt_all$Project_ID), AGG_MIN))
  if (length(keep) >= 2) {
    sm        <- t(as.matrix(mt_agg[, raw_cols]))
    corrected <- limma::removeBatchEffect(sm, batch = mt_agg$Project_ID)
    corr_df   <- as.data.frame(t(corrected))
    corr_cols <- paste0(c(up_name, down_name, "Composite"), SCORE_SUFFIX, "_Corrected")
    colnames(corr_df) <- corr_cols
    mt_agg <- dplyr::bind_cols(mt_agg, corr_df)
    for (sc in corr_cols) {
      rows[[length(rows) + 1]] <- run_all_stats(mt_agg, sc, "OS",  "PanCancer", "Pan-Cancer", sc, "Surv_3yr")
      rows[[length(rows) + 1]] <- run_all_stats(mt_agg, sc, "RFS", "PanCancer", "Pan-Cancer", sc, "Recurrence_3yr")
    }
  } else {
    message("  [!] < 2 cohorts meet the treated floor — skipping batch-corrected aggregate.")
    mt_agg <- NULL
  }

  # ---- FDR (same Family x Test scheme as the all-patient module) -------------
  stats_df <- do.call(rbind, rows) %>%
    dplyr::mutate(Family = ifelse(Project == "Pan-Cancer", "PanCancer", "PerCohort"))
  stats_df <- apply_fdr(stats_df, by = c("Family", "Test"))
  write.csv(stats_df, file.path(out_dir, "FDR_Stats_Summary_Treated.csv"), row.names = FALSE)

  # ---- figures: the same set the all-patient section builds (reuse composites) -
  # Fig 6 A/B/C + composite, Fig 7 violins, and the aggregate appendix forest —
  # every builder reused unchanged, only the save names get a _Treated suffix.
  for (fn in c("build_pancan_forest", "build_pancan_aggregate",
               "build_pancan_violins", ".pancan_fig6_composite"))
    if (!exists(fn)) stop(fn, "() not found — source modules/composites.R first.")
  guard <- function(what, expr) tryCatch(expr, error = function(e) {
    msg <- paste0("[pancan treated] ", what, " FAILED: ", conditionMessage(e))
    message("  ", msg); warning(msg, call. = FALSE); NULL })

  # Per-cohort forests read Family == "PerCohort" and need the per-project score
  # SDs for the per-SD HR rescaling, so they take mt_all (every treated cohort),
  # not the batch-correction subset mt_agg. Each returns a make_fig closure.
  # order_score = up_col, not the default "Up_ssGSEA": up_name is caller-supplied, and
  # a mismatch would silently drop the shared ascending-Up-OS-HR row ordering.
  pA <- guard("Fig6A (Up forest)", build_pancan_forest(
    stats_df, mt_all, out_dir, up_col, "DTP Up (treated)", "Fig6A_pancan_forest_Up_Treated", up_col))
  pB <- guard("Fig6B (Composite forest)", build_pancan_forest(
    stats_df, mt_all, out_dir, comp_col, "DTP Composite (treated)", "Fig6B_pancan_forest_Composite_Treated", up_col))
  pC <- guard("Fig6C (Down forest)", build_pancan_forest(
    stats_df, mt_all, out_dir, down_col, "DTP Down (treated)", "Fig6C_pancan_forest_Down_Treated", up_col))
  guard("Fig6 composite",
        .pancan_fig6_composite(pA, pB, pC, out_dir, "Fig6_pancan_survival_composite_Treated"))

  # Fig 7 violins and the aggregate forest both re-fit removeBatchEffect on the
  # master they are handed, so both take mt_agg to stay consistent with the
  # Pan-Cancer stats rows they annotate themselves with.
  if (!is.null(mt_agg)) {
    guard("Fig7 violins", build_pancan_violins(
      stats_df, mt_agg, out_dir, save_name = "Fig7_pancan_violins_Treated"))
    guard("aggregate forest", build_pancan_aggregate(stats_df, mt_agg, out_dir))
  }

  message("Pan-cancer treated section complete -> ", out_dir)
  invisible(stats_df)
}

# Auto-run only when executed as a script (skipped when source()'d by main.R).
if (sys.nframe() == 0L) {
  suppressPackageStartupMessages({
    library(dplyr); library(ggplot2); library(patchwork); library(limma)
  })
  for (f in c("core/config.R", "core/io.R", "core/stats.R",
              "core/plotting.R", "core/gsea_plots.R")) source(f)
  source("modules/composites.R")
  .args <- commandArgs(trailingOnly = TRUE)
  root  <- if (length(.args) >= 1) .args[[1]] else {
    roots <- sort(Sys.glob("CRC_DTP_*"), decreasing = TRUE)
    roots <- roots[dir.exists(roots)]
    if (!length(roots)) stop("no CRC_DTP_* run root found; pass one as an argument.")
    roots[[1]]
  }
  run_pancan_treated(root)
}
