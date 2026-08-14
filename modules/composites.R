#!/usr/bin/env Rscript
# ==============================================================================
# modules/composites.R
# ------------------------------------------------------------------------------
# Standalone composite-figure module. Consolidates ALL composite-figure
# generation for the thesis into one file: the CRC-survival composites (5 groups,
# build_crc_composites), the pan-cancer composites (build_pancan_composites), and
# the metastasis composites (build_mets_gsea_composite + build_mets_liver_composite).
# PRESENTATION LAYER ONLY: reads a completed run's tabular / RDS outputs and draws
# the figures. It never re-runs an analysis module and never re-downloads.
#
# Entry point run_composites(out_root) auto-detects the newest CRC_DTP_* run root
# (or takes one) and regenerates every composite + publication figure for whichever
# analysis modules (crc_survival / pancan_survival / mets_de) produced output there.
# main.R sources this after the analysis modules and calls run_composites(out_root)
# once, after the three analyses and before finalize_run().
#
#   * in the pipeline:  sourced by main.R; run_composites(out_root) after the runs
#   * standalone:       Rscript modules/composites.R [OUT_ROOT]
#                       (OUT_ROOT defaults to the newest CRC_DTP_* run root)
#
# TITLED vs PUBLICATION copies. Every composite is saved twice: a titled copy to
# <run>/<module>/composites/ carrying the FULL statistics on every testable cell
# (native effect size + raw p + BH-FDR q, via the `detailed = TRUE` branch and the
# .stat_annot() helper in core/plotting.R), and a title/subtitle-stripped copy to
# <...>/composites/publication/ that keeps the clean asterisk-only convention
# (`detailed = FALSE`, byte-for-behaviour identical to the pre-consolidation output).
#
# Shared helpers come from core/plotting.R (sourced first): .base_theme, .sig_stars,
# .save_fig, .png_dev, .strip_titles, .fmt_p, .stat_annot, and the canonical
# DATASET_FULL / DATASET_SHORT / ENDPOINT_LABEL maps. The metastasis GSEA enrichment
# primitive plot_gsea_enrichment() stays in core/gsea_plots.R (it is a run-time
# helper used by modules/mets_de.R); build_mets_gsea_composite() calls it from here.
# ==============================================================================

suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr); library(readr)
  library(stringr); library(forcats); library(patchwork); library(scales)
  library(survival)
})

# ==============================================================================
# PART 1 — CRC-survival composites (build_crc_composites, 5 groups)
# ==============================================================================
# ---- Shared constants --------------------------------------------------------
# MIN_KM_GROUP_N / MIN_EVENTS live in core/config.R (sourced before this file in
# main.R and in the standalone entry at the bottom); referenced only inside
# function bodies, never at top-level parse time.
X_CAP          <- 36     # months (follow-up truncated/administratively censored at 36 mo)
FOREST_XLIM    <- c(0.1, 10)   # subgroup HR display window (log scale)

# Curated label + row-order dictionary for KNOWN signatures (groups 1 & 2).
# This is NOT a whitelist: signatures not listed here are still displayed by the
# CRC figures — .crc_modules() derives the display set from the score columns
# actually present in the run and labels any newcomer by its column name. Listing
# a signature here just gives it a nicer label and a fixed position; order below
# is the row order (inverse-signal modules DTP Down / MYC / Columnar grouped last).
# Keys MUST be the score columns the run actually writes ("<Signature>_ssGSEA",
# one per column of SIG_FILE); a key that matches nothing silently demotes that
# signature to the "novel" bucket (raw label, appended last).
MODULES <- c(
  Up_ssGSEA          = "DTP Up",
  Composite_ssGSEA   = "DTP Composite",
  Fetal_ssGSEA       = "Foetal intestinal",
  RSC_ssGSEA         = "Regenerative stem cell",
  revSC_ssGSEA       = "Revival revCSC",
  IBD_ssGSEA         = "IBD",
  Down_ssGSEA        = "DTP Down",
  MYC_ssGSEA         = "MYC module",
  CSC_ssGSEA         = "Columnar stem cell"
)

# Score columns present in a run ("<Sig>_ssGSEA") -> display labels + row order.
# Signatures listed in MODULES keep their curated label and order; any signature
# not in MODULES is still shown, labelled by its column name (minus _ssGSEA) and
# appended after the curated rows (deterministic, alphabetical). Returns:
#   cols   - score columns to iterate, in display order
#   labels - named vector score_col -> display label
#   levels - display labels in row order (factor levels)
.crc_modules <- function(score_cols) {
  score_cols <- unique(score_cols)
  known  <- intersect(names(MODULES), score_cols)   # curated, in MODULES order
  novel  <- sort(setdiff(score_cols, names(MODULES)))
  labels <- c(MODULES[known], setNames(sub("_ssGSEA$", "", novel), novel))
  list(cols = c(known, novel), labels = labels, levels = unname(labels))
}

# ---- Shared helpers ----------------------------------------------------------
# Bind dplyr verbs into a builder's local environment. In the full pipeline
# session the Bioconductor stack (AnnotationDbi/GSEABase via GSVA) masks
# dplyr::select() with an S4 generic, which breaks tbl_df pipelines; qualifying
# the verbs locally makes each builder immune to load-order masking.
.dplyr_local <- function(env) {
  for (fn in c("select","filter","mutate","transmute","summarise","group_by","arrange",
               "left_join","semi_join","distinct","pull","bind_rows","rename","count"))
    assign(fn, getExportedValue("dplyr", fn), envir = env)
}

# Null/NA-coalescing helper for optional-row column access (a Wilcoxon or KM row
# may be absent for a given cohort/score/metric).
`%||%` <- function(a, b)
  if (is.null(a) || length(a) == 0 || (length(a) == 1 && is.na(a))) b else a

.find_run_dir <- function() {
  cand <- sort(Sys.glob(file.path("CRC_DTP_*", "crc_survival")), decreasing = TRUE)
  if (!length(cand)) stop("No CRC_DTP_*/crc_survival directory found; pass RUN_DIR explicitly.")
  cand[[1]]
}
# Memoized CSV reader: repeated (.rd(run_dir, f)) calls in one session return a
# copy of the cached data.frame instead of re-reading from disk. Keyed by full
# path. Always return a row-slice copy so callers cannot mutate the cache entry.
.rd <- local({
  .cache <- new.env(parent = emptyenv())
  function(run_dir, f) {
    key <- file.path(run_dir, f)
    if (!exists(key, envir = .cache, inherits = FALSE))
      .cache[[key]] <- read.csv(key, check.names = FALSE, stringsAsFactors = FALSE)
    .cache[[key]][TRUE, , drop = FALSE]
  }
})
# Within-subset score SD for rescaling continuous Cox unit-change HRs to
# per-1-SD HRs via exp(log(HR) * sd). `x` is the score vector already restricted
# to the analysis rows (e.g. non-NA score / time / event inside a subgroup).
.score_sd <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  x <- x[!is.na(x)]
  if (length(x) < 2L) return(NA_real_)
  stats::sd(x)
}
# .sig_stars(), .save_fig(), .png_dev, .strip_titles(), .base_theme, and the
# canonical DATASET_FULL / DATASET_SHORT / ENDPOINT_LABEL maps live in
# core/plotting.R (sourced first) so this file, core/gsea_plots.R and
# core/pancan_composites.R share one definition.

# =============================================================================
# GROUP 1 — ssGSEA score vs 3-year outcome (Wilcoxon)
# =============================================================================
build_outcome_figures <- function(run_dir, out_dir, violin_modules = "all") {
  message("[group 1] outcome figures")
  .dplyr_local(environment())
  PAL <- c(Poor = "#E64B35", Good = "#4DBBD5")

  CONTRASTS <- tribble(
    ~key,          ~dataset,      ~cohort,             ~endpoint, ~stats_ds,    ~stats_cohort,        ~clin,   ~group_col,       ~poor,       ~good,
    "M_whole_OS",  "Marisa (GSE39582)", "Whole cohort",     "OS",  "GSE39582",  "GSE_All_Patients",  "gse",  "Surv_3yr",       "Dead_3yr", "Alive_3yr",
    "M_whole_RFS", "Marisa (GSE39582)", "Whole cohort",     "RFS", "GSE39582",  "GSE_All_Patients",  "gse",  "Recurrence_3yr", "Recurred", "Recurrence-Free",
    "M_trt_OS",    "Marisa (GSE39582)", "Adjuvant-treated", "OS",  "GSE39582",  "GSE_Treated",       "gse",  "Surv_3yr",       "Dead_3yr", "Alive_3yr",
    "M_trt_RFS",   "Marisa (GSE39582)", "Adjuvant-treated", "RFS", "GSE39582",  "GSE_Treated",       "gse",  "Recurrence_3yr", "Recurred", "Recurrence-Free",
    "T_whole_OS",  "TCGA-COAD",         "Whole cohort",     "OS",  "TCGA-COAD", "TCGA_All_Patients", "tcga", "Surv_3yr",       "Dead_3yr", "Alive_3yr",
    "T_whole_RFS", "TCGA-COAD",         "Whole cohort",     "RFS", "TCGA-COAD", "TCGA_All_Patients", "tcga", "Recurrence_3yr", "Recurred", "Recurrence-Free"
  )
  CONTRASTS$key     <- factor(CONTRASTS$key, levels = CONTRASTS$key)
  CONTRASTS$dataset <- factor(CONTRASTS$dataset, levels = unique(CONTRASTS$dataset))
  # Canonical endpoint labels (ENDPOINT_LABEL, core/plotting.R): "3-yr OS"/"3-yr RFS".
  CONTRASTS$col_lab <- paste0(CONTRASTS$cohort, "\n(", ENDPOINT_LABEL[CONTRASTS$endpoint], ")")
  col_lab_map <- setNames(CONTRASTS$col_lab, CONTRASTS$key)

  clin  <- list(gse = .rd(run_dir, "GSE39582_clinical.csv"),
                tcga = .rd(run_dir, "TCGA_COAD_clinical.csv"))
  stats <- .rd(run_dir, "CRC_Statistical_Summary.csv")

  # Signatures to display are derived from the score columns present in the run,
  # not a hardcoded list — a newly scored signature appears automatically.
  score_cols <- unique(unlist(lapply(clin, function(d)
                  grep("_ssGSEA$", names(d), value = TRUE))))
  MOD <- .crc_modules(score_cols)

  cohort_subset <- function(row) {
    d <- clin[[row$clin]]
    if (row$cohort == "Adjuvant-treated") d <- d[d$Chemo_adj == "Y" & !is.na(d$Chemo_adj), ]
    d
  }
  stat_lookup <- function(row, score) {
    s <- stats[stats$Dataset == row$stats_ds & stats$Project == row$stats_cohort &
               stats$Test == "Wilcoxon" & stats$Metric == row$endpoint &
               stats$Score == score, ]
    if (!nrow(s)) return(list(fdr = NA, r = NA, raw = NA, n = NA, sig = FALSE))
    list(fdr = s$FDR_P[1], r = s$Effect_r[1], raw = s$Raw_P[1], n = s$N[1],
         sig = isTRUE(as.logical(s$Is_Significant[1])))
  }

  long_rows <- list(); cell_rows <- list()
  for (i in seq_len(nrow(CONTRASTS))) {
    row <- CONTRASTS[i, ]
    d   <- cohort_subset(row)
    g   <- d[[row$group_col]]
    keep <- g %in% c(row$poor, row$good)
    d <- d[keep, ]; g <- g[keep]
    grp <- factor(ifelse(g == row$poor, "Poor", "Good"), levels = c("Poor", "Good"))
    for (sc in MOD$cols) {
      if (!sc %in% names(d)) next
      v  <- suppressWarnings(as.numeric(d[[sc]])); ok <- !is.na(v)
      long_rows[[length(long_rows) + 1]] <- tibble(
        key = row$key, module = MOD$labels[[sc]], score_col = sc, group = grp[ok], value = v[ok])
      st <- stat_lookup(row, sc)
      # signed rank-biserial r: pipeline sign is vs first level (poor), so flip
      # so + = higher in poor outcome for the red-poor palette.
      cell_rows[[length(cell_rows) + 1]] <- tibble(
        key = row$key, module = MOD$labels[[sc]], score_col = sc,
        fdr = st$fdr, raw = st$raw, signed_r = -st$r, sig = st$sig, stars = .sig_stars(st$fdr, "ns"))
    }
  }
  long  <- bind_rows(long_rows) %>%
    mutate(module = factor(module, MOD$levels), key = factor(key, levels(CONTRASTS$key)))
  cells <- bind_rows(cell_rows) %>%
    mutate(module = factor(module, MOD$levels), key = factor(key, levels(CONTRASTS$key))) %>%
    left_join(select(CONTRASTS, key, dataset), by = "key")

  # -- Fig 1A: effect-size matrix --
  # Two label variants per cell: clean (asterisk only, for publication) and full
  # (rank-biserial r + raw p + FDR q, for the titled report copy). make_matrix()
  # picks by `detailed`; publication copies stay identical to the clean output.
  cells$lab_clean <- ifelse(cells$sig, cells$stars, "")
  .eff_r <- ifelse(is.na(cells$signed_r), NA_character_, sprintf("r=%+.2f", cells$signed_r))
  cells$lab_full <- mapply(.stat_annot, detailed = TRUE, clean = cells$lab_clean,
                           effect = .eff_r, raw_p = cells$raw, fdr = cells$fdr)
  lim <- max(abs(cells$signed_r), na.rm = TRUE)
  make_matrix <- function(detailed) {
    d <- cells; d$lab <- if (detailed) d$lab_full else d$lab_clean
    ggplot(d, aes(x = key, y = fct_rev(module), fill = signed_r)) +
      geom_tile(colour = "white", linewidth = 0.6) +
      geom_text(aes(label = lab), size = if (detailed) 2.2 else 4.2, fontface = "bold",
                colour = "grey15", vjust = 0.5, lineheight = 0.82) +
      facet_grid(. ~ dataset, scales = "free_x", space = "free_x") +
      scale_fill_gradient2(low = "#3B6EA5", mid = "white", high = "#E64B35", midpoint = 0,
        limits = c(-lim, lim), name = "Rank-biserial r",
        guide = guide_colourbar(barwidth = 9, barheight = 0.5, title.vjust = 1)) +
      scale_x_discrete(labels = col_lab_map) +
      labs(title = "DTP-module scores vs 3-year outcome across cohorts",
           subtitle = if (detailed)
             "Colour = effect size (red: higher in poor-outcome group). Cells: rank-biserial r, raw p, FDR q (BH); * FDR<0.05."
           else
             "Colour = effect size (red: higher in poor-outcome group).  * marks FDR-significant Wilcoxon (BH, p<0.05).",
           x = NULL, y = NULL) +
      .base_theme + theme(axis.text.x = element_text(size = 8, lineheight = 0.85),
                          panel.spacing.x = unit(6, "pt"))
  }

  # -- Fig 1B: violin grid --
  if (identical(violin_modules, "all")) {
    curated_levels <- MOD$levels
  } else {
    keep_cols <- union(cells %>% filter(sig) %>% pull(score_col) %>% unique(),
                       c("Up_ssGSEA", "Composite_ssGSEA"))
    keep_lab  <- MOD$labels[intersect(keep_cols, names(MOD$labels))]
    curated_levels <- MOD$levels[MOD$levels %in% keep_lab]
  }
  vl <- long  %>% filter(module %in% curated_levels) %>% mutate(module = factor(module, curated_levels))
  vc <- cells %>% filter(module %in% curated_levels) %>% mutate(module = factor(module, curated_levels))
  yr <- vl %>% group_by(module, key) %>%
    summarise(ymin = min(value), ymax = max(value), .groups = "drop") %>% mutate(rng = ymax - ymin)
  vc <- vc %>% left_join(yr, by = c("module", "key"))
  # Clean vs full annotation for the violin panels (see Fig 1A note).
  vc$lab_clean_v <- ifelse(vc$stars == "", "ns", vc$stars)
  .eff_rv <- ifelse(is.na(vc$signed_r), NA_character_, sprintf("r=%+.2f", vc$signed_r))
  vc$lab_full_v <- mapply(.stat_annot, detailed = TRUE, clean = vc$lab_clean_v,
                          effect = .eff_rv, raw_p = vc$raw, fdr = vc$fdr)
  nlab <- vl %>% group_by(module, key, group) %>% summarise(n = n(), .groups = "drop") %>%
    left_join(select(yr, module, key, ymin, rng), by = c("module", "key"))

  build_violin_block <- function(keys, title, show_y, detailed) {
    dl <- filter(vl, key %in% keys); dc <- filter(vc, key %in% keys); dn <- filter(nlab, key %in% keys)
    dc$lab <- if (detailed) dc$lab_full_v else dc$lab_clean_v
    ggplot(dl, aes(x = group, y = value, fill = group)) +
      geom_violin(alpha = 0.85, trim = FALSE, linewidth = 0.3, colour = "grey30") +
      geom_boxplot(width = 0.16, fill = "white", alpha = 0.9, outlier.shape = NA, linewidth = 0.3) +
      geom_text(data = dc, inherit.aes = FALSE,
                aes(x = 1.5, y = ymax + 0.12 * rng, label = lab),
                fontface = "bold", size = if (detailed) 2.1 else 3.2, lineheight = 0.82) +
      geom_text(data = dn, inherit.aes = FALSE,
                aes(x = group, y = ymin - 0.10 * rng, label = paste0("n=", n)),
                size = 2.4, colour = "grey35") +
      facet_grid(module ~ key, scales = "free_y", switch = "y",
                 labeller = labeller(key = col_lab_map, module = label_wrap_gen(14))) +
      scale_fill_manual(values = PAL, name = "3-year outcome",
        labels = c(Poor = "Deceased / recurred", Good = "Alive / recurrence-free")) +
      scale_x_discrete(labels = c(Poor = "Poor", Good = "Good")) +
      scale_y_continuous(expand = expansion(mult = 0.16)) +
      labs(title = title, x = NULL, y = if (show_y) "ssGSEA score" else NULL) +
      .base_theme +
      theme(plot.title = element_text(face = "bold", size = 11, hjust = 0.5),
            strip.placement = "outside",
            strip.text.y.left = if (show_y) element_text(angle = 0, face = "bold", size = 8.5) else element_blank(),
            strip.text.x = element_text(size = 8, lineheight = 0.85),
            axis.text.x = element_text(size = 8), panel.spacing = unit(4, "pt"))
  }
  make_violin <- function(detailed) {
    block_m <- build_violin_block(c("M_whole_OS","M_whole_RFS","M_trt_OS","M_trt_RFS"), "Marisa (GSE39582)", TRUE, detailed)
    block_t <- build_violin_block(c("T_whole_OS","T_whole_RFS"), "TCGA-COAD", FALSE, detailed)
    block_m + block_t +
      plot_layout(widths = c(4, 2), guides = "collect") +
      plot_annotation(
        title = if (identical(violin_modules, "all")) "Signature-score distributions by 3-year outcome"
                else "Score distributions for outcome-associated modules",
        subtitle = if (detailed)
          "Full cohort, true n shown; y-axis free within each dataset. Panels: rank-biserial r, raw p, FDR q (BH); * FDR<0.05, ns = not."
        else
          "Full cohort, true n shown; y-axis free within each dataset (ssGSEA scale differs GSE vs TCGA).  * = FDR-significant; ns = not.",
        theme = theme(plot.title = element_text(face = "bold", size = 14, hjust = 0),
                      plot.subtitle = element_text(size = 10, colour = "grey30", hjust = 0))) &
      theme(legend.position = "bottom")
  }

  n_rows <- length(curated_levels); viol_h <- 1.5 * n_rows + 1.4
  # Titled report copies carry the full statistics (detailed = TRUE); the
  # publication/ copies are the clean, title-stripped versions (detailed = FALSE).
  .save_fig(make_matrix(TRUE), "Fig1A_outcome_effect_matrix", 9.5, 4.6, out_dir)
  .save_fig(make_violin(TRUE), "Fig1B_outcome_violins", 11.0, viol_h, out_dir)
  mk <- function(m, v) wrap_elements(m) / wrap_elements(v) +
    plot_layout(heights = c(4.6, viol_h - 0.4)) + plot_annotation(tag_levels = "A") &
    theme(plot.tag = element_text(face = "bold", size = 16))
  H <- 4.6 + viol_h + 0.6
  .save_fig(mk(make_matrix(TRUE), make_violin(TRUE)), "Fig1_outcome_composite", 11.5, H, out_dir)
  .save_fig(mk(.strip_titles(make_matrix(FALSE)), .strip_titles(make_violin(FALSE))),
            "Fig1_outcome_composite", 11.5, H, file.path(out_dir, "publication"))
}

# =============================================================================
# GROUP 2 — Kaplan-Meier by median High/Low score
# =============================================================================
build_km_figures <- function(run_dir, out_dir, km_modules = "all") {
  message("[group 2] KM figures")
  .dplyr_local(environment())
  PAL <- c(Low = "#4DBBD5", High = "#E64B35")

  CONTRASTS <- tribble(
    ~key,          ~dataset,             ~cohort,             ~metric, ~stats_ds,   ~stats_cohort,       ~clin,
    "M_whole_OS",  "Marisa (GSE39582)", "Whole cohort",      "OS",   "GSE39582",  "GSE_All_Patients",  "gse",
    "M_whole_RFS", "Marisa (GSE39582)", "Whole cohort",      "RFS",  "GSE39582",  "GSE_All_Patients",  "gse",
    "M_trt_OS",    "Marisa (GSE39582)", "Adjuvant-treated",  "OS",   "GSE39582",  "GSE_Treated",       "gse",
    "M_trt_RFS",   "Marisa (GSE39582)", "Adjuvant-treated",  "RFS",  "GSE39582",  "GSE_Treated",       "gse",
    "T_whole_OS",  "TCGA-COAD",         "Whole cohort",      "OS",   "TCGA-COAD", "TCGA_All_Patients", "tcga",
    "T_whole_RFS", "TCGA-COAD",         "Whole cohort",      "RFS",  "TCGA-COAD", "TCGA_All_Patients", "tcga"
  )
  CONTRASTS$key     <- factor(CONTRASTS$key, levels = CONTRASTS$key)
  CONTRASTS$dataset <- factor(CONTRASTS$dataset, levels = unique(CONTRASTS$dataset))
  # Canonical endpoint labels (ENDPOINT_LABEL, core/plotting.R): "3-yr OS"/"3-yr RFS".
  CONTRASTS$col_lab <- paste0(CONTRASTS$cohort, "\n(", ENDPOINT_LABEL[CONTRASTS$metric], ")")
  col_lab_map <- setNames(CONTRASTS$col_lab, CONTRASTS$key)

  # metric_cols() is the canonical helper in core/stats.R (sourced before this file).
  clin  <- list(gse = .rd(run_dir, "GSE39582_clinical.csv"),
                tcga = .rd(run_dir, "TCGA_COAD_clinical.csv"))
  stats <- .rd(run_dir, "CRC_Statistical_Summary.csv")

  # Display set derived from the score columns present (see .crc_modules) — a
  # newly scored signature appears automatically, no hardcoded list.
  score_cols <- unique(unlist(lapply(clin, function(d)
                  grep("_ssGSEA$", names(d), value = TRUE))))
  MOD <- .crc_modules(score_cols)

  cohort_subset <- function(row) {
    d <- clin[[row$clin]]
    if (row$cohort == "Adjuvant-treated") d <- d[d$Chemo_adj == "Y" & !is.na(d$Chemo_adj), ]
    d
  }
  cox_stat <- function(row, score) {   # continuous-score Cox: inferential effect
    s <- stats[stats$Dataset == row$stats_ds & stats$Project == row$stats_cohort &
               stats$Test == "Cox" & stats$Metric == row$metric & stats$Score == score, ]
    if (!nrow(s) || !isTRUE(as.logical(s$Is_Testable[1])))
      return(list(hr_unit = NA_real_, fdr = NA_real_, raw = NA_real_, sig = FALSE, testable = FALSE))
    list(hr_unit = as.numeric(s$HR[1]), fdr = as.numeric(s$FDR_P[1]),
         raw = as.numeric(s$Raw_P[1]), sig = isTRUE(as.logical(s$Is_Significant[1])), testable = TRUE)
  }
  logrank_stat <- function(row, score) {   # log-rank: shown on the curves
    s <- stats[stats$Dataset == row$stats_ds & stats$Project == row$stats_cohort &
               stats$Test == "KM" & stats$Metric == row$metric & stats$Score == score, ]
    if (!nrow(s)) return(list(fdr = NA_real_, raw = NA_real_, sig = FALSE))
    list(fdr = as.numeric(s$FDR_P[1]), raw = as.numeric(s$Raw_P[1]),
         sig = isTRUE(as.logical(s$Is_Significant[1])))
  }

  curve_rows <- list(); cell_rows <- list()
  for (i in seq_len(nrow(CONTRASTS))) {
    row <- CONTRASTS[i, ]; d0 <- cohort_subset(row)
    mc <- metric_cols(row$metric); tcol <- mc[["t"]]; ecol <- mc[["e"]]
    for (sc in MOD$cols) {
      if (!all(c(sc, tcol, ecol) %in% names(d0))) next
      d <- d0[!is.na(d0[[sc]]) & !is.na(d0[[tcol]]) & !is.na(d0[[ecol]]), ]
      score <- as.numeric(d[[sc]])
      d <- data.frame(time = pmin(as.numeric(d[[tcol]]), X_CAP), event = as.numeric(d[[ecol]]),
                      Group = factor(ifelse(score > median(score), "High", "Low"), levels = c("Low", "High")))
      tb <- table(d$Group)
      drawable <- length(tb) == 2 && all(tb >= MIN_KM_GROUP_N) && sum(d$event) >= MIN_EVENTS
      cx <- cox_stat(row, sc); lr <- logrank_stat(row, sc)
      hr_sd <- if (cx$testable && is.finite(cx$hr_unit)) exp(log(cx$hr_unit) * .score_sd(score)) else NA_real_
      cell_rows[[length(cell_rows) + 1]] <- tibble(
        key = row$key, module = MOD$labels[[sc]], score_col = sc, testable = cx$testable,
        hr = hr_sd, log2hr = log2(hr_sd), fdr = cx$fdr, raw = cx$raw, sig = cx$sig,
        stars = .sig_stars(cx$fdr, "ns"), lr_fdr = lr$fdr, lr_raw = lr$raw, lr_sig = lr$sig)
      if (!drawable) next
      sf <- survfit(Surv(time, event) ~ Group, data = d)
      strata <- rep(sub(".*=", "", names(sf$strata)), sf$strata)
      cd <- rbind(data.frame(time = 0, surv = 1, ncens = 0, Group = levels(d$Group)),
                  data.frame(time = sf$time, surv = sf$surv, ncens = sf$n.censor, Group = strata))
      cd$key <- row$key; cd$module <- MOD$labels[[sc]]; cd$small <- nrow(d) < 500
      curve_rows[[length(curve_rows) + 1]] <- cd
    }
  }
  cells <- bind_rows(cell_rows) %>%
    mutate(module = factor(module, MOD$levels), key = factor(key, levels(CONTRASTS$key))) %>%
    left_join(select(CONTRASTS, key, dataset), by = "key")
  curves <- bind_rows(curve_rows) %>%
    mutate(module = factor(module, MOD$levels), key = factor(key, levels(CONTRASTS$key)),
           Group = factor(Group, levels = c("Low", "High")))

  # -- Fig 2A: Cox HR/SD matrix --
  # Clean label = asterisk only; full label = per-SD HR + raw p + FDR q.
  cells$lab_clean <- ifelse(cells$sig, cells$stars, "")
  .eff_hr <- ifelse(is.na(cells$hr), NA_character_, sprintf("HR=%.2f", cells$hr))
  cells$lab_full <- mapply(.stat_annot, detailed = TRUE, clean = cells$lab_clean,
                           effect = .eff_hr, raw_p = cells$raw, fdr = cells$fdr)
  lim <- min(max(abs(cells$log2hr[is.finite(cells$log2hr)]), na.rm = TRUE), 2)
  make_matrix <- function(detailed) {
    d <- cells; d$lab <- if (detailed) d$lab_full else d$lab_clean
    ggplot(d, aes(x = key, y = fct_rev(module), fill = pmax(pmin(log2hr, lim), -lim))) +
      geom_tile(colour = "white", linewidth = 0.6) +
      geom_text(aes(label = lab), size = if (detailed) 2.2 else 4.2, fontface = "bold",
                colour = "grey15", vjust = 0.5, lineheight = 0.82) +
      facet_grid(. ~ dataset, scales = "free_x", space = "free_x") +
      scale_fill_gradient2(low = "#3B6EA5", mid = "white", high = "#E64B35", midpoint = 0,
        limits = c(-lim, lim), name = expression(log[2]~"HR (per 1 SD of score)"),
        guide = guide_colourbar(barwidth = 9, barheight = 0.5, title.vjust = 1), na.value = "grey85") +
      scale_x_discrete(labels = col_lab_map) +
      labs(title = "Survival association of DTP-module scores across cohorts",
           subtitle = if (detailed)
             "Univariable continuous-score Cox PH, follow-up truncated at 36 mo. Cells: per-SD HR, raw p, FDR q (BH); * FDR<0.05."
           else
             "Univariable continuous-score Cox PH, follow-up truncated at 36 mo.  Colour = log2 HR per 1 SD (red = higher score worse).  * = FDR<0.05 (BH).",
           x = NULL, y = NULL) +
      .base_theme + theme(axis.text.x = element_text(size = 8, lineheight = 0.85),
                          panel.spacing.x = unit(6, "pt"))
  }

  # -- Fig 2B: KM curve grid --
  km_levels <- if (identical(km_modules, "all")) MOD$levels else {
    keep_cols <- union(cells %>% filter(sig) %>% pull(score_col) %>% unique(),
                       c("Up_ssGSEA", "Composite_ssGSEA"))
    keep_lab  <- MOD$labels[intersect(keep_cols, names(MOD$labels))]
    MOD$levels[MOD$levels %in% keep_lab]
  }
  cv <- curves %>% filter(module %in% km_levels) %>% mutate(module = factor(module, km_levels))
  drawn <- distinct(cv, key, module)
  cc <- cells %>% filter(module %in% km_levels) %>% semi_join(drawn, by = c("key", "module")) %>%
    mutate(module = factor(module, km_levels),
           # Clean = "log-rank <stars>"; full adds the raw log-rank p + FDR q.
           lab_clean = ifelse(is.na(lr_fdr), "log-rank n/a",
                              paste0("log-rank ", .sig_stars(lr_fdr, "n.s."))),
           lab_full = mapply(.stat_annot, detailed = TRUE, clean = lab_clean,
                             effect = "log-rank", raw_p = lr_raw, fdr = lr_fdr))

  build_km_block <- function(keys, title, show_y, detailed) {
    dc <- filter(cv, key %in% keys); da <- filter(cc, key %in% keys)
    da$lab <- if (detailed) da$lab_full else da$lab_clean
    ggplot(dc, aes(x = time, y = surv, colour = Group)) +
      geom_step(linewidth = 0.6) +
      geom_point(data = subset(dc, ncens > 0 & small), shape = 3, size = 1.1, show.legend = FALSE) +
      geom_text(data = da, inherit.aes = FALSE, aes(x = 0.5, y = 0.06, label = lab),
                hjust = 0, vjust = 0, size = 2.5, colour = "grey25", lineheight = 0.9) +
      facet_grid(module ~ key, switch = "y",
                 labeller = labeller(key = col_lab_map, module = label_wrap_gen(14))) +
      scale_colour_manual(values = PAL, name = "ssGSEA score",
        labels = c(Low = "Low (< median)", High = "High (> median)")) +
      scale_x_continuous(limits = c(0, X_CAP), breaks = c(0, 12, 24, 36), expand = expansion(mult = c(0, 0.02))) +
      scale_y_continuous(limits = c(0, 1), labels = percent_format(accuracy = 1), expand = expansion(mult = c(0, 0.02))) +
      labs(title = title, x = "Months from diagnosis (follow-up truncated at 36)", y = if (show_y) "Survival probability" else NULL) +
      .base_theme +
      theme(plot.title = element_text(face = "bold", size = 11, hjust = 0.5),
            strip.placement = "outside",
            strip.text.y.left = if (show_y) element_text(angle = 0, face = "bold", size = 8.5) else element_blank(),
            strip.text.x = element_text(size = 8, lineheight = 0.85),
            panel.spacing.x = unit(11, "pt"), panel.spacing.y = unit(5, "pt"),
            panel.grid.major.y = element_line(colour = "grey93", linewidth = 0.3))
  }
  make_km <- function(detailed) {
    block_m <- build_km_block(c("M_whole_OS","M_whole_RFS","M_trt_OS","M_trt_RFS"), "Marisa (GSE39582)", TRUE, detailed)
    block_t <- build_km_block(c("T_whole_OS","T_whole_RFS"), "TCGA-COAD", FALSE, detailed)
    block_m + block_t +
      plot_layout(widths = c(4, 2), guides = "collect") +
      plot_annotation(
        title = if (identical(km_modules, "all")) "Kaplan-Meier survival by signature score (median High/Low split)"
                else "Kaplan-Meier survival by score, outcome-associated modules",
        subtitle = if (detailed)
          "Median High/Low split; curves annotated with log-rank raw p and BH-FDR q. Continuous-Cox per-SD HR effect sizes are in the matrix panel. Time runs from diagnosis with follow-up truncated at 36 mo (NOT a 3-yr landmark)."
        else
          "Median High/Low split; log-rank BH-FDR significance shown as stars (* <0.05, ** <0.01, *** <0.001; n.s. = not significant). Continuous-Cox HR effect sizes are in the matrix panel. Time runs from diagnosis with follow-up truncated at 36 mo (NOT a 3-yr landmark).",
        theme = theme(plot.title = element_text(face = "bold", size = 14, hjust = 0),
                      plot.subtitle = element_text(size = 10, colour = "grey30", hjust = 0))) &
      theme(legend.position = "bottom")
  }

  n_rows <- length(km_levels); km_h <- 2.0 * n_rows + 1.4
  .save_fig(make_matrix(TRUE), "Fig2A_km_hr_matrix", 9.5, 4.6, out_dir)
  .save_fig(make_km(TRUE),     "Fig2B_km_curves",    11.0, km_h, out_dir)
  mk <- function(m, k) wrap_elements(m) / wrap_elements(k) +
    plot_layout(heights = c(4.6, km_h - 0.4)) + plot_annotation(tag_levels = "A") &
    theme(plot.tag = element_text(face = "bold", size = 16))
  H <- 4.6 + km_h + 0.6
  .save_fig(mk(make_matrix(TRUE), make_km(TRUE)), "Fig2_km_composite", 11.5, H, out_dir)
  .save_fig(mk(.strip_titles(make_matrix(FALSE)), .strip_titles(make_km(FALSE))),
            "Fig2_km_composite", 11.5, H, file.path(out_dir, "publication"))
}

# =============================================================================
# GROUP 3 — effect modification (subgroup forests)
# =============================================================================
build_subgroup_figures <- function(run_dir, out_dir, primary_score = "Up_ssGSEA") {
  message("[group 3] subgroup forests")
  .dplyr_local(environment())
  SCORES <- c(Up_ssGSEA = "DTP Up", Composite_ssGSEA = "DTP Composite", Down_ssGSEA = "DTP Down")
  MOD_LABEL  <- c(CMS = "CMS subtype", PDS = "PDS subtype", Stage_bin = "Stage", MSI_group = "MSI status")
  MOD_LEVELS <- names(MOD_LABEL)
  LEVEL_ORDER <- c("CMS1","CMS2","CMS3","CMS4","PDS1","PDS2","PDS3","Early","Late","MSS","MSI")
  COLS <- tibble(ckey = c("GSE_OS","GSE_RFS","GSEt_OS","GSEt_RFS","TCGA_OS","TCGA_RFS"),
                 Dataset = c("GSE39582","GSE39582","GSE39582 (treated)","GSE39582 (treated)","TCGA-COAD","TCGA-COAD"),
                 Metric  = c("OS","RFS","OS","RFS","OS","RFS"))
  # Canonical short cohort + endpoint labels (DATASET_SHORT / ENDPOINT_LABEL, core/plotting.R).
  COLS$clab <- paste0(DATASET_SHORT[COLS$Dataset], "\n(", ENDPOINT_LABEL[COLS$Metric], ")")
  col_lab_map <- setNames(COLS$clab, COLS$ckey)
  DS_PREFIX <- c("GSE39582" = "GSE", "GSE39582 (treated)" = "GSEt", "TCGA-COAD" = "TCGA")
  DIR_PAL <- c("Higher hazard (HR>1)" = "#E64B35", "Lower hazard (HR<1)" = "#3B6EA5")

  subg  <- .rd(run_dir, "Subgroup_Score_HRs.csv")
  inter <- .rd(run_dir, "Interaction_Cox_Summary.csv")
  clin  <- list(GSE39582 = .rd(run_dir, "GSE39582_clinical.csv"),
                "TCGA-COAD" = .rd(run_dir, "TCGA_COAD_clinical.csv"))
  # Harmonised modifier columns (CMS / PDS / Stage_bin / MSI_group) rebuilt with the
  # SAME rule the analysis module used, so sd_of() can restrict to the subgroup Level.
  clin$GSE39582       <- harmonize_crc_modifiers(clin$GSE39582, "GSE39582")
  clin[["TCGA-COAD"]] <- harmonize_crc_modifiers(clin[["TCGA-COAD"]], "TCGA-COAD")
  # Treated-Marisa subset so sd_of() rescales the treated columns on the treated SD.
  clin[["GSE39582 (treated)"]] <- clin$GSE39582[which(clin$GSE39582$Chemo_adj == "Y"), , drop = FALSE]
  metric_evt <- c(OS = "OS3Y_event", RFS = "RFS3Y_event")
  metric_tim <- c(OS = "OS3Y_delay", RFS = "RFS3Y_delay")
  # WITHIN-SUBGROUP score SD, over the same rows the subgroup Cox modelled (non-NA
  # score + time + event inside that modifier Level). The HR being rescaled is the
  # score effect INSIDE the level, so 1 SD must be that level's SD — using the
  # cohort-wide SD would print a different "HR per 1 SD" than the §3.3 figure/table
  # (build_subtype_survival_figures) report for the identical Cox fit.
  # Core SD of the filtered score vector is .score_sd() (file-local helper above).
  sd_of <- function(ds, metric, score, modifier, level) {
    d <- clin[[ds]]; ev <- metric_evt[[metric]]; tm <- metric_tim[[metric]]
    if (!all(c(score, ev, tm, modifier) %in% names(d))) return(NA_real_)
    d <- d[which(as.character(d[[modifier]]) == level), , drop = FALSE]
    v <- suppressWarnings(as.numeric(d[[score]]))
    ok <- !is.na(v) & !is.na(d[[tm]]) & !is.na(d[[ev]])
    .score_sd(v[ok])
  }

  subg <- subg %>%
    mutate(ckey = paste0(unname(DS_PREFIX[Dataset]), "_", Metric),
           # Modifier/Level are still character here (factored further down).
           sd_sc = mapply(sd_of, Dataset, Metric, Score, Modifier, Level),
           hr_sd = exp(log(HR) * sd_sc), lo_sd = exp(log(HR_lower) * sd_sc), hi_sd = exp(log(HR_upper) * sd_sc),
           Modifier = factor(Modifier, MOD_LEVELS), Level = factor(Level, LEVEL_ORDER),
           ckey = factor(ckey, COLS$ckey),
           # Non-convergent / failed Cox: non-finite per-SD point or CI bounds.
           # (Pure NA = missing/not modelled — not flagged here.) Do not squish these
           # to FOREST_XLIM boundaries as if they were real off-scale effects.
           nonconv = (!is.na(hr_sd) & !is.finite(hr_sd)) |
                     (!is.na(lo_sd) & !is.finite(lo_sd)) |
                     (!is.na(hi_sd) & !is.finite(hi_sd)),
           offscale = !nonconv & !is.na(hr_sd) &
                      (hr_sd < FOREST_XLIM[1] | hr_sd > FOREST_XLIM[2]),
           hr_d = ifelse(nonconv, NA_real_, squish(hr_sd, FOREST_XLIM)),
           lo_d = ifelse(nonconv, NA_real_, squish(lo_sd, FOREST_XLIM)),
           hi_d = ifelse(nonconv, NA_real_, squish(hi_sd, FOREST_XLIM)),
           # Per-row left-edge label: n / e (clean) plus the per-SD HR + CI (full).
           ne_clean = sprintf("n=%d, e=%d", as.integer(N), as.integer(N_events)),
           ne_full  = ifelse(is.na(hr_sd) | nonconv, ne_clean,
                             sprintf("%s\nHR=%.2f [%.2f-%.2f]", ne_clean, hr_sd, lo_sd, hi_sd)))
  if (any(subg$nonconv))
    message(sprintf(paste0("[group 3] excluding %d non-convergent Cox row(s) from ",
                           "forest display (non-finite per-SD HR or CI; not squished ",
                           "to axis boundary)"), sum(subg$nonconv)))
  inter <- inter %>%
    mutate(ckey = paste0(unname(DS_PREFIX[Dataset]), "_", Metric),
           Modifier = factor(Modifier, MOD_LEVELS), ckey = factor(ckey, COLS$ckey),
           # Interaction LRT: clean = "interaction <stars>"; full adds raw p + FDR q
           # on ONE line (kept single-line so it never overruns the top data row).
           int_lab = ifelse(is.na(FDR_P), "interaction n/a",
                            paste0("interaction ", .sig_stars(FDR_P, "n.s."))),
           int_full = mapply(function(rp, fp) {
                        if (is.na(fp)) return("interaction n/a")
                        paste0("interaction: ",
                               paste(stats::na.omit(c(.fmt_p("p", rp),
                                     paste0(.fmt_p("q", fp), .sig_stars(fp, "")))), collapse = " "))
                      }, Raw_P, FDR_P),
           int_sig = !is.na(Is_Significant) & as.logical(Is_Significant))

  # -- Fig 3A: interaction matrix --
  mat <- inter %>% filter(Score %in% names(SCORES)) %>%
    mutate(Score = factor(SCORES[Score], unname(SCORES)),
           Modifier = fct_rev(factor(MOD_LABEL[as.character(Modifier)], unname(MOD_LABEL))),
           mlog = pmin(-log10(FDR_P), 4),
           lab_clean = .sig_stars(FDR_P),
           lab_full  = mapply(.stat_annot, detailed = TRUE, clean = lab_clean,
                              raw_p = Raw_P, fdr = FDR_P))
  make_matrix <- function(detailed) {
    d <- mat; d$lab <- if (detailed) d$lab_full else d$lab_clean
    ggplot(d, aes(x = ckey, y = Modifier, fill = mlog)) +
      geom_tile(colour = "white", linewidth = 0.6) +
      geom_text(aes(label = lab), fontface = "bold", size = if (detailed) 2.4 else 4.2,
                vjust = 0.5, lineheight = 0.82) +
      facet_grid(. ~ Score) +
      scale_fill_gradient(low = "white", high = "#762A83", limits = c(0, 4),
        name = expression(-log[10]~"interaction FDR"),
        guide = guide_colourbar(barwidth = 9, barheight = 0.5, title.vjust = 1)) +
      scale_x_discrete(labels = col_lab_map) +
      labs(title = "Effect modification of the DTP signature (score-by-subgroup interaction)",
           subtitle = if (detailed)
             "LRT of score*modifier vs score+modifier (Cox, follow-up truncated at 36 mo). Cells: raw p, FDR q (BH); * FDR<0.05. Colour = -log10 FDR."
           else
             "Likelihood-ratio test of score*modifier vs score+modifier (Cox, follow-up truncated at 36 mo).  * marks FDR-significant (BH, p<0.05).",
           x = NULL, y = NULL) +
      .base_theme + theme(axis.text.x = element_text(size = 8, lineheight = 0.85),
                          panel.spacing.x = unit(6, "pt"))
  }

  # -- Fig 3B: per-score forests --
  build_forest <- function(score_col, detailed) {
    # Drop non-convergent rows so they are not drawn at a fabricated boundary.
    ds <- subg  %>% filter(Score == score_col, !nonconv)
    di <- inter %>% filter(Score == score_col)
    ds$dir     <- ifelse(ds$hr_sd > 1, "Higher hazard (HR>1)", "Lower hazard (HR<1)")
    ds$ne_lab  <- if (detailed) ds$ne_full else ds$ne_clean   # +per-row HR/CI when detailed
    di$int_txt <- if (detailed) di$int_full else di$int_lab   # +interaction raw p / FDR q
    ggplot(ds, aes(x = hr_d, y = fct_rev(Level))) +
      geom_vline(xintercept = 1, linetype = "dashed", colour = "grey55") +
      geom_errorbarh(aes(xmin = lo_d, xmax = hi_d, colour = dir), height = 0.28, linewidth = 0.5) +
      geom_point(aes(colour = dir, shape = offscale), size = 2.1) +
      # Per-cohort n / events (+ per-SD HR & CI when detailed), drawn INSIDE each
      # facet (cohort x endpoint) at the left edge — the counts differ by cohort and
      # by OS/RFS, so they cannot live on the single shared y-axis label.
      geom_text(inherit.aes = FALSE,
                aes(x = FOREST_XLIM[1], y = fct_rev(Level), label = ne_lab),
                hjust = 0, vjust = -0.5, size = 1.95, colour = "grey45", lineheight = 0.82) +
      geom_text(data = di, inherit.aes = FALSE, aes(x = FOREST_XLIM[1], y = Inf, label = int_txt),
                hjust = 0, vjust = 1.25, size = 2.5, colour = ifelse(di$int_sig, "#B2182B", "grey35"),
                lineheight = 0.82) +
      scale_y_discrete() +
      facet_grid(Modifier ~ ckey, scales = "free_y", space = "free_y", switch = "y",
                 labeller = labeller(ckey = col_lab_map, Modifier = function(x) MOD_LABEL[x])) +
      scale_x_log10(limits = FOREST_XLIM, breaks = c(0.1, 0.3, 1, 3, 10),
                    labels = c("0.1","0.3","1","3","10")) +
      scale_colour_manual(values = DIR_PAL, name = "Score effect", drop = FALSE) +
      scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 21),
                         labels = c(`FALSE` = "in range", `TRUE` = "off-scale (squished)"), name = NULL) +
      labs(title = paste0("Per-subgroup hazard ratio of ", SCORES[[score_col]], " score"),
           subtitle = if (detailed)
             "Continuous-score Cox HR per 1 SD within each subgroup (95% CI); x truncated to [0.1, 10]. Rows: n/e + per-SD HR [CI]; interaction raw p / FDR q per facet."
           else
             "Continuous-score Cox HR per 1 SD within each subgroup (95% CI); x truncated to [0.1, 10]. n / e = subgroup size / events in that cohort.",
           x = "Hazard ratio per 1 SD (log scale)", y = NULL) +
      .base_theme +
      theme(strip.placement = "outside",
            strip.text.y.left = element_text(angle = 0, face = "bold", size = 8.5),
            strip.text.x = element_text(size = 8, lineheight = 0.85),
            axis.text.y = element_text(size = 7.5), panel.spacing = unit(4, "pt"),
            panel.grid.major.x = element_line(colour = "grey93", linewidth = 0.3))
  }
  .save_fig(make_matrix(TRUE), "Fig3A_interaction_matrix", 12.5, 3.4, out_dir)
  for (sc in names(SCORES))
    .save_fig(build_forest(sc, TRUE), paste0("Fig3B_subgroup_forest_", sub("_ssGSEA", "", sc)), 14.0, 8.0, out_dir)
  mk <- function(m, f) wrap_elements(m) / wrap_elements(f) +
    plot_layout(heights = c(3.4, 8.2)) + plot_annotation(tag_levels = "A") &
    theme(plot.tag = element_text(face = "bold", size = 16))
  .save_fig(mk(make_matrix(TRUE), build_forest(primary_score, TRUE)),
            "Fig3_subgroup_composite", 14.0, 12.2, out_dir)
  .save_fig(mk(.strip_titles(make_matrix(FALSE)), .strip_titles(build_forest(primary_score, FALSE))),
            "Fig3_subgroup_composite", 14.0, 12.2, file.path(out_dir, "publication"))
}

# =============================================================================
# GROUP 4 — univariable DTP-to-outcome WITHIN molecular/clinical subgroups
#           (thesis Section 3.3). Distinct from group 3 (§3.4 multivariable
#           effect-modification forest) and from groups 1-2 (§3.2 whole-cohort
#           and adjuvant-treated contrasts). Reads only the run's CSVs.
# =============================================================================
# Core DTP scores (mirrors CORE_DTP_SCORES in core/config.R; defined locally so
# this file still runs standalone without sourcing config).
.CORE_SCORES  <- c("Up_ssGSEA", "Down_ssGSEA", "Composite_ssGSEA")
.CORE_LABEL   <- c(Up_ssGSEA = "DTP Up", Down_ssGSEA = "DTP Down",
                   Composite_ssGSEA = "DTP Composite")
# Gate constants (match core/config.R; used only for the not-testable annotation
# — the authoritative gate result is the Is_Testable column already in the CSV).
if (!exists("MIN_COX_N")) MIN_COX_N <- 10

# Classify a Project name into a subgroup category, or drop it. Matches only PURE
# cohorts: *_Stage[1-4] (NOT the treatment/MSS-qualified stage subsets) and
# *_All_MSS/MSI. Compound and whole-cohort projects match nothing and are excluded
# (they belong to §3.2 / §3.4).
.subgroup_category <- function(project) {
  if (grepl("_CMS[1-4]$", project))                      return("CMS")
  if (grepl("_PDS[1-3]$", project) || grepl("_Mixed$", project)) return("PDS")
  if (grepl("_Stage[1-4]$", project))                    return("Stage")   # pure stage only
  if (grepl("_All_MSS$|_All_MSI$", project))             return("MSI")     # pure MSI only
  NA_character_
}
.subgroup_label <- function(project) {
  s <- gsub("_", " ", sub("^(GSE|TCGA)_", "", project))
  for (r in list(c("Stage1","Stage I"), c("Stage2","Stage II"),
                 c("Stage3","Stage III"), c("Stage4","Stage IV")))
    s <- sub(r[[1]], r[[2]], s, fixed = TRUE)   # "Stage3" -> "Stage III"
  s
}

# Reconstruct subgroup membership from a clinical frame, mirroring the cohort
# definitions in modules/crc_survival.R, so the within-subgroup score SD can be
# computed to rescale the Cox HR per 1 SD (as in Fig 2). Fails loudly for a kept
# Project it does not recognise rather than silently returning nothing.
.subgroup_members <- function(clin, project, dataset) {
  na0 <- function(x) { x[is.na(x)] <- FALSE; x }
  if (grepl("_CMS[1-4]$", project)) return(na0(clin$CMS == sub(".*_", "", project)))
  if (grepl("_PDS[1-3]$", project)) return(na0(clin$PDS == sub(".*_", "", project)))
  if (grepl("_Mixed$", project))    return(na0(clin$PDS == "Mixed"))
  m <- if (dataset == "GSE39582") switch(project,
      "GSE_All_MSS" = clin$MMR == "pMMR",
      "GSE_All_MSI" = clin$MMR == "dMMR",
      "GSE_Stage1"  = clin$TNM_stage == 1,
      "GSE_Stage2"  = clin$TNM_stage == 2,
      "GSE_Stage3"  = clin$TNM_stage == 3,
      "GSE_Stage4"  = clin$TNM_stage == 4,
      NULL)
    else switch(project,
      "TCGA_All_MSS" = clin$paper_MSI_status == "MSS",
      "TCGA_All_MSI" = clin$paper_MSI_status %in% c("MSI-H", "MSI-L"),
      "TCGA_Stage1"  = clin$Stage == "I",
      "TCGA_Stage2"  = clin$Stage == "II",
      "TCGA_Stage3"  = clin$Stage == "III",
      "TCGA_Stage4"  = clin$Stage == "IV",
      NULL)
  if (is.null(m))
    stop("[group 4] no membership rule for kept subgroup '", project,
         "'; add it to .subgroup_members() (mirror modules/crc_survival.R).")
  na0(m)
}

# Score-across-subtype eps^2 heatmap (Kruskal-Wallis). Rows = core DTP scores,
# column facets = CMS/PDS subtype axis, x within facet = Dataset (all cohorts in
# `sss`, e.g. GSE39582 + TCGA-COAD); fill = eps^2, label = eps^2 (+ raw p / FDR q
# when detailed). Reads the Subtype_Score_Stats rows the caller already loaded.
# Shared by the §3.3 survival figure (standalone Fig3_3B) and Panel B of the
# subtype-violin composite (build_subtype_violin_composite).
.subtype_eps2_heatmap <- function(sss, detailed) {
  .dplyr_local(environment())
  pb <- sss %>% filter(Score %in% .CORE_SCORES, Subtype_Axis %in% c("CMS", "PDS")) %>%
    mutate(Score = factor(.CORE_LABEL[Score], rev(unname(.CORE_LABEL))),
           Axis  = factor(paste0(Subtype_Axis, " subtype"), c("CMS subtype", "PDS subtype")),
           eps2  = suppressWarnings(as.numeric(Eps2)),
           star  = .sig_stars(FDR_P, ""),
           lab_clean = ifelse(star == "", sprintf("%.02f", eps2),
                              paste0(sprintf("%.02f", eps2), " ", star)),
           lab_full  = mapply(.stat_annot, detailed = TRUE, clean = lab_clean,
                              effect = sprintf("eps2=%.2f", eps2),
                              raw_p = suppressWarnings(as.numeric(Raw_P)), fdr = FDR_P))
  pb$lab <- if (detailed) pb$lab_full else pb$lab_clean
  ggplot(pb, aes(x = Dataset, y = Score, fill = eps2)) +
    geom_tile(colour = "white", linewidth = 0.6) +
    geom_text(aes(label = lab), size = if (detailed) 2.4 else 3.0, fontface = "bold",
              colour = "grey10", lineheight = 0.82) +
    facet_grid(. ~ Axis) +
    scale_fill_gradient(low = "white", high = "#762A83", limits = c(0, NA),
      name = expression(epsilon^2~"(effect size)"),
      guide = guide_colourbar(barwidth = 9, barheight = 0.5, title.vjust = 1)) +
    labs(title = "Does the DTP score itself differ across subtypes? (Kruskal-Wallis)",
         subtitle = if (detailed)
           "Epsilon-squared effect size across subtype levels; cells: eps2, raw p, FDR q (BH); * FDR<0.05."
         else
           "Epsilon-squared effect size across subtype levels; * = FDR<0.05 (BH).",
         x = NULL, y = NULL) +
    .base_theme + theme(axis.text.x = element_text(size = 9),
                        panel.spacing.x = unit(6, "pt"))
}

build_subtype_survival_figures <- function(run_dir, out_dir) {
  message("[group 4] subtype/subgroup univariable survival (§3.3)")
  .dplyr_local(environment())

  stats <- .rd(run_dir, "CRC_Statistical_Summary.csv")
  sss   <- .rd(run_dir, "Subtype_Score_Stats.csv")
  need_stats <- c("Dataset","Project","Test","Metric","Score","Raw_P","Effect_r","HR",
                  "HR_lower","HR_upper","C_index","N","N_events","Is_Testable","FDR_P","Is_Significant")
  need_sss   <- c("Dataset","Subtype_Axis","Score","Raw_P","Eps2","N","K","FDR_P","Is_Significant")
  miss1 <- setdiff(need_stats, names(stats)); miss2 <- setdiff(need_sss, names(sss))
  if (length(miss1)) stop("CRC_Statistical_Summary.csv missing column(s): ", paste(miss1, collapse = ", "))
  if (length(miss2)) stop("Subtype_Score_Stats.csv missing column(s): ",   paste(miss2, collapse = ", "))

  clin  <- list(GSE39582 = .rd(run_dir, "GSE39582_clinical.csv"),
                "TCGA-COAD" = .rd(run_dir, "TCGA_COAD_clinical.csv"))
  # metric_cols() is the canonical helper in core/stats.R (sourced before this file).
  # Within-subgroup score SD (non-NA score/time/event), computed once per member set.
  # Core SD of the filtered score vector is .score_sd() (file-local helper above).
  sd_of <- function(ds, project, score, metric) {
    d <- clin[[ds]]; mc <- metric_cols(metric)
    if (!all(c(score, mc[["t"]], mc[["e"]]) %in% names(d))) return(NA_real_)
    mem <- .subgroup_members(d, project, ds); d <- d[mem, , drop = FALSE]
    v <- suppressWarnings(as.numeric(d[[score]]))
    ok <- !is.na(v) & !is.na(d[[mc[["t"]]]]) & !is.na(d[[mc[["e"]]]])
    .score_sd(v[ok])
  }
  # Pull one stat cell (returns NA-filled list if the row is absent).
  cell <- function(ds, project, test, metric, score) {
    s <- stats[stats$Dataset == ds & stats$Project == project & stats$Test == test &
               stats$Metric == metric & stats$Score == score, ]
    if (!nrow(s)) return(NULL)
    s[1, ]
  }

  # Subgroup universe: every Project that classifies into a subgroup category.
  proj <- unique(stats[, c("Dataset", "Project")])
  proj$Category <- vapply(proj$Project, .subgroup_category, character(1))
  proj <- proj[!is.na(proj$Category), , drop = FALSE]
  if (!nrow(proj)) stop("[group 4] no subgroup cohorts found in CRC_Statistical_Summary.csv")
  proj$Subgroup <- vapply(proj$Project, .subgroup_label, character(1))
  proj$dtag     <- ifelse(proj$Dataset == "GSE39582", "GSE", "TCGA")   # sort key only
  # Canonical short cohort name (DATASET_SHORT, core/plotting.R): "Marisa" / "TCGA".
  proj$ylab     <- paste0(unname(DATASET_SHORT[proj$Dataset]), " · ", proj$Subgroup)
  CAT_LEVELS <- intersect(c("CMS","PDS","Stage","MSI"), unique(proj$Category))

  # Endpoint labels: use canonical ENDPOINT_LABEL from core/plotting.R (do not
  # re-declare a local END_LAB duplicate).
  metrics   <- c("OS", "RFS")

  # ---- Assemble long cell table (all core scores; table uses all, Panel A a subset) ----
  rows <- list()
  for (i in seq_len(nrow(proj))) {
    p <- proj[i, ]
    for (sc in .CORE_SCORES) for (mt in metrics) {
      cx <- cell(p$Dataset, p$Project, "Cox",      mt, sc)
      wl <- cell(p$Dataset, p$Project, "Wilcoxon", mt, sc)
      km <- cell(p$Dataset, p$Project, "KM",       mt, sc)
      if (is.null(cx)) next
      testable <- isTRUE(as.logical(cx$Is_Testable))
      sd_sc <- if (testable) sd_of(p$Dataset, p$Project, sc, mt) else NA_real_
      to_sd <- function(x) { x <- suppressWarnings(as.numeric(x));
        if (!testable || is.na(x) || is.na(sd_sc)) NA_real_ else exp(log(x) * sd_sc) }
      hr_sd <- to_sd(cx$HR); lo_sd <- to_sd(cx$HR_lower); hi_sd <- to_sd(cx$HR_upper)
      rows[[length(rows) + 1]] <- tibble(
        Dataset = p$Dataset, Category = p$Category, Subgroup = p$Subgroup, ylab = p$ylab,
        Score = sc, Endpoint = mt, N = as.integer(cx$N), N_events = as.integer(cx$N_events),
        Is_Testable = testable,
        Effect_r = suppressWarnings(as.numeric(wl$Effect_r %||% NA)),
        Wilcoxon_FDR_P = suppressWarnings(as.numeric(wl$FDR_P %||% NA)),
        HR_perSD = hr_sd, HR_lower_perSD = lo_sd, HR_upper_perSD = hi_sd,
        log2HR = ifelse(is.na(hr_sd), NA_real_, log2(hr_sd)),
        Cox_Raw_P = suppressWarnings(as.numeric(cx$Raw_P)),
        Cox_FDR_P = suppressWarnings(as.numeric(cx$FDR_P)),
        Cox_sig = isTRUE(as.logical(cx$Is_Significant)),
        C_index = suppressWarnings(as.numeric(cx$C_index)),
        KM_logrank_FDR_P = suppressWarnings(as.numeric(km$FDR_P %||% NA)))
    }
  }
  cells <- bind_rows(rows)
  if (!nrow(cells)) stop("[group 4] no Cox cells assembled for subgroup cohorts")
  cells$Category <- factor(cells$Category, CAT_LEVELS)
  cells$ylab     <- factor(cells$ylab, levels = rev(unique(proj$ylab[order(proj$Category, proj$dtag, proj$Subgroup)])))

  # =========================== Panel A: tile grid ===========================
  A_SCORES <- c("Up_ssGSEA", "Composite_ssGSEA")   # primary display (Down in table)
  ck_lab <- c(`Up_ssGSEA|OS` = "DTP Up\n(3-yr OS)", `Up_ssGSEA|RFS` = "DTP Up\n(3-yr RFS)",
              `Composite_ssGSEA|OS` = "DTP Comp.\n(3-yr OS)", `Composite_ssGSEA|RFS` = "DTP Comp.\n(3-yr RFS)")
  pa <- cells %>% filter(Score %in% A_SCORES) %>%
    mutate(ckey = factor(paste0(Score, "|", Endpoint), levels = names(ck_lab)),
           tile_clean = ifelse(!Is_Testable, "n/t",
                        ifelse(Cox_sig, .sig_stars(Cox_FDR_P, ""), "")),
           # Full: per-SD HR + raw p + FDR q on every TESTABLE cell; gated stays "n/t".
           .effA = ifelse(is.na(HR_perSD), NA_character_, sprintf("HR=%.2f", HR_perSD)),
           tile_full = ifelse(!Is_Testable, "n/t",
                       mapply(.stat_annot, detailed = TRUE, clean = tile_clean,
                              effect = .effA, raw_p = Cox_Raw_P, fdr = Cox_FDR_P)))
  lim <- min(max(abs(pa$log2HR[is.finite(pa$log2HR)]), na.rm = TRUE), 2)
  make_A <- function(detailed) {
    d <- pa; d$tile <- if (detailed) d$tile_full else d$tile_clean
    ggplot(d, aes(x = ckey, y = ylab, fill = pmax(pmin(log2HR, lim), -lim))) +
      geom_tile(colour = "white", linewidth = 0.6) +
      geom_text(aes(label = tile), colour = ifelse(d$tile == "n/t", "grey45", "grey10"),
                fontface = "bold", size = if (detailed) 2.2 else 3.1, vjust = 0.5, lineheight = 0.82) +
      facet_grid(Category ~ ., scales = "free_y", space = "free_y", switch = "y") +
      scale_fill_gradient2(low = "#3B6EA5", mid = "white", high = "#E64B35", midpoint = 0,
        limits = c(-lim, lim), name = expression(log[2]~"HR (per 1 SD of score)"),
        guide = guide_colourbar(barwidth = 9, barheight = 0.5, title.vjust = 1), na.value = "grey85") +
      scale_x_discrete(labels = ck_lab) +
      labs(title = "DTP score vs 3-year outcome within molecular and clinical subgroups",
           subtitle = if (detailed)
             paste0("Univariable Cox HR per 1 SD (follow-up truncated at 36 mo). Cells: per-SD HR, raw p, FDR q (BH); * FDR<0.05.\n",
                    "Grey \"n/t\" = not testable (< ", MIN_EVENTS, " events / ", MIN_COX_N,
                    " n). Most subgroups are underpowered: a null or n/t cell is NOT evidence of absence.")
           else
             paste0("Univariable Cox HR per 1 SD (follow-up truncated at 36 mo). * = FDR<0.05.\n",
                    "Grey \"n/t\" = not testable (< ", MIN_EVENTS, " events / ", MIN_COX_N,
                    " n), distinct from tested-but-null.\n",
                    "Most subgroups are underpowered: a non-significant or not-testable cell is NOT evidence of absence."),
           x = NULL, y = NULL) +
      .base_theme + theme(axis.text.x = element_text(size = 8, lineheight = 0.85),
                          axis.text.y = element_text(size = 8),
                          strip.placement = "outside",
                          strip.text.y.left = element_text(angle = 0, face = "bold", size = 9),
                          panel.spacing.y = unit(4, "pt"))
  }

  # ---- Save panels (standalone) ----
  # Panel A = subgroup HR matrix (Fig3_3A). Panel B = score-across-subtype eps^2
  # heatmap (Fig3_3B), built by the shared .subtype_eps2_heatmap() helper — that
  # same heatmap is now Panel B of the subtype-violin composite
  # (build_subtype_violin_composite). The former combined
  # "Fig3_subtype_survival_composite" (A over B) was retired: Panel A already
  # stands alone here as Fig3_3A, and the eps^2 heatmap now travels with the
  # violins it summarises.
  n_sub <- nlevels(cells$ylab); a_h <- max(4.5, 0.32 * n_sub + 2.6)
  .save_fig(make_A(TRUE), "Fig3_3A_subgroup_hr_matrix", 8.8, a_h, out_dir)
  .save_fig(.strip_titles(make_A(FALSE)), "Fig3_3A_subgroup_hr_matrix", 8.8, a_h,
            file.path(out_dir, "publication"))
  .save_fig(.subtype_eps2_heatmap(sss, TRUE), "Fig3_3B_score_across_subtype", 7.5, 3.2, out_dir)

  # ================================ Table 3.3 ================================
  rd2 <- function(x) round(x, 2); rd3 <- function(x) round(x, 3)
  tbl <- cells %>%
    transmute(Dataset, Category = as.character(Category), Subgroup, Score, Endpoint,
              N, N_events, Is_Testable,
              Effect_r = rd2(Effect_r), Wilcoxon_FDR_P = rd3(Wilcoxon_FDR_P),
              HR_perSD = rd2(HR_perSD), HR_lower = rd2(HR_lower_perSD), HR_upper = rd2(HR_upper_perSD),
              Cox_FDR_P = rd3(Cox_FDR_P), C_index = rd2(C_index),
              KM_logrank_FDR_P = rd3(KM_logrank_FDR_P)) %>%
    arrange(factor(Category, CAT_LEVELS), Subgroup,
            factor(Score, .CORE_SCORES), factor(Endpoint, metrics))
  readr::write_csv(tbl, file.path(out_dir, "Table_3_3_subtype_survival.csv"))
  message("  wrote Table_3_3_subtype_survival.csv (", nrow(tbl), " rows)")
}

# =============================================================================
# GROUP 4b — subtype-score violins (thesis Section 3.3, Fig3_3C). One 3-row x
#   4-column grid: rows = DTP core score (Up / Down / Composite), columns =
#   cohort x subtype-axis (Marisa CMS, Marisa PDS, TCGA-COAD CMS, TCGA-COAD PDS).
#   Each cell is the per-sample ssGSEA score distribution across that axis's
#   subtype levels for that cohort. CMS shows CMS1-4 (unclassified NA dropped);
#   PDS shows PDS1-3 ("Mixed" dropped, matching the Cox convention elsewhere).
#
#   Unlike the retired eps2/KW panel this composite used to carry (that omnibus
#   test is now standalone as Fig3_3B_score_across_subtype via
#   .subtype_eps2_heatmap(), unchanged), every subtype PAIR within an axis is
#   tested directly: pairwise Mann-Whitney with rank-biserial r
#   (get_wilcox_stats(), core/stats.R), BH-corrected within each
#   (Dataset, Axis, Score) family — the same family the single omnibus KW test
#   used to cover — and written in full to Table_3_3C_subtype_pairwise.csv. Only
#   FDR<0.05 pairs are bracketed on the figure, so the figure and the CSV make
#   identical claims by construction. Rebuildable from a completed run's
#   GSE39582 + TCGA-COAD clinical CSVs alone.
# =============================================================================
.SUBTYPE_LEVELS <- list(CMS = c("CMS1", "CMS2", "CMS3", "CMS4"),
                        PDS = c("PDS1", "PDS2", "PDS3"))

# Pairwise Mann-Whitney + rank-biserial r for one (Dataset, Axis, Score), one row
# per unordered pair of that axis's .SUBTYPE_LEVELS (CMS: 6 pairs; PDS: 3 pairs).
# Pairs are factored in .SUBTYPE_LEVELS order, so get_wilcox_stats() (Group_A =
# first/"x" argument to wilcox.test(), Group_B = second/"y") yields: NEGATIVE r =
# higher score in Group_A (the earlier-numbered subtype); POSITIVE r = higher
# score in Group_B (the later-numbered one) — verified against Median_A/Median_B
# in the emitted table. FDR is left to the caller (family-level BH).
.subtype_pair_stats <- function(clin, dataset, axis, score) {
  lv <- .SUBTYPE_LEVELS[[axis]]
  d  <- clin[clin[[axis]] %in% lv, c(axis, score)]
  d[[score]] <- suppressWarnings(as.numeric(d[[score]]))
  d  <- d[!is.na(d[[score]]), , drop = FALSE]
  rows <- lapply(utils::combn(lv, 2, simplify = FALSE), function(pr) {
    dd <- d[d[[axis]] %in% pr, , drop = FALSE]
    dd[[axis]] <- factor(dd[[axis]], levels = pr)
    a  <- dd[[score]][dd[[axis]] == pr[1]]; b <- dd[[score]][dd[[axis]] == pr[2]]
    ws <- get_wilcox_stats(dd, score, axis)
    data.frame(Dataset = dataset, Subtype_Axis = axis, Score = score,
               Group_A = pr[1], Group_B = pr[2], N_A = length(a), N_B = length(b),
               Median_A = if (length(a)) median(a) else NA_real_,
               Median_B = if (length(b)) median(b) else NA_real_,
               Effect_r = ws$r, Raw_P = ws$p, Is_Testable = !is.na(ws$p),
               stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}

# Nested significance brackets (narrowest span first, stacked upward) for the
# FDR<0.05 pairs of one cell. `lv` gives x-position order (match(level, lv));
# `y0`/`span` are the cell's base data range top / width, before any bracket
# headroom. Returns the annotate() layers plus the ylim top the cell needs to
# fit every bracket it draws — the source of the plan's per-cell headroom.
.violin_brackets <- function(sig_df, lv, y0, span, detailed) {
  if (!nrow(sig_df)) return(list(layers = list(), ylim_top = y0 + span * 0.06))
  sig_df$xa <- match(sig_df$Group_A, lv); sig_df$xb <- match(sig_df$Group_B, lv)
  sig_df <- sig_df[order(abs(sig_df$xb - sig_df$xa)), ]
  step   <- span * 0.11; tick_h <- step * 0.22
  layers <- list()
  for (i in seq_len(nrow(sig_df))) {
    r <- sig_df[i, ]; y <- y0 + i * step
    lab <- if (detailed) paste0(.sig_stars(r$FDR_P, ""), "\n", .fmt_p("q", r$FDR_P))
           else .sig_stars(r$FDR_P, "")
    layers[[length(layers) + 1]] <- ggplot2::annotate("segment",
      x = r$xa, xend = r$xb, y = y, yend = y, linewidth = 0.4, colour = "grey20")
    layers[[length(layers) + 1]] <- ggplot2::annotate("segment",
      x = c(r$xa, r$xb), xend = c(r$xa, r$xb), y = y, yend = y - tick_h,
      linewidth = 0.4, colour = "grey20")
    layers[[length(layers) + 1]] <- ggplot2::annotate("text",
      x = (r$xa + r$xb) / 2, y = y, label = lab, vjust = -0.15,
      size = if (detailed) 2.5 else 3, fontface = "bold", colour = "grey10",
      lineheight = 0.85)
  }
  list(layers = layers, ylim_top = y0 + nrow(sig_df) * step + step * 0.9)
}

build_subtype_violin_composite <- function(run_dir, out_dir) {
  message("[group 4b] subtype-score violins (Section 3.3)")
  .dplyr_local(environment())

  clin_list <- list(GSE39582 = .rd(run_dir, "GSE39582_clinical.csv"),
                    `TCGA-COAD` = .rd(run_dir, "TCGA_COAD_clinical.csv"))
  axes <- names(.SUBTYPE_LEVELS)                    # CMS, PDS (within-cohort column order)
  for (ds in names(clin_list)) {
    miss <- setdiff(c(.CORE_SCORES, axes), names(clin_list[[ds]]))
    if (length(miss))
      stop(ds, "_clinical.csv missing column(s): ", paste(miss, collapse = ", "))
  }

  # ---- Pairwise Mann-Whitney table (the figure can only bracket what's here) ----
  pw <- do.call(rbind, lapply(names(clin_list), function(ds)
    do.call(rbind, lapply(axes, function(ax)
      do.call(rbind, lapply(.CORE_SCORES, function(sc)
        .subtype_pair_stats(clin_list[[ds]], ds, ax, sc)))))))
  pw <- pw %>%
    group_by(Dataset, Subtype_Axis, Score) %>%
    mutate(FDR_P = p.adjust(Raw_P, method = "BH")) %>%
    ungroup() %>%
    mutate(FDR_P_global = p.adjust(Raw_P, method = "BH"),
           Is_Significant = !is.na(FDR_P) & FDR_P < 0.05)
  out_tbl <- pw %>%
    transmute(Dataset, Subtype_Axis, Score, Group_A, Group_B, N_A, N_B,
              Median_A = round(Median_A, 4), Median_B = round(Median_B, 4),
              Effect_r = round(Effect_r, 3), Raw_P, FDR_P, FDR_P_global,
              Is_Testable, Is_Significant)
  readr::write_csv(out_tbl, file.path(out_dir, "Table_3_3C_subtype_pairwise.csv"))
  message("  wrote Table_3_3C_subtype_pairwise.csv (", nrow(out_tbl), " rows)")

  # Grid columns, Axis varying fastest within a Dataset: (Marisa,CMS), (Marisa,PDS),
  # (TCGA-COAD,CMS), (TCGA-COAD,PDS) — the 4-column layout the plan specifies.
  cols <- expand.grid(Axis = axes, Dataset = names(clin_list), stringsAsFactors = FALSE)

  # Base y-range shared PER (Dataset, Score): the two cohorts sit on different
  # ssGSEA scales, so a grid-wide (or row-wide) window would flatten one of them;
  # the two axis columns of the SAME cohort+score are directly comparable, so they
  # share this base range. Per-cell bracket headroom is layered on in make_cell().
  base_ylim <- function(dataset, score) {
    v <- suppressWarnings(as.numeric(clin_list[[dataset]][[score]])); v <- v[!is.na(v)]
    rng <- range(v); pad <- diff(rng) * 0.04
    c(rng[1] - pad, rng[2] + pad)
  }

  # One violin cell. Column header (top row only) is a plot.title so it is
  # dropped by .strip_titles for the publication copy — the CMS1-4 / PDS1-3 x
  # ticks self-identify the column; the row identity is the left-most column's
  # y-axis title, which survives stripping. Brackets are annotate() layers, so
  # they are untouched by .strip_titles and appear identically in both copies.
  make_cell <- function(score, dataset, axis, detailed, is_top_row, is_left_col) {
    lv <- .SUBTYPE_LEVELS[[axis]]
    cl <- clin_list[[dataset]]
    d  <- cl[cl[[axis]] %in% lv, c(axis, score)]
    d[[score]] <- suppressWarnings(as.numeric(d[[score]]))
    d  <- d[!is.na(d[[score]]), , drop = FALSE]
    d[[axis]] <- factor(d[[axis]], levels = lv)
    tb   <- table(d[[axis]])
    lbls <- setNames(paste0(names(tb), "\n(n=", as.integer(tb), ")"), names(tb))

    bylim <- base_ylim(dataset, score)
    sig   <- pw[pw$Dataset == dataset & pw$Subtype_Axis == axis & pw$Score == score &
                pw$Is_Significant, , drop = FALSE]
    br    <- .violin_brackets(sig, lv, bylim[2], diff(bylim), detailed)

    ggplot(d, aes(x = .data[[axis]], y = .data[[score]], fill = .data[[axis]])) +
      geom_violin(alpha = 0.8, trim = FALSE, scale = "width",
                  linewidth = 0.3, colour = "grey30") +
      geom_boxplot(width = 0.15, fill = "white", alpha = 0.9,
                   outlier.shape = NA, linewidth = 0.3) +
      br$layers +
      scale_fill_brewer(palette = "Set2") +
      scale_x_discrete(labels = lbls) +
      coord_cartesian(ylim = c(bylim[1], br$ylim_top)) +
      labs(title = if (is_top_row) sprintf("%s · %s", DATASET_FULL[[dataset]], axis) else NULL,
           x = NULL,
           y = if (is_left_col) paste0(.CORE_LABEL[[score]], " ssGSEA") else NULL) +
      .base_theme +
      theme(legend.position = "none",
            plot.title  = element_text(face = "bold", size = 11, hjust = 0.5),
            axis.text.x = element_text(size = 7.5, lineheight = 0.85))
  }

  # The 3x4 violin grid (a patchwork of the 12 cells) — no A/B panel tags; this
  # is now the whole figure, not one panel of a taller composite.
  make_grid <- function(detailed) {
    cells <- list()
    for (sc in .CORE_SCORES) {
      is_top <- sc == .CORE_SCORES[[1]]
      for (i in seq_len(nrow(cols))) {
        is_left <- cols$Dataset[i] == names(clin_list)[[1]] && cols$Axis[i] == axes[[1]]
        cells[[length(cells) + 1]] <- make_cell(
          sc, cols$Dataset[i], cols$Axis[i], detailed, is_top, is_left)
      }
    }
    wrap_plots(cells, ncol = nrow(cols), byrow = TRUE)
  }

  CAP <- paste0(
    "Marisa (GSE39582) and TCGA-COAD; CMS and PDS subtype axes. CMS-unclassified and PDS \"Mixed\" are excluded.\n",
    "Pairwise Mann-Whitney tests with rank-biserial r, BH-corrected within each cohort x axis x score;\n",
    "only FDR<0.05 pairs are bracketed. Full pairwise table: Table_3_3C_subtype_pairwise.csv.")

  make_composite <- function(detailed) {
    grid <- make_grid(detailed)
    if (!detailed) grid <- .strip_titles(grid)
    grid +
      plot_annotation(
        title = if (detailed)
          "DTP signature scores across molecular subtypes: Marisa (GSE39582) and TCGA-COAD" else NULL,
        caption = CAP,
        theme = theme(plot.title   = element_text(face = "bold", size = 15, hjust = 0),
                      plot.caption = element_text(size = 8, colour = "grey30", hjust = 0)))
  }

  .save_fig(make_composite(TRUE),  "Fig3_3C_subtype_violins", 13, 15, out_dir)
  .save_fig(make_composite(FALSE), "Fig3_3C_subtype_violins", 13, 15, file.path(out_dir, "publication"))
}

# =============================================================================
# GROUP 5 — confounding: does the DTP score stay prognostic after adjustment?
#           (thesis Section 3.4, the confounding half — Figure 4 covers effect
#           modification). Reads only Adjusted_Cox_Summary.csv from the run.
#           ONE heatmap (Dataset x endpoint rows, adjustment-model columns):
#             fill   = % change in the score's log-HR after adjustment
#                      (diverging; blue attenuated toward null, red strengthened;
#                       Greenland confounding magnitude — within +/-10% ~ white),
#             number = adjusted-model BH-FDR, likelihood-ratio test of the score
#                      vs the covariate-only model (residual significance;
#                      * FDR<0.05, grey n/t = not testable). Colour and number
#                      encode different quantities — see the legend + subtitle.
# =============================================================================
build_confounding_figure <- function(run_dir, out_dir, primary_score = "Up_ssGSEA") {
  message("[group 5] confounding figure")
  .dplyr_local(environment())

  MODEL_LEVELS  <- c("Clinicopath", "CMS_adjusted", "PDS_adjusted")
  MODEL_LABEL   <- c(Clinicopath = "Clinicopathology", CMS_adjusted = "+CMS", PDS_adjusted = "+PDS")
  DS_LEVELS     <- c("GSE39582", "GSE39582 (treated)", "TCGA-COAD")
  # Row-strip labels = canonical DATASET_FULL; y-axis = canonical ENDPOINT_LABEL (core/plotting.R).
  METRIC_LEVELS <- c("OS", "RFS")
  score_lab     <- sub("_ssGSEA$", "", primary_score)   # "Up" / "Composite" / "Down"

  adj <- .rd(run_dir, "Adjusted_Cox_Summary.csv")
  num <- function(x) suppressWarnings(as.numeric(x))
  d <- adj %>%
    filter(Score == primary_score, Model %in% MODEL_LEVELS,
           Dataset %in% DS_LEVELS, Metric %in% METRIC_LEVELS) %>%
    mutate(Model   = factor(Model,   MODEL_LEVELS),
           Dataset = factor(Dataset, DS_LEVELS),
           Metric  = factor(Metric,  METRIC_LEVELS),
           FDR_P = num(FDR_P), LRT_score_P = num(LRT_score_P),
           Delta_logHR_pct = num(Delta_logHR_pct),
           Is_Testable = as.logical(Is_Testable))
  if (!nrow(d)) stop("[group 5] no rows for score '", primary_score,
                     "' in Adjusted_Cox_Summary.csv")

  # ------------------------------- the heatmap -------------------------------
  # Tile colour = % change in the score's log-HR after adjustment (Greenland
  # confounding magnitude); cell number = adjusted-model FDR (residual
  # significance). Untestable models (events-per-variable gate) show grey "n/t"
  # rather than a number; none occur in the current data.
  pa <- d %>%
    mutate(lab_clean = ifelse(!Is_Testable, "n/t",
                        paste0(sprintf("%.3f", FDR_P), .sig_stars(FDR_P))),
           # Full: the adjustment %-shift is the native, interpretable effect size
           # (the raw adjusted Cox HR is per unit of log2 ssGSEA and uninterpretable,
           # so it is NOT shown); plus the raw LRT p and the BH-FDR q. Gated -> "n/t".
           .effF = ifelse(is.na(Delta_logHR_pct), NA_character_,
                          sprintf("d%+.0f%%", Delta_logHR_pct)),
           lab_full = ifelse(!Is_Testable, "n/t",
                      mapply(.stat_annot, detailed = TRUE, clean = "",
                             effect = .effF, raw_p = LRT_score_P, fdr = FDR_P)))
  # Symmetric limits so 0% sits at white and +/-10% renders near-white. Computed
  # BEFORE the labels because the text colour has to track this scale, not a fixed
  # percentage: with a data-driven `lim`, an absolute cutoff would put white text on
  # a near-white tile as soon as the run's largest shift grew.
  lim <- ceiling(max(abs(d$Delta_logHR_pct), na.rm = TRUE) / 10) * 10
  dark_at <- 0.65 * lim   # fraction of the ramp beyond which the fill is dark enough
  pa <- pa %>%
    # Dark tiles occur at BOTH colour extremes (strong attenuation OR strong
    # strengthening), so switch to white text on either dark end.
    mutate(lab_col = ifelse(!Is_Testable, "grey45",
                     ifelse(!is.na(Delta_logHR_pct) & abs(Delta_logHR_pct) > dark_at,
                            "white", "grey10")))
  make_fig <- function(detailed) {
    dd <- pa; dd$lab <- if (detailed) dd$lab_full else dd$lab_clean
    ggplot(dd, aes(x = Model, y = fct_rev(Metric), fill = Delta_logHR_pct)) +
      geom_tile(colour = "white", linewidth = 0.6) +
      geom_text(aes(label = lab), colour = dd$lab_col, fontface = "bold",
                size = if (detailed) 2.1 else 2.9, lineheight = 0.8) +
      facet_grid(Dataset ~ ., switch = "y", labeller = labeller(Dataset = DATASET_FULL)) +
      scale_fill_gradient2(low = "#3B6EA5", mid = "white", high = "#E64B35", midpoint = 0,
        limits = c(-lim, lim), na.value = "grey85",
        name = "Change in score log-HR after adjustment (%)  [-] attenuated  <->  [+] strengthened",
        guide = guide_colourbar(barwidth = 12, barheight = 0.5, title.position = "top",
                                title.hjust = 0.5)) +
      scale_x_discrete(labels = MODEL_LABEL) +
      scale_y_discrete(labels = ENDPOINT_LABEL) +
      labs(title = "Independence of the DTP score after adjustment",
           subtitle = if (detailed)
             paste0("Tile colour = % change in the DTP ", score_lab,
                    " score log-HR after adjustment (blue attenuated, red strengthened).\n",
                    "Cells: %-shift (effect), raw LRT p, BH-FDR q; * FDR<0.05; n/t = not testable.\n",
                    "Per 1 SD, follow-up truncated at 36 mo.  |shift| > 10% ~ meaningful confounding (Greenland).  Exploratory.")
           else
             paste0("Tile colour = % change in the DTP ", score_lab,
                    " score log-HR after adjustment (blue = attenuated toward null, red = strengthened).\n",
                    "Number = adjusted-model LRT BH-FDR (score vs covariate-only model); * marks FDR < 0.05.\n",
                    "Per 1 SD, follow-up truncated at 36 mo.  |shift| > 10% ~ meaningful confounding (Greenland).  Exploratory."),
           x = NULL, y = NULL) +
      .base_theme +
      theme(axis.text.x = element_text(size = 9),
            legend.title = element_text(size = 8.5),
            plot.title.position = "plot",   # title/subtitle span the full width, not just the panel
            strip.placement = "outside",
            strip.text.y.left = element_text(angle = 0, face = "bold", size = 9),
            panel.spacing.y = unit(4, "pt"))
  }

  .save_fig(make_fig(TRUE), "Fig5_confounding_summary", 9.5, 5.6, out_dir)
  .save_fig(.strip_titles(make_fig(FALSE)), "Fig5_confounding_summary", 9.5, 5.6,
            file.path(out_dir, "publication"))

  # Caption block alongside the figure (mirrors the §3.4 exploratory framing).
  cap <- c(
    "Figure 5. Independence of the DTP score after adjustment (Section 3.4; exploratory).",
    "",
    paste0("Adjusted-Cox confounding heatmap. Each cell is the DTP ", score_lab,
           " score term from a Cox model (follow-up truncated at 36 months) of overall"),
    "survival (OS) or recurrence-free survival (RFS), adjusted for clinicopathology (stage + microsatellite",
    "status), CMS subtype, or PDS subtype, in GSE39582 (whole cohort and adjuvant-treated subset) and TCGA-COAD.",
    "",
    "Tile colour is the percent change in the score's log-hazard ratio (the Cox coefficient) after adjustment,",
    "relative to the unadjusted model: blue = attenuated toward the null (the adjustor absorbs part of the",
    "signal), red = strengthened, white = unchanged. Following Greenland's rule of thumb a shift beyond +/-10%",
    "is treated as materially confounded; cells within +/-10% therefore render near-white.",
    "",
    "The cell number is the Benjamini-Hochberg FDR of that score term's likelihood-ratio test (the adjusted model",
    "with vs. without the score; residual significance after adjustment) -- the same nested-model test the",
    "effect-modification analysis uses, and more robust than the Wald test for the large per-SD hazard ratios",
    "here. An asterisk marks FDR < 0.05. The FDR family is the three DTP scores (Up / Composite / Down) within",
    "each cohort x endpoint x model. A grey \"n/t\" cell denotes a model the events-per-variable gate (>=10 events",
    "per parameter, n>=10, >=5 events) never fitted (not modelled, as opposed to fitted-but-null); all cells here",
    "were testable.",
    "",
    "Exploratory analysis: adjustment removes only measured confounders, so residual confounding remains",
    "possible; multiple-testing correction reduces but does not eliminate false positives.")
  writeLines(cap, file.path(out_dir, "Fig5_confounding_summary_caption.txt"))
  message("  wrote Fig5_confounding_summary_caption.txt")
}

# =============================================================================
# Orchestrator — build all figure groups from one run directory
# =============================================================================
build_crc_composites <- function(run_dir, out_dir = file.path(run_dir, "composites")) {
  stopifnot(dir.exists(run_dir))
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  message("Building CRC composite figures\n  run: ", run_dir, "\n  out: ", out_dir)
  # Each group is guarded so a missing input for one does not sink the others.
  # A failure is BOTH messaged (in-line) and warned (so it lands in the end-of-run
  # "Warning messages:" summary) — a silent skip must not slip past unnoticed.
  groups <- list("1 (outcomes)"       = build_outcome_figures,
                 "2 (KM)"              = build_km_figures,
                 "3 (subgroups)"       = build_subgroup_figures,
                 "4 (subtype survival)" = build_subtype_survival_figures,
                 "4b (subtype violins)" = build_subtype_violin_composite,
                 "5 (confounding)"     = build_confounding_figure)
  status <- setNames(logical(length(groups)), names(groups))
  for (nm in names(groups)) {
    status[[nm]] <- tryCatch({ groups[[nm]](run_dir, out_dir); TRUE },
      error = function(e) {
        msg <- sprintf("[composite figures] group %s FAILED: %s", nm, conditionMessage(e))
        message("  ", msg); warning(msg, call. = FALSE)
        FALSE
      })
  }
  n_ok <- sum(status); n <- length(status)
  summ <- sprintf("Composite figures: %d/%d groups OK%s", n_ok, n,
    if (n_ok < n) paste0(" — FAILED: ", paste(names(status)[!status], collapse = ", ")) else "")
  message(summ)
  if (n_ok < n) warning(summ, call. = FALSE)   # surface partial failure at run end
  message("CRC composite figures done -> ", out_dir)
  invisible(status)
}

# ==============================================================================
# PART 2 — Pan-cancer composites (build_pancan_composites)
# ==============================================================================
# Reuse the CRC forest display window (log-scale HR truncation) if not present.
if (!exists("FOREST_XLIM")) FOREST_XLIM <- c(0.1, 10)

.PANCAN_DIR_PAL <- c("Higher hazard (HR>1)" = "#E64B35", "Lower hazard (HR<1)" = "#3B6EA5")
.PANCAN_CORR_LAB <- c(Up_ssGSEA_Corrected = "DTP Up", Down_ssGSEA_Corrected = "DTP Down",
                      Composite_ssGSEA_Corrected = "DTP Composite")
.PANCAN_HIGHLIGHT <- "TCGA-COAD"          # colorectal tissue — locate it among the 24
.PANCAN_ENDPT_LAB <- c(OS = "Overall survival", RFS = "Recurrence")
.PANCAN_METRIC_EVT <- c(OS = "OS3Y_event", RFS = "RFS3Y_event")
.PANCAN_METRIC_TIM <- c(OS = "OS3Y_delay", RFS = "RFS3Y_delay")

# ---- helpers -----------------------------------------------------------------
.pancan_num <- function(x) suppressWarnings(as.numeric(x))

# Within-cohort SD of one score over the SAME subset get_cox_stats() models:
# non-NA score AND non-NA truncated follow-up time AND event for that metric.
.pancan_sd <- function(master, project, score, metric) {
  ev <- .PANCAN_METRIC_EVT[[metric]]; tm <- .PANCAN_METRIC_TIM[[metric]]
  d  <- master[master$Project_ID == project, , drop = FALSE]
  v  <- .pancan_num(d[[score]]); t <- .pancan_num(d[[tm]]); e <- .pancan_num(d[[ev]])
  ok <- !is.na(v) & !is.na(t) & !is.na(e)
  if (sum(ok) < 2) return(NA_real_)
  stats::sd(v[ok])
}

# Add per-SD HR columns (point + both bounds) + direction + off-scale squish,
# mirroring the CRC forest. `sd_vec` is aligned row-for-row with `d`.
# Non-convergent fits (non-finite per-SD HR or CI) are NOT squished to the
# FOREST_XLIM boundary — display coords stay NA so geoms drop them (na.rm).
.pancan_rescale <- function(d, sd_vec) {
  hr <- .pancan_num(d$HR); lo <- .pancan_num(d$HR_lower); hi <- .pancan_num(d$HR_upper)
  ok <- !is.na(hr) & !is.na(sd_vec)
  d$hr_sd <- ifelse(ok, exp(log(hr) * sd_vec), NA_real_)
  d$lo_sd <- ifelse(ok, exp(log(lo) * sd_vec), NA_real_)
  d$hi_sd <- ifelse(ok, exp(log(hi) * sd_vec), NA_real_)
  nonconv <- (!is.na(d$hr_sd) & !is.finite(d$hr_sd)) |
             (!is.na(d$lo_sd) & !is.finite(d$lo_sd)) |
             (!is.na(d$hi_sd) & !is.finite(d$hi_sd))
  if (any(nonconv))
    message(sprintf(paste0("[pancan forest] %d non-convergent Cox row(s) not ",
                           "squished to axis boundary (non-finite per-SD HR or CI)"),
                    sum(nonconv)))
  d$dir      <- ifelse(is.na(d$hr_sd) | nonconv, NA_character_,
                       ifelse(d$hr_sd > 1, "Higher hazard (HR>1)", "Lower hazard (HR<1)"))
  d$offscale <- !nonconv & !is.na(d$hr_sd) &
                (d$hr_sd < FOREST_XLIM[1] | d$hr_sd > FOREST_XLIM[2])
  d$hr_d <- ifelse(nonconv, NA_real_, scales::squish(d$hr_sd, FOREST_XLIM))
  d$lo_d <- ifelse(nonconv, NA_real_, scales::squish(d$lo_sd, FOREST_XLIM))
  d$hi_d <- ifelse(nonconv, NA_real_, scales::squish(d$hi_sd, FOREST_XLIM))
  d$nonconv <- nonconv
  d
}

# Cohort row order shared by Panel A (Up forest) and Panel B (Composite/Down tile)
# so the panels align: ascending Up overall-survival per-SD HR, not-testable OS
# (NA) sorted to the bottom via na.last = FALSE.
# `order_score` is the Up column NAME, which is caller-supplied (run_pancan_treated
# takes an `up_name`), so it must not be hardcoded: with no matching rows the order
# silently degrades to CSV order while the subtitle and caption still claim the
# ascending-Up-OS-HR ordering. Warn loudly rather than fail quietly.
.pancan_cohort_order <- function(fdr, master, order_score = "Up_ssGSEA") {
  .dplyr_local(environment())
  d <- fdr %>% filter(Family == "PerCohort", Test == "Cox",
                      Score == order_score, Metric == "OS")
  if (!nrow(d))
    warning("[pancan composites] no PerCohort/Cox/OS rows for ordering score '",
            order_score, "'; cohort row order falls back to input order and is NOT ",
            "the ascending Up-OS per-SD HR order the figure claims.", call. = FALSE)
  d <- .pancan_rescale(d, mapply(function(pr, mt) .pancan_sd(master, pr, order_score, mt),
                                 d$Project, d$Metric))
  ord <- d$Project[order(d$hr_sd, na.last = FALSE)]
  unique(c(ord, setdiff(unique(fdr$Project[fdr$Family == "PerCohort"]), ord)))
}

# ==============================================================================
# PER-COHORT FOREST for one score (per 1 SD), OS | RFS facets. Used for all three
# scores (Fig 6 panels A = Up, B = Composite, C = Down); rows are ALWAYS ordered by
# the Up OS per-SD HR (shared helper) so the three panels line up tissue-for-tissue.
# ==============================================================================
build_pancan_forest <- function(fdr, master, out_dir, score = "Up_ssGSEA",
                                score_lab = "DTP Up", save_name = "Fig6A_pancan_forest_Up",
                                order_score = "Up_ssGSEA") {
  .dplyr_local(environment())
  SCORE <- score
  d <- fdr %>%
    filter(Family == "PerCohort", Test == "Cox", Score == SCORE, Metric %in% c("OS", "RFS")) %>%
    mutate(Is_Testable = as.logical(Is_Testable), Is_Significant = as.logical(Is_Significant),
           FDR_P = .pancan_num(FDR_P), Raw_P = .pancan_num(Raw_P),
           N = .pancan_num(N), N_events = .pancan_num(N_events))
  sd_vec <- mapply(function(pr, mt) .pancan_sd(master, pr, SCORE, mt), d$Project, d$Metric)
  d <- .pancan_rescale(d, sd_vec)
  # Single-level row facet carrying the score name, so every panel is
  # self-identifying (as a bold grey strip down the right edge) even after
  # .strip_titles() drops the plot title/subtitle for the publication copy.
  d$ScoreLab <- factor(score_lab)

  # Row order = ascending Up OS per-SD HR (shared helper; ggplot puts the first
  # level at the bottom, so the highest-hazard cohort lands at the top). ALL three
  # score panels reuse this same order so the tissues line up across A/B/C.
  d$Project <- factor(d$Project, levels = .pancan_cohort_order(fdr, master, order_score))
  d$Metric  <- factor(d$Metric, c("OS", "RFS"))

  # TCGA-COAD highlight: pale band behind the row (both facets) + coloured/bold label.
  coad_pos <- match(.PANCAN_HIGHLIGHT, levels(d$Project))
  lab_cols <- ifelse(levels(d$Project) == .PANCAN_HIGHLIGHT, "#B2182B", "grey20")
  lab_face <- ifelse(levels(d$Project) == .PANCAN_HIGHLIGHT, "bold", "plain")
  star_x   <- FOREST_XLIM[2] * 0.9

  # Left-edge label: n / e (clean) plus the per-SD HR + CI (full). Right-margin:
  # bare "*" for FDR<0.05 (clean) vs raw p + FDR q on every testable row (full).
  d$ne_clean <- sprintf("n=%d, e=%d", as.integer(d$N), as.integer(d$N_events))
  d$ne_full  <- ifelse(is.na(d$hr_sd), d$ne_clean,
                       sprintf("%s\nHR=%.2f [%.2f-%.2f]", d$ne_clean, d$hr_sd, d$lo_sd, d$hi_sd))
  d$rm_full  <- mapply(function(rp, fp)
                   paste(stats::na.omit(c(.fmt_p("p", rp),
                         paste0(.fmt_p("q", fp), .sig_stars(fp, "")))), collapse = "  "),
                   d$Raw_P, d$FDR_P)

  make_fig <- function(detailed) {
    dd <- d; dd$ne <- if (detailed) dd$ne_full else dd$ne_clean
    p <- ggplot(dd, aes(x = hr_d, y = Project))
    if (!is.na(coad_pos))
      p <- p + annotate("rect", xmin = FOREST_XLIM[1], xmax = FOREST_XLIM[2],
                        ymin = coad_pos - 0.5, ymax = coad_pos + 0.5,
                        fill = "#FDE7C9", alpha = 0.7)
    p <- p +
      geom_vline(xintercept = 1, linetype = "dashed", colour = "grey55") +
      geom_errorbarh(aes(xmin = lo_d, xmax = hi_d, colour = dir), height = 0.3,
                     linewidth = 0.5, na.rm = TRUE) +
      geom_point(aes(colour = dir, shape = offscale), size = 2.0, na.rm = TRUE) +
      # per-cohort n / events (+ per-SD HR & CI when detailed) at the left edge
      geom_text(data = subset(dd, Is_Testable & !is.na(N)), inherit.aes = FALSE,
                aes(x = FOREST_XLIM[1], y = Project, label = ne),
                hjust = 0, vjust = -0.6, size = 1.9, colour = "grey45", lineheight = 0.82) +
      # grey n/t marker for not-testable cohorts (TCGA-TGCT OS today)
      geom_text(data = subset(dd, !Is_Testable), inherit.aes = FALSE,
                aes(x = 1, y = Project, label = "n/t"),
                colour = "grey55", size = 2.6, fontface = "italic")
    if (detailed)   # raw p + FDR q on every testable row, right margin
      p <- p + geom_text(data = subset(dd, Is_Testable & !is.na(FDR_P)), inherit.aes = FALSE,
                         aes(x = FOREST_XLIM[2], y = Project, label = rm_full),
                         hjust = 1, vjust = -0.6, size = 1.7, colour = "grey30")
    else            # bare FDR<0.05 star, right margin
      p <- p + geom_text(data = subset(dd, !is.na(FDR_P) & FDR_P < 0.05), inherit.aes = FALSE,
                         aes(x = star_x, y = Project, label = "*"),
                         size = 4.6, fontface = "bold", colour = "grey15")
    p +
      facet_grid(ScoreLab ~ Metric, labeller = labeller(Metric = .PANCAN_ENDPT_LAB)) +
      scale_x_log10(limits = FOREST_XLIM, breaks = c(0.1, 0.3, 1, 3, 10),
                    labels = c("0.1", "0.3", "1", "3", "10")) +
      scale_y_discrete(labels = function(x) ifelse(x == .PANCAN_HIGHLIGHT,
                                                   paste0(x, " (CRC)"), x)) +
      scale_colour_manual(values = .PANCAN_DIR_PAL, name = "Score effect",
                          na.value = "grey70", na.translate = FALSE, drop = FALSE) +
      scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 21),
                         labels = c(`FALSE` = "in range", `TRUE` = "off-scale (squished)"),
                         name = NULL) +
      labs(title = sprintf("%s score vs survival across %d TCGA cohorts", score_lab, nlevels(d$Project)),
           subtitle = if (detailed) paste0(
             "Univariable continuous-score Cox HR per 1 SD within each cohort (95% CI); x truncated to [0.1, 10].\n",
             "Left: n/e + per-SD HR [CI]; right: raw p, FDR q (BH); * FDR<0.05; grey n/t = not testable. TCGA-COAD (CRC) highlighted.")
           else paste0(
             "Univariable continuous-score Cox HR per 1 SD within each cohort (95% CI); x truncated to [0.1, 10].\n",
             "Red HR>1 / blue HR<1; * FDR<0.05 (BH); grey n/t = not testable. TCGA-COAD (CRC tissue) highlighted."),
           x = "Hazard ratio per 1 SD (log scale)", y = NULL) +
      .base_theme +
      theme(axis.text.y = element_text(size = 7.5, colour = lab_cols, face = lab_face),
            panel.grid.major.x = element_line(colour = "grey93", linewidth = 0.3),
            panel.spacing.x = unit(9, "pt"))
  }

  .save_fig(make_fig(TRUE), save_name, 9.0, 8.0, out_dir)
  make_fig   # return the closure so the composite can build titled vs clean copies
}

# ==============================================================================
# APPENDIX — batch-corrected pan-cancer aggregate Cox forest (Up/Down/Composite,
#            per 1 SD). Moved out of the Fig 6 composite; saved standalone.
# ==============================================================================
build_pancan_aggregate <- function(fdr, master, out_dir) {
  .dplyr_local(environment())

  # Reproduce the module's corrected scores (Methods 2.8): removeBatchEffect with
  # cancer type as the batch, on the master [Up, Down, Composite] matrix. Verified
  # to reproduce the module exactly (corrected Up pooled SD ~ 0.064).
  raw_cols <- c("Up_ssGSEA", "Down_ssGSEA", "Composite_ssGSEA")
  sm <- t(as.matrix(sapply(master[, raw_cols], .pancan_num)))
  corrected <- limma::removeBatchEffect(sm, batch = master$Project_ID)
  corr <- as.data.frame(t(corrected))
  names(corr) <- paste0(raw_cols, "_Corrected")
  message(sprintf("[group B] corrected Up pooled SD = %.4f (module target ~0.064)",
                  stats::sd(corr$Up_ssGSEA_Corrected)))
  for (mc in c("OS3Y_event", "OS3Y_delay", "RFS3Y_event", "RFS3Y_delay"))
    corr[[mc]] <- .pancan_num(master[[mc]])

  agg_sd <- function(score, metric) {
    ev <- .PANCAN_METRIC_EVT[[metric]]; tm <- .PANCAN_METRIC_TIM[[metric]]
    v <- corr[[score]]; ok <- !is.na(v) & !is.na(corr[[tm]]) & !is.na(corr[[ev]])
    if (sum(ok) < 2) NA_real_ else stats::sd(v[ok])
  }

  d <- fdr %>%
    filter(Family == "PanCancer", Test == "Cox", Metric %in% c("OS", "RFS")) %>%
    mutate(C_index = .pancan_num(C_index), FDR_P = .pancan_num(FDR_P), Raw_P = .pancan_num(Raw_P),
           Is_Significant = as.logical(Is_Significant),
           ScoreLab = .PANCAN_CORR_LAB[Score])
  sd_vec <- mapply(function(sc, mt) agg_sd(sc, mt), d$Score, d$Metric)
  d <- .pancan_rescale(d, sd_vec)
  d$Metric   <- factor(d$Metric, c("OS", "RFS"))
  d$ScoreLab <- factor(d$ScoreLab, rev(c("DTP Up", "DTP Down", "DTP Composite")))

  # Verification print: aggregate per-SD HRs
  for (i in order(d$Metric, d$ScoreLab))
    message(sprintf("  [aggregate] %-4s %-14s per-unit HR=%9.3f  ->  per-SD HR=%.3f  C=%.3f  FDR=%.3f%s",
                    as.character(d$Metric[i]), as.character(d$ScoreLab[i]),
                    .pancan_num(d$HR[i]), d$hr_sd[i], d$C_index[i], d$FDR_P[i],
                    ifelse(isTRUE(d$Is_Significant[i]), " *", "")))

  # Direction-consistency tally for the Up score (from the per-cohort family).
  up <- fdr %>% filter(Family == "PerCohort", Test == "Cox", Score == "Up_ssGSEA") %>%
    mutate(Is_Testable = as.logical(Is_Testable), Is_Significant = as.logical(Is_Significant),
           HR = .pancan_num(HR))
  tal <- function(mt) { x <- up %>% filter(Metric == mt, Is_Testable, !is.na(HR))
    list(gt = sum(x$HR > 1), n = nrow(x), sig = sum(x$Is_Significant, na.rm = TRUE)) }
  os <- tal("OS"); rf <- tal("RFS")
  tally_str <- sprintf("DTP Up direction across cohorts - HR > 1: %d of %d (OS), %d of %d (RFS); FDR-significant (OS): %d",
                       os$gt, os$n, rf$gt, rf$n, os$sig)
  message("  [tally] ", tally_str)

  lab_x <- FOREST_XLIM[2]
  p <- ggplot(d, aes(x = hr_d, y = ScoreLab)) +
    geom_vline(xintercept = 1, linetype = "dashed", colour = "grey55") +
    geom_errorbarh(aes(xmin = lo_d, xmax = hi_d, colour = dir), height = 0.22,
                   linewidth = 0.6, na.rm = TRUE) +
    geom_point(aes(colour = dir, shape = offscale), size = 2.8, na.rm = TRUE) +
    # per-SD HR + C-index + raw p + FDR q next to each aggregate estimate (the pooled
    # OS effect is significant but its concordance is ~0.53 — do not oversell it).
    geom_text(aes(x = lab_x, label = sprintf("HR %.2f  C=%.2f\n%s  %s%s", hr_sd, C_index,
                                              ifelse(is.na(Raw_P), "", ifelse(Raw_P < 0.001, "p<.001", sprintf("p=%.3f", Raw_P))),
                                              ifelse(FDR_P < 0.001, "q<.001", sprintf("q=%.3f", FDR_P)),
                                              .sig_stars(FDR_P, ""))),
              hjust = 1, vjust = -0.6, size = 2.4, colour = "grey20", lineheight = 0.82) +
    facet_grid(. ~ Metric, labeller = labeller(Metric = .PANCAN_ENDPT_LAB)) +
    scale_x_log10(limits = FOREST_XLIM, breaks = c(0.1, 0.3, 1, 3, 10),
                  labels = c("0.1", "0.3", "1", "3", "10")) +
    scale_colour_manual(values = .PANCAN_DIR_PAL, name = "Score effect", drop = FALSE) +
    scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 21), guide = "none") +
    labs(title = "Pan-cancer aggregate (cancer type as batch)",
         subtitle = "Pooled continuous-score Cox HR per 1 SD of the batch-corrected score (95% CI). C = concordance; cells add raw p, FDR q.",
         caption = tally_str,
         x = "Hazard ratio per 1 SD (log scale)", y = NULL) +
    .base_theme +
    theme(plot.caption = element_text(hjust = 0, size = 9, colour = "grey25"),
          panel.grid.major.x = element_line(colour = "grey93", linewidth = 0.3),
          panel.spacing.x = unit(9, "pt"))

  .save_fig(p, "Fig_Pancan_aggregate_forest", 8.5, 3.2, out_dir)
  p
}

# ==============================================================================
# TABLE — Table_3_5_pancan_survival.csv (Appendix Table A4)
# ==============================================================================
build_pancan_table <- function(fdr, master, out_dir) {
  .dplyr_local(environment())
  rd3 <- function(x) round(x, 3); rd4 <- function(x) round(x, 4)

  # Per-cohort rows (all three scores, both endpoints), per-SD rescaled.
  per <- fdr %>%
    filter(Family == "PerCohort", Test == "Cox", Metric %in% c("OS", "RFS")) %>%
    mutate(Is_Testable = as.logical(Is_Testable))
  sd_vec <- mapply(function(pr, sc, mt) .pancan_sd(master, pr, sc, mt),
                   per$Project, per$Score, per$Metric)
  per <- .pancan_rescale(per, sd_vec)

  # Aggregate rows (batch-corrected), per-SD rescaled by the corrected pooled SD.
  raw_cols <- c("Up_ssGSEA", "Down_ssGSEA", "Composite_ssGSEA")
  sm   <- t(as.matrix(sapply(master[, raw_cols], .pancan_num)))
  corr <- as.data.frame(t(limma::removeBatchEffect(sm, batch = master$Project_ID)))
  names(corr) <- paste0(raw_cols, "_Corrected")
  for (mc in c("OS3Y_event", "OS3Y_delay", "RFS3Y_event", "RFS3Y_delay"))
    corr[[mc]] <- .pancan_num(master[[mc]])
  agg_sd <- function(score, metric) {
    ev <- .PANCAN_METRIC_EVT[[metric]]; tm <- .PANCAN_METRIC_TIM[[metric]]
    v <- corr[[score]]; ok <- !is.na(v) & !is.na(corr[[tm]]) & !is.na(corr[[ev]])
    if (sum(ok) < 2) NA_real_ else stats::sd(v[ok])
  }
  agg <- fdr %>% filter(Family == "PanCancer", Test == "Cox", Metric %in% c("OS", "RFS")) %>%
    mutate(Is_Testable = as.logical(Is_Testable))
  agg <- .pancan_rescale(agg, mapply(function(sc, mt) agg_sd(sc, mt), agg$Score, agg$Metric))

  mk_rows <- function(df, project_label) df %>%
    transmute(Project = if (is.null(project_label)) Project else project_label,
              Score, Endpoint = Metric, N = as.integer(.pancan_num(N)),
              N_events = as.integer(.pancan_num(N_events)), Is_Testable,
              HR_perSD = rd3(hr_sd), HR_lower = rd3(lo_sd), HR_upper = rd3(hi_sd),
              C_index = rd3(.pancan_num(C_index)), Cox_FDR_P = rd4(.pancan_num(FDR_P)),
              Is_Significant = as.logical(Is_Significant))

  tbl <- bind_rows(
    mk_rows(per, NULL) %>% arrange(Project, factor(Score, paste0(c("Up","Down","Composite"), "_ssGSEA")),
                                   factor(Endpoint, c("OS","RFS"))),
    mk_rows(agg, "Pan-Cancer (batch-corrected)") %>%
      arrange(factor(Score, paste0(c("Up","Down","Composite"), "_ssGSEA_Corrected")),
              factor(Endpoint, c("OS","RFS"))))
  readr::write_csv(tbl, file.path(out_dir, "Table_3_5_pancan_survival.csv"))
  message(sprintf("  wrote Table_3_5_pancan_survival.csv (%d rows: %d per-cohort + %d aggregate)",
                  nrow(tbl), nrow(per), nrow(agg)))
  invisible(tbl)
}

# ==============================================================================
# FIGURE 7 — pooled pan-cancer association as violins. Batch-corrected score
#            (cancer type as batch, recomputed from the master exactly as the
#            aggregate) split by 3-year outcome group, for the three scores x two
#            endpoints. Each panel annotated with the Pan-Cancer Wilcoxon
#            rank-biserial r + BH-FDR (negative r = higher score in the event group).
# ==============================================================================
build_pancan_violins <- function(fdr, master, out_dir, save_name = "Fig7_pancan_violins") {
  .dplyr_local(environment())
  SCORES <- c(Up_ssGSEA_Corrected = "DTP Up", Down_ssGSEA_Corrected = "DTP Down",
              Composite_ssGSEA_Corrected = "DTP Composite")
  SCORE_LEVELS <- c("DTP Up", "DTP Down", "DTP Composite")
  PAL <- c("No event" = "#4DBBD5", "Event" = "#E64B35")

  # Recompute the batch-corrected scores (removeBatchEffect, cancer type as batch),
  # identical to the aggregate; then attach the 3-year outcome-group labels.
  n_coh <- length(unique(master$Project_ID))
  raw_cols <- c("Up_ssGSEA", "Down_ssGSEA", "Composite_ssGSEA")
  sm   <- t(as.matrix(sapply(master[, raw_cols], .pancan_num)))
  corr <- as.data.frame(t(limma::removeBatchEffect(sm, batch = master$Project_ID)))
  names(corr) <- paste0(raw_cols, "_Corrected")
  corr$Surv_3yr       <- master$Surv_3yr
  corr$Recurrence_3yr <- master$Recurrence_3yr

  # Long frame: one row per patient x score x endpoint, with the outcome group.
  mk_long <- function(score_col, score_lab) dplyr::bind_rows(
    data.frame(ScoreLab = score_lab, Endpoint = "OS", value = corr[[score_col]],
               group = c(Dead_3yr = "Event", Alive_3yr = "No event")[corr$Surv_3yr],
               stringsAsFactors = FALSE),
    data.frame(ScoreLab = score_lab, Endpoint = "RFS", value = corr[[score_col]],
               group = c(Recurred = "Event", `Recurrence-Free` = "No event")[corr$Recurrence_3yr],
               stringsAsFactors = FALSE))
  long <- dplyr::bind_rows(Map(mk_long, names(SCORES), unname(SCORES))) %>%
    filter(!is.na(group), !is.na(value)) %>%
    mutate(ScoreLab = factor(ScoreLab, SCORE_LEVELS),
           Endpoint = factor(Endpoint, c("OS", "RFS")),
           group    = factor(group, c("No event", "Event")))

  # Per-panel Wilcoxon annotation (Pan-Cancer family) + per-group n. Clean label =
  # rank-biserial r + FDR; full label additionally shows the raw Wilcoxon p.
  w <- fdr %>% filter(Project == "Pan-Cancer", Test == "Wilcoxon", Metric %in% c("OS", "RFS")) %>%
    mutate(Effect_r = .pancan_num(Effect_r), FDR_P = .pancan_num(FDR_P), Raw_P = .pancan_num(Raw_P),
           ScoreLab = factor(SCORES[Score], SCORE_LEVELS),
           Endpoint = factor(Metric, c("OS", "RFS")),
           lab_clean = sprintf("r = %+.3f\nFDR %s%s", Effect_r,
                               ifelse(FDR_P < 0.001, "< 0.001", sprintf("= %.3f", FDR_P)),
                               .sig_stars(FDR_P)),
           lab_full = sprintf("r = %+.3f\n%s  FDR %s%s", Effect_r,
                              ifelse(is.na(Raw_P), "", ifelse(Raw_P < 0.001, "p<.001", sprintf("p=%.3f", Raw_P))),
                              ifelse(FDR_P < 0.001, "< 0.001", sprintf("= %.3f", FDR_P)),
                              .sig_stars(FDR_P)))
  nlab <- long %>% group_by(ScoreLab, Endpoint, group) %>%
    summarise(n = dplyr::n(), .groups = "drop")

  message("  [violins] Pan-Cancer Wilcoxon r / FDR per score x endpoint:")
  for (i in order(w$ScoreLab, w$Endpoint))
    message(sprintf("    %-14s %-3s  r=%+.3f  FDR=%.4f%s", as.character(w$ScoreLab[i]),
                    as.character(w$Endpoint[i]), w$Effect_r[i], w$FDR_P[i],
                    ifelse(w$FDR_P[i] < 0.05, " *", "")))

  make_fig <- function(detailed) {
    ww <- w; ww$lab <- if (detailed) ww$lab_full else ww$lab_clean
    ggplot(long, aes(x = group, y = value, fill = group)) +
      geom_violin(alpha = 0.85, trim = FALSE, scale = "width", linewidth = 0.3, colour = "grey30") +
      geom_boxplot(width = 0.16, fill = "white", alpha = 0.9, outlier.shape = NA, linewidth = 0.3) +
      geom_text(data = ww, inherit.aes = FALSE, aes(x = 1.5, y = Inf, label = lab),
                vjust = 1.15, size = 2.7, colour = "grey20", lineheight = 0.9) +
      geom_text(data = nlab, inherit.aes = FALSE, aes(x = group, y = -Inf, label = paste0("n=", n)),
                vjust = -0.6, size = 2.3, colour = "grey40") +
      facet_grid(ScoreLab ~ Endpoint, scales = "free_y",
                 labeller = labeller(Endpoint = .PANCAN_ENDPT_LAB)) +
      scale_fill_manual(values = PAL, name = "3-year outcome",
                        labels = c("No event" = "No event", "Event" = "Event (died / recurred)")) +
      scale_y_continuous(expand = expansion(mult = 0.15)) +
      labs(title = "Pooled pan-cancer DTP score by 3-year outcome",
           subtitle = paste0(
             sprintf("Batch-corrected score (cancer type as batch), all %d cohorts pooled. ", n_coh),
             if (detailed) "Wilcoxon rank-biserial r, raw p and BH-FDR q (Pan-Cancer family).\n"
             else "Wilcoxon rank-biserial r and BH-FDR (Pan-Cancer family).\n",
             "Negative r = higher score in the event group. * FDR<0.05, ** <0.01, *** <0.001."),
           caption = "Negative rank-biserial r = higher batch-corrected score in the event (died / recurred) group.",
           x = NULL, y = "Batch-corrected ssGSEA score") +
      .base_theme +
      theme(plot.caption = element_text(hjust = 0, size = 8.5, colour = "grey30"),
            panel.spacing = unit(7, "pt"),
            panel.grid.major.y = element_line(colour = "grey93", linewidth = 0.3))
  }

  .save_fig(make_fig(TRUE), save_name, 8.5, 9.0, out_dir)
  .save_fig(.strip_titles(make_fig(FALSE)), save_name, 8.5, 9.0, file.path(out_dir, "publication"))
  invisible(make_fig(TRUE))
}

# Assemble the Fig 6 composite from three build_pancan_forest CLOSURES: A = Up on
# top (full width), B = Composite and C = Down side by side beneath it (identical
# row order). Only Panel A keeps the shared legend. Flat design string so
# tag_levels = "A" numbers the patches in add-order. The titled composite uses
# detailed = TRUE; the publication copy uses the clean detailed = FALSE panels
# (title-stripped BEFORE wrap_elements). Shared by the all-patient and treated
# sections so the two figures stay laid out identically.
.pancan_fig6_composite <- function(pA, pB, pC, out_dir,
                                   save_name = "Fig6_pancan_survival_composite") {
  if (is.null(pA) || is.null(pB) || is.null(pC))
    stop("panel A/B/C did not all build; composite skipped")
  nl     <- theme(legend.position = "none")
  design <- "AA\nBC"
  mk <- function(a, b, cc) wrap_elements(a) + wrap_elements(b + nl) + wrap_elements(cc + nl) +
    plot_layout(design = design, heights = c(1, 1)) + plot_annotation(tag_levels = "A") &
    theme(plot.tag = element_text(face = "bold", size = 16))
  .save_fig(mk(pA(TRUE), pB(TRUE), pC(TRUE)), save_name, 15, 16.5, out_dir)
  .save_fig(mk(.strip_titles(pA(FALSE)), .strip_titles(pB(FALSE)), .strip_titles(pC(FALSE))),
            save_name, 15, 16.5, file.path(out_dir, "publication"))
  invisible(TRUE)
}

# ==============================================================================
# Orchestrator — mirrors build_crc_composites (guarded panels + composite)
# ==============================================================================
build_pancan_composites <- function(run_dir, out_dir = file.path(run_dir, "composites")) {
  stopifnot(dir.exists(run_dir))
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  message("Building pan-cancer composite figures\n  run: ", run_dir, "\n  out: ", out_dir)

  fdr    <- .rd(run_dir, "FDR_Stats_Summary.csv")
  master <- .rd(run_dir, "MASTER_Pan_Cancer_Clinical.csv")
  n_coh  <- length(unique(master$Project_ID))   # cohort count (data-driven, not hardcoded)

  panels <- new.env(parent = emptyenv())
  groups <- list(
    "A (Up forest)"        = function() panels$A <- build_pancan_forest(
                                fdr, master, out_dir, "Up_ssGSEA", "DTP Up", "Fig6A_pancan_forest_Up"),
    "B (Composite forest)" = function() panels$B <- build_pancan_forest(
                                fdr, master, out_dir, "Composite_ssGSEA", "DTP Composite", "Fig6B_pancan_forest_Composite"),
    "C (Down forest)"      = function() panels$C <- build_pancan_forest(
                                fdr, master, out_dir, "Down_ssGSEA", "DTP Down", "Fig6C_pancan_forest_Down"),
    "aggregate (appendix)" = function() build_pancan_aggregate(fdr, master, out_dir),
    "violins (Fig7)"       = function() build_pancan_violins(fdr, master, out_dir),
    "table"                = function() build_pancan_table(fdr, master, out_dir))
  status <- setNames(logical(length(groups)), names(groups))
  for (nm in names(groups)) {
    status[[nm]] <- tryCatch({ groups[[nm]](); TRUE },
      error = function(e) {
        msg <- sprintf("[pancan composites] group %s FAILED: %s", nm, conditionMessage(e))
        message("  ", msg); warning(msg, call. = FALSE); FALSE
      })
  }

  # Fig 6 composite (guarded); layout lives in .pancan_fig6_composite so the
  # treated section builds an identically laid-out figure.
  status[["composite"]] <- tryCatch({
    .pancan_fig6_composite(panels$A, panels$B, panels$C, out_dir)
    TRUE
  }, error = function(e) {
    msg <- sprintf("[pancan composites] group composite FAILED: %s", conditionMessage(e))
    message("  ", msg); warning(msg, call. = FALSE); FALSE
  })

  # ---- run-specific numbers for the caption blocks --------------------------
  # These captions are REWRITTEN on every run and sit next to a freshly drawn
  # figure, so every quantitative claim in them has to be read out of `fdr` — a
  # hardcoded cohort list or r value would be indistinguishable from a current one
  # after the signature panel or a cohort set changes.
  .num  <- .pancan_num
  .tv   <- function(x) { b <- suppressWarnings(as.logical(x)); !is.na(b) & b }
  .cohort_endpoints <- function(sc, want_sig = TRUE) {
    r <- fdr[fdr$Family == "PerCohort" & fdr$Test == "Cox" & fdr$Score == sc &
             fdr$Metric %in% c("OS", "RFS"), , drop = FALSE]
    r <- if (want_sig) r[.tv(r$Is_Significant), , drop = FALSE]
         else          r[!.tv(r$Is_Testable), , drop = FALSE]
    if (!nrow(r)) return(character(0))
    r <- r[order(r$Project, r$Metric), , drop = FALSE]
    paste(r$Project, r$Metric)
  }
  .listify <- function(x, none = "none") if (!length(x)) none else paste(x, collapse = ", ")
  nt_up   <- .cohort_endpoints("Up_ssGSEA", want_sig = FALSE)
  sig_up  <- .cohort_endpoints("Up_ssGSEA")
  sig_cmp <- .cohort_endpoints("Composite_ssGSEA")
  sig_dn  <- .cohort_endpoints("Down_ssGSEA")
  # Pan-Cancer Wilcoxon r / FDR per score x endpoint, in the caption's own wording.
  wc <- fdr[fdr$Family == "PanCancer" & fdr$Test == "Wilcoxon" &
            fdr$Metric %in% c("OS", "RFS"), , drop = FALSE]
  wc_txt <- if (!nrow(wc)) "no Pan-Cancer Wilcoxon rows in this run" else {
    lab <- .PANCAN_CORR_LAB[wc$Score]
    lab <- ifelse(is.na(lab), wc$Score, sub("^DTP ", "", lab))
    f   <- .num(wc$FDR_P)
    .listify(sprintf("%s %s r = %+.3f (FDR %s)", lab, wc$Metric, .num(wc$Effect_r),
                     ifelse(is.na(f), "n/a",
                       ifelse(f < 0.001, "< 0.001",
                         paste0(sprintf("%.3f", f), ifelse(f >= 0.05, ", n.s.", ""))))))
  }

  # Caption block alongside the figure.
  cap <- c(
    "Figure 6. Per-cohort DTP-score survival association across the pan-cancer TCGA cohort (Section 3.5).",
    "Per-cohort forests, one per DTP score, on a shared cohort order: DTP Up on top, DTP Composite and DTP Down below.",
    "",
    "(A) DTP Up, (B) DTP Composite, (C) DTP Down. Each panel is a forest of the univariable continuous-score Cox",
    sprintf("hazard ratio per 1 standard deviation of that score within each of %d solid-tumour TCGA cohorts, with 95%%", n_coh),
    "confidence interval on a log scale (truncated to [0.1, 10]), faceted by endpoint: overall survival (OS) and",
    "recurrence (RFS, derived from the TCGA-CDR disease-free interval). Red = HR > 1 (higher score, higher",
    "hazard), blue = HR < 1; an asterisk marks FDR < 0.05; a grey \"n/t\" marks a cohort that failed the Cox gate",
    sprintf("(n >= 10 and >= 5 events; in this run: %s). All three panels share the SAME row order -- ascending DTP Up",
            .listify(nt_up, "no cohort failed the gate")),
    "overall-survival per-SD HR -- so the tissues line up; TCGA-COAD (the colorectal tissue) is highlighted in each.",
    sprintf("Up is FDR-significant in %d cohort-endpoint(s) (%s);", length(sig_up), .listify(sig_up)),
    sprintf("Composite in %d (%s);", length(sig_cmp), .listify(sig_cmp)),
    sprintf("Down in %d (%s).", length(sig_dn), .listify(sig_dn)),
    "",
    "Hazard ratios are per 1 SD; the raw Cox HRs are per unit of the (log2) ssGSEA score and are not shown. FDR is",
    "Benjamini-Hochberg within each family (per-cohort vs. batch-corrected aggregate) and test, pooling OS and RFS.",
    "The pooled association is shown as violins in Figure 7; the batch-corrected aggregate Cox forest is an appendix",
    "figure (Fig_Pancan_aggregate_forest). Full per-cohort and aggregate statistics are in Appendix Table A4",
    "(Table_3_5_pancan_survival.csv).")
  writeLines(cap, file.path(out_dir, "Fig6_pancan_survival_caption.txt"))

  # Figure 7 caption (pooled violins).
  cap7 <- c(
    "Figure 7. Pooled pan-cancer DTP-score association by 3-year outcome (Section 3.5).",
    "",
    "Batch-corrected ssGSEA scores (limma::removeBatchEffect with cancer type as batch, recomputed from the master",
    sprintf("exactly as the aggregate) for all %d solid-tumour TCGA cohorts pooled, shown as violins split by 3-year", n_coh),
    "outcome group -- no event vs event (died for OS; recurred for RFS). Six panels: the DTP Up, Down and Composite",
    "scores (rows) by overall survival and recurrence (columns). Each panel is annotated with the Pan-Cancer",
    "Wilcoxon rank-biserial effect size r and its Benjamini-Hochberg FDR; a NEGATIVE r means a higher score in the",
    paste0("event group. In this run: ", wc_txt, "."),
    "Violin widths are scaled equal so the unequal group sizes (annotated as n) are comparable.")
  writeLines(cap7, file.path(out_dir, "Fig7_pancan_violins_caption.txt"))

  n_ok <- sum(status); n <- length(status)
  summ <- sprintf("Pan-cancer composites: %d/%d groups OK%s", n_ok, n,
    if (n_ok < n) paste0(" — FAILED: ", paste(names(status)[!status], collapse = ", ")) else "")
  message(summ)
  if (n_ok < n) warning(summ, call. = FALSE)
  message("Pan-cancer composite figures done -> ", out_dir)
  invisible(status)
}

# ==============================================================================
# PART 3 — Metastasis GSEA composite (build_mets_gsea_composite)
# ==============================================================================
# ---------------------------------------------------------------------------
# .mets_nes_heatmap()  —  one-column NES heatmap (fill = NES, label = NES +
# FDR stars) over all signatures for a single GSEA contrast. Signatures are
# ordered top->bottom by `sig_order`; fill limits `lim` are shared across the
# two contrasts so the colours are directly comparable. The NES sign convention
# lives on the x-axis title (not the panel title) so it survives title-stripping
# in the publication copy. Significance uses the shared .sig_stars (plotting.R).
# ---------------------------------------------------------------------------
.mets_nes_heatmap <- function(gsea_df, contrast_lab, pos_label, sig_order, lim, detailed = FALSE) {
  d <- gsea_df
  d$Signature <- factor(d$ID, levels = rev(sig_order))
  d$Contrast  <- contrast_lab
  clean <- paste0(sprintf("%.2f", d$NES), .sig_stars(d$p.adjust))
  # Titled copy: NES + raw p + FDR q; publication copy keeps NES + stars only.
  d$lab <- if (detailed)
    mapply(.stat_annot, detailed = TRUE, clean = clean, effect = sprintf("%.2f", d$NES),
           raw_p = d$pvalue, fdr = d$p.adjust)
  else clean
  ggplot2::ggplot(d, ggplot2::aes(x = Contrast, y = Signature, fill = NES)) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.6) +
    ggplot2::geom_text(ggplot2::aes(label = lab), size = if (detailed) 2.3 else 3,
                       fontface = "bold", colour = "grey10", lineheight = 0.82) +
    ggplot2::scale_fill_gradient2(low = "#3B6EA5", mid = "white", high = "#E64B35",
      midpoint = 0, limits = c(-lim, lim), name = "NES") +
    ggplot2::labs(title = "NES + FDR",
                  x = sprintf("+NES -> up in %s", pos_label), y = NULL) +
    pub_theme +
    ggplot2::theme(legend.position = "right",
                   panel.grid  = ggplot2::element_blank(),
                   axis.ticks  = ggplot2::element_blank(),
                   axis.title.x = ggplot2::element_text(size = 9, face = "bold"),
                   plot.title  = ggplot2::element_text(size = 11, face = "bold", hjust = 0))
}

# ---------------------------------------------------------------------------
# build_mets_gsea_composite()
#
# Four-panel composite for the mets GSEA analysis (thesis mets figure):
#   A  Primary-vs-Normal                       all-signature enrichment plot
#   B  Primary-vs-Normal                       NES heatmap (NES + FDR stars)
#   C  Primary-vs-Metastasis, PURITY-ADJUSTED  all-signature enrichment plot
#   D  Primary-vs-Metastasis, PURITY-ADJUSTED  NES heatmap
#
# The whole figure is anchored on the primary tumour, matching the thesis
# narrative (primary -> metastasis): panels A/B read Primary vs Normal (positive
# NES = up in Primary) and panels C/D read Primary vs Metastasis. Panels C/D read
# the liver-score-adjusted ranking: the metastasis arm carries heavy hepatic
# signal (liver-validation figure), so the metastasis-direction enrichment must
# be read off the adjusted fit. Panels A/B (Primary-vs-Normal) have no liver
# tissue and stay unadjusted. The "purity-adjusted" distinction lives in the C
# title + the D column/axis labels so it survives title-stripping in the
# publication copy.
#
# Rebuildable from a completed run: reads the saved gseaResult objects
# (results/gsea_PrimaryVsNormal.rds + gsea_PrimaryVsMetastasis_adjusted.rds
# written by run_mets_de). NES / FDR for the heatmaps come from as.data.frame(obj).
# ---------------------------------------------------------------------------
build_mets_gsea_composite <- function(mets_dir,
                                      out_dir = file.path(mets_dir, "plots", "GSEA")) {
  pn_rds <- file.path(mets_dir, "results", "gsea_PrimaryVsNormal.rds")
  # Primary-vs-Metastasis: the purity-adjusted ranking is the primary figure.
  pm_rds <- file.path(mets_dir, "results", "gsea_PrimaryVsMetastasis_adjusted.rds")
  if (!file.exists(pn_rds) || !file.exists(pm_rds))
    stop("build_mets_gsea_composite(): missing gsea RDS in ",
         file.path(mets_dir, "results"), " (need gsea_PrimaryVsNormal.rds + ",
         "gsea_PrimaryVsMetastasis_adjusted.rds; the Normal-vs-Primary contrast was ",
         "re-oriented to Primary-vs-Normal, so run_mets_de() must be re-run against ",
         "this run root before the composite can be rebuilt).")
  gsea_pn <- readRDS(pn_rds); gsea_pm <- readRDS(pm_rds)

  df_pn <- as.data.frame(gsea_pn); df_pm <- as.data.frame(gsea_pm)
  if (!nrow(df_pn) || !nrow(df_pm))
    stop("build_mets_gsea_composite(): empty GSEA result(s).")

  # Shared signature order (top->bottom) = descending Primary-vs-Metastasis NES,
  # the headline contrast; any signature absent there is appended.
  ord       <- df_pm$ID[order(-df_pm$NES)]
  sig_order <- c(ord, setdiff(union(df_pn$ID, df_pm$ID), ord))
  glim      <- max(abs(c(df_pn$NES, df_pm$NES)), na.rm = TRUE)

  # sig_order is the UNION of the two results, so it can name a signature that one
  # object lacks (a set can drop out of one contrast and not the other). Each
  # enrichment plot may only be asked for IDs its own object holds —
  # .gsea_running_data() stop()s otherwise, which would take the whole composite AND
  # build_mets_liver_composite() (same tryCatch in run_composites) down with it. The
  # shared order is preserved; only the missing rows are omitted, and flagged.
  pn_ids <- intersect(sig_order, df_pn$ID)
  pm_ids <- intersect(sig_order, df_pm$ID)
  if (length(pn_ids) < length(sig_order) || length(pm_ids) < length(sig_order))
    message("  [mets GSEA] signature(s) absent from one contrast, omitted from that panel: ",
            paste(setdiff(sig_order, intersect(pn_ids, pm_ids)), collapse = ", "))

  # A / C enrichment plots have no per-cell stats to add, so build once. B / D NES
  # heatmaps carry the per-signature stats: rebuilt inside mk() with detailed = TRUE
  # for the titled copy, detailed = FALSE for the clean publication copy.
  A <- plot_gsea_enrichment(gsea_pn, gene_set_ids = pn_ids,
         title = "Primary vs Normal: all signatures",
         label_left = "Up in Primary", label_right = "Up in Normal")
  C <- plot_gsea_enrichment(gsea_pm, gene_set_ids = pm_ids,
         title = "Primary vs Metastasis, purity-adjusted for liver score: all signatures",
         label_left = "Up in Metastasis", label_right = "Up in Primary")

  # Assemble twice: titled original (detailed heatmaps), and a clean title-stripped
  # copy for publication/. Strip the panels (A–D) BEFORE wrap_elements, and drop the
  # overall annotation title, so nothing survives frozen inside a wrap_elements() cell.
  mk <- function(strip) {
    detailed <- !strip
    b <- .mets_nes_heatmap(df_pn, "P / N", pos_label = "Primary", sig_order, glim, detailed)
    d <- .mets_nes_heatmap(df_pm, "P / M (adj)", pos_label = "Metastasis, purity-adjusted",
                           sig_order, glim, detailed)
    a <- A; cc <- C
    if (strip) { a <- .strip_titles(A); cc <- .strip_titles(C)
                 b <- .strip_titles(b); d <- .strip_titles(d) }
    patchwork::wrap_plots(
      patchwork::wrap_elements(a), patchwork::wrap_elements(b),
      patchwork::wrap_elements(cc), patchwork::wrap_elements(d),
      ncol = 2, widths = c(3, 1)) +
    patchwork::plot_annotation(
      tag_levels = "A",
      title = if (strip) NULL else
        "Mets GSEA: primary tumour vs normal mucosa and vs liver metastasis (metastasis arm purity-adjusted for liver score)",
      theme = ggplot2::theme(plot.title = if (strip) ggplot2::element_blank() else
        ggplot2::element_text(face = "bold", size = 15))) &
      ggplot2::theme(plot.tag = ggplot2::element_text(face = "bold", size = 16))
  }

  n_sig <- length(sig_order); h <- max(14, 2 * (n_sig * 0.55 + 4))
  .save_fig(mk(FALSE), "Fig_Mets_GSEA_composite", 16, h, out_dir)
  .save_fig(mk(TRUE),  "Fig_Mets_GSEA_composite", 16, h, file.path(out_dir, "publication"))
  message("  (", n_sig, " signatures, ", round(h, 1), " in tall)")
  invisible(mk(FALSE))
}

# ==============================================================================
# METS COMPOSITE #2 — liver purity in GSE50760 metastases and its effect on the
#            Primary-vs-Metastasis GSEA (thesis mets figure Fig_Mets_LiverPurity).
#            (A) per-patient liver-score trajectory Normal -> Primary -> Metastasis
#            (metastases carry heavy hepatic signal); (B) unadjusted-vs-adjusted
#            Primary-vs-Metastasis NES heatmap. Rebuildable from the run's
#            PurityAdjusted outputs; guarded on file existence.
# ==============================================================================
build_mets_liver_composite <- function(mets_dir,
                                       out_dir = file.path(mets_dir, "plots", "GSEA")) {
  liver_csv <- file.path(mets_dir, "results", "PurityAdjusted", "liver_validation_scores_3way.csv")
  cmp_csv   <- file.path(mets_dir, "results", "GSEA_PvM_adjusted_vs_unadjusted.csv")
  if (!file.exists(liver_csv) || !file.exists(cmp_csv)) {
    message("[mets liver composite] skip (need ", basename(liver_csv), " + ",
            basename(cmp_csv), " under ", file.path(mets_dir, "results"), ")")
    return(invisible(NULL))
  }
  .dplyr_local(environment())

  # -- Panel A: per-patient liver-score trajectory (the purity problem) --------
  liver_df <- read.csv(liver_csv, stringsAsFactors = FALSE)
  liver_df$Tissue <- factor(liver_df$Tissue, levels = c("Normal", "Primary", "Metastasis"))
  pA <- ggplot(liver_df, aes(Tissue, Liver, group = Patient)) +
    geom_line(colour = "gray60", alpha = 0.6) +
    geom_point(aes(colour = Tissue), size = 3) +
    scale_colour_manual(values = pub_palette) +
    labs(title = "Liver signal across Normal -> Primary -> Metastasis",
         subtitle = "HSIAO_LIVER_SPECIFIC_GENES ssGSEA per sample; lines join each patient's three tissues.",
         x = NULL, y = "Liver score") +
    .base_theme + theme(legend.position = "none")

  # -- Panel B: unadjusted vs purity-adjusted P/M NES (what adjustment changes) -
  # This comparison CSV carries NES + FDR only (no raw p), so the titled copy shows
  # NES + FDR q; the publication copy keeps NES + stars.
  cmp <- read.csv(cmp_csv, stringsAsFactors = FALSE)
  hm <- bind_rows(
      transmute(cmp, ID, Model = "Unadjusted", NES = NES_unadjusted, FDR = FDR_unadjusted),
      transmute(cmp, ID, Model = "Adjusted",   NES = NES_adjusted,   FDR = FDR_adjusted)) %>%
    mutate(Model = factor(Model, c("Unadjusted", "Adjusted")),
           Signature = factor(ID, levels = rev(cmp$ID[order(cmp$NES_adjusted)])),
           lab_clean = paste0(sprintf("%.2f", NES), .sig_stars(FDR)),
           lab_full  = mapply(.stat_annot, detailed = TRUE, clean = lab_clean,
                              effect = sprintf("%.2f", NES), fdr = FDR))
  lim <- max(abs(hm$NES), na.rm = TRUE)
  make_B <- function(detailed) {
    dd <- hm; dd$lab <- if (detailed) dd$lab_full else dd$lab_clean
    ggplot(dd, aes(Model, Signature, fill = NES)) +
      geom_tile(colour = "white", linewidth = 0.6) +
      geom_text(aes(label = lab), size = if (detailed) 2.5 else 3, fontface = "bold",
                colour = "grey10", lineheight = 0.82) +
      scale_fill_gradient2(low = "#3B6EA5", mid = "white", high = "#E64B35",
                           midpoint = 0, limits = c(-lim, lim), name = "NES") +
      labs(title = "Primary vs Metastasis NES: unadjusted vs purity-adjusted",
           subtitle = if (detailed) "NES with BH-FDR q per cell; * <0.05, ** <0.01, *** <0.001."
                      else "FDR stars * <0.05, ** <0.01, *** <0.001 (BH).",
           x = "Primary vs Metastasis ranking", y = NULL) +
      .base_theme +
      theme(panel.grid = element_blank(), axis.ticks = element_blank(), legend.position = "right")
  }

  CAP <- paste0(
    "(A) HSIAO_LIVER_SPECIFIC_GENES ssGSEA score per sample; lines join each patient's three tissues - metastases carry heavy hepatic signal.\n",
    "(B) Primary-vs-Metastasis GSEA NES before vs after adjusting the ranking for the liver score; FDR stars * <0.05, ** <0.01, *** <0.001 (BH).")

  mk <- function(strip) {
    a <- if (strip) .strip_titles(pA) else pA
    b <- make_B(!strip); if (strip) b <- .strip_titles(b)
    (wrap_elements(a) + wrap_elements(b)) +
      plot_layout(widths = c(1.25, 1)) +
      plot_annotation(
        tag_levels = "A",
        title   = if (strip) NULL else "Liver purity in GSE50760 metastases and its effect on Primary-vs-Metastasis GSEA",
        caption = CAP,
        theme = theme(
          plot.title   = if (strip) element_blank() else element_text(face = "bold", size = 14, hjust = 0),
          plot.caption = element_text(hjust = 0, size = 8.5, colour = "grey30"))) &
      theme(plot.tag = element_text(face = "bold", size = 16))
  }

  .save_fig(mk(FALSE), "Fig_Mets_LiverPurity_composite", 12.5, 6.2, out_dir)
  .save_fig(mk(TRUE),  "Fig_Mets_LiverPurity_composite", 12.5, 6.2, file.path(out_dir, "publication"))
  invisible(mk(FALSE))
}

# ==============================================================================
# run_composites() — standalone entry point. Detect the newest CRC_DTP_* run root
# (or take one), then regenerate every composite for whichever analysis modules
# produced output there. Each module is guarded so one failure never sinks the rest.
# ==============================================================================
run_composites <- function(out_root = NULL) {
  if (is.null(out_root)) {
    roots <- sort(Sys.glob("CRC_DTP_*"), decreasing = TRUE)
    roots <- roots[dir.exists(roots)]
    if (!length(roots)) stop("run_composites(): no CRC_DTP_* run root found; pass out_root.")
    out_root <- roots[[1]]
  }
  stopifnot(dir.exists(out_root))
  message("Composites for run root: ", out_root)
  builders <- list(
    crc_survival    = function(d) build_crc_composites(d),
    pancan_survival = function(d) build_pancan_composites(d),
    mets_de         = function(d) { build_mets_gsea_composite(d); build_mets_liver_composite(d) })
  for (m in names(builders)) {
    d <- file.path(out_root, m)
    if (dir.exists(d)) {
      tryCatch(builders[[m]](d),
        error = function(e) {
          msg <- sprintf("[composites] %s FAILED: %s", m, conditionMessage(e))
          message("  ", msg); warning(msg, call. = FALSE)
        })
    } else {
      message("[composites] skip ", m, " (no ", d, ")")
    }
  }
  message("Composites done for ", out_root)
  invisible(out_root)
}

# Auto-run only when executed as a script (skipped when source()'d by main.R).
if (sys.nframe() == 0L) {
  suppressPackageStartupMessages({
    library(dplyr); library(ggplot2); library(patchwork); library(limma)
  })
  # core/stats.R supplies harmonize_crc_modifiers(), the single definition of the
  # CMS / PDS / Stage_bin / MSI_group modifier columns group 3 rebuilds.
  for (f in c("core/config.R", "core/stats.R", "core/plotting.R", "core/gsea_plots.R")) source(f)
  .args <- commandArgs(trailingOnly = TRUE)
  run_composites(if (length(.args) >= 1) .args[[1]] else NULL)
}
