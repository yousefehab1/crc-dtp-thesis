#!/usr/bin/env Rscript
# ==============================================================================
# core/composite_figures.R
# ------------------------------------------------------------------------------
# Publication composite figures for the CRC-survival part of the thesis, built
# from a completed run's tabular outputs (no ssGSEA re-computation). Sourced by
# main.R alongside the other core files; modules/crc_survival.R calls
# build_crc_composites() at the end of its run. Three figure groups:
#
#   Group 1  build_outcome_figures()   ssGSEA score vs 3-yr OS/RFS outcome
#            Fig1A effect-size (Wilcoxon rank-biserial) matrix + Fig1B violins
#   Group 2  build_km_figures()        Kaplan-Meier by median High/Low score
#            Fig2A continuous-Cox HR/SD matrix + Fig2B KM curve grid (log-rank p)
#   Group 3  build_subgroup_figures()  effect modification by CMS/PDS/Stage/MSI
#            Fig3A interaction-FDR matrix + Fig3B per-score subgroup forests
#   Group 4  build_subtype_survival_figures()  §3.3 univariable DTP-to-outcome
#            WITHIN molecular/clinical subgroups: Fig3_3A per-1-SD Cox-HR tile
#            grid (gated cells marked not-testable) + Fig3_3B score-across-subtype
#            (Kruskal-Wallis eps^2), plus Table_3_3_subtype_survival.csv
#
# Statistics follow the main analysis (report §"continuous score for inference,
# split for display"): Wilcoxon + signed rank-biserial r for outcomes; univariable
# continuous-score Cox HR (rescaled per 1 SD) for survival/subgroups; log-rank FDR
# on the median-split KM curves; interaction-LRT FDR for effect modification.
#
# Reads from <RUN_DIR>:  *_clinical.csv, CRC_Statistical_Summary.csv,
#   Subgroup_Score_HRs.csv, Interaction_Cox_Summary.csv, Subtype_Score_Stats.csv.
#   Writes PDF + 300-dpi PNG to <RUN_DIR>/composites/.
#
# Use:
#   * in the pipeline:  sourced by main.R; run_crc_survival() calls it automatically
#   * standalone:       Rscript core/composite_figures.R [RUN_DIR]
#                       (RUN_DIR defaults to newest CRC_DTP_*/crc_survival)
# ==============================================================================

suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr); library(readr)
  library(stringr); library(forcats); library(patchwork); library(scales)
  library(survival)
})

# ---- Shared constants --------------------------------------------------------
MIN_KM_GROUP_N <- 5      # per High/Low group (matches core/config.R)
MIN_EVENTS     <- 5
X_CAP          <- 36     # months (3-yr landmark truncation)
FOREST_XLIM    <- c(0.1, 10)   # subgroup HR display window (log scale)

# Signature column -> publication label. Order = row order (inverse-signal
# modules DTP Down / MYC / Columnar grouped last). Shared by groups 1 & 2.
MODULES <- c(
  Up_ssGSEA          = "DTP Up",
  Composite_ssGSEA   = "DTP Composite",
  Foetal_ssGSEA      = "Foetal intestinal",
  regStemCell_ssGSEA = "Regenerative stem cell",
  revCSC_ssGSEA      = "Revival revCSC",
  IBD_ssGSEA         = "IBD",
  Down_ssGSEA        = "DTP Down",
  Myc_ssGSEA         = "MYC module",
  CSC_ssGSEA         = "Columnar stem cell"
)
MODULE_LEVELS <- unname(MODULES)

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
.rd <- function(run_dir, f) read.csv(file.path(run_dir, f), check.names = FALSE,
                                     stringsAsFactors = FALSE)
# stars from a p-value; `ns` is the label for testable-but-not-significant.
.sig_stars <- function(p, ns = "") ifelse(is.na(p), "",
  ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", ns))))

.base_theme <- theme_classic(base_size = 11) +
  theme(
    plot.title      = element_text(face = "bold", size = 14, hjust = 0),
    plot.subtitle   = element_text(size = 10, colour = "grey30", hjust = 0),
    strip.background = element_rect(fill = "grey92", colour = NA),
    strip.text      = element_text(face = "bold", size = 9),
    legend.position = "bottom",
    plot.tag        = element_text(face = "bold", size = 16))

.png_dev <- if (requireNamespace("ragg", quietly = TRUE)) ragg::agg_png else NULL
.save_fig <- function(p, name, w, h, out_dir) {
  # Vector PDF via base device (no cairo/X11 dependency) + 300-dpi PNG.
  ggsave(file.path(out_dir, paste0(name, ".pdf")), p, width = w, height = h,
         device = pdf, limitsize = FALSE)
  ggsave(file.path(out_dir, paste0(name, ".png")), p, width = w, height = h,
         dpi = 300, bg = "white", limitsize = FALSE, device = .png_dev)
  message("  wrote ", name, " (", w, "x", h, " in)")
}

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
  END_LABEL <- c(OS = "3-yr survival", RFS = "3-yr recurrence")
  CONTRASTS$col_lab <- paste0(CONTRASTS$cohort, "\n(", END_LABEL[CONTRASTS$endpoint], ")")
  col_lab_map <- setNames(CONTRASTS$col_lab, CONTRASTS$key)

  clin  <- list(gse = .rd(run_dir, "GSE39582_clinical.csv"),
                tcga = .rd(run_dir, "TCGA_COAD_clinical.csv"))
  stats <- .rd(run_dir, "CRC_Statistical_Summary.csv")

  cohort_subset <- function(row) {
    d <- clin[[row$clin]]
    if (row$cohort == "Adjuvant-treated") d <- d[d$Chemo_adj == "Y" & !is.na(d$Chemo_adj), ]
    d
  }
  stat_lookup <- function(row, score) {
    s <- stats[stats$Dataset == row$stats_ds & stats$Project == row$stats_cohort &
               stats$Test == "Wilcoxon" & stats$Metric == row$endpoint &
               stats$Score == score, ]
    if (!nrow(s)) return(list(fdr = NA, r = NA, n = NA, sig = FALSE))
    list(fdr = s$FDR_P[1], r = s$Effect_r[1], n = s$N[1],
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
    for (sc in names(MODULES)) {
      if (!sc %in% names(d)) next
      v  <- suppressWarnings(as.numeric(d[[sc]])); ok <- !is.na(v)
      long_rows[[length(long_rows) + 1]] <- tibble(
        key = row$key, module = MODULES[[sc]], score_col = sc, group = grp[ok], value = v[ok])
      st <- stat_lookup(row, sc)
      # signed rank-biserial r: pipeline sign is vs first level (poor), so flip
      # so + = higher in poor outcome for the red-poor palette.
      cell_rows[[length(cell_rows) + 1]] <- tibble(
        key = row$key, module = MODULES[[sc]], score_col = sc,
        fdr = st$fdr, signed_r = -st$r, sig = st$sig, stars = .sig_stars(st$fdr, "ns"))
    }
  }
  long  <- bind_rows(long_rows) %>%
    mutate(module = factor(module, MODULE_LEVELS), key = factor(key, levels(CONTRASTS$key)))
  cells <- bind_rows(cell_rows) %>%
    mutate(module = factor(module, MODULE_LEVELS), key = factor(key, levels(CONTRASTS$key))) %>%
    left_join(select(CONTRASTS, key, dataset), by = "key")

  # -- Fig 1A: effect-size matrix --
  lim <- max(abs(cells$signed_r), na.rm = TRUE)
  fig_matrix <- ggplot(cells, aes(x = key, y = fct_rev(module), fill = signed_r)) +
    geom_tile(colour = "white", linewidth = 0.6) +
    geom_text(aes(label = ifelse(sig, stars, "")), size = 4.2, fontface = "bold",
              colour = "grey15", vjust = 0.78) +
    facet_grid(. ~ dataset, scales = "free_x", space = "free_x") +
    scale_fill_gradient2(low = "#3B6EA5", mid = "white", high = "#E64B35", midpoint = 0,
      limits = c(-lim, lim), name = "Rank-biserial r",
      guide = guide_colourbar(barwidth = 9, barheight = 0.5, title.vjust = 1)) +
    scale_x_discrete(labels = col_lab_map) +
    labs(title = "DTP-module scores vs 3-year outcome across cohorts",
         subtitle = "Colour = effect size (red: higher in poor-outcome group).  * marks FDR-significant Wilcoxon (BH, p<0.05).",
         x = NULL, y = NULL) +
    .base_theme + theme(axis.text.x = element_text(size = 8, lineheight = 0.85),
                        panel.spacing.x = unit(6, "pt"))

  # -- Fig 1B: violin grid --
  if (identical(violin_modules, "all")) {
    curated_levels <- MODULE_LEVELS
  } else {
    sig_cols <- cells %>% filter(sig) %>% pull(score_col) %>% unique()
    keep <- MODULES[names(MODULES) %in% union(sig_cols, c("Up_ssGSEA", "Composite_ssGSEA"))]
    curated_levels <- MODULE_LEVELS[MODULE_LEVELS %in% unname(keep)]
  }
  vl <- long  %>% filter(module %in% curated_levels) %>% mutate(module = factor(module, curated_levels))
  vc <- cells %>% filter(module %in% curated_levels) %>% mutate(module = factor(module, curated_levels))
  yr <- vl %>% group_by(module, key) %>%
    summarise(ymin = min(value), ymax = max(value), .groups = "drop") %>% mutate(rng = ymax - ymin)
  vc <- vc %>% left_join(yr, by = c("module", "key"))
  nlab <- vl %>% group_by(module, key, group) %>% summarise(n = n(), .groups = "drop") %>%
    left_join(select(yr, module, key, ymin, rng), by = c("module", "key"))

  build_violin_block <- function(keys, title, show_y) {
    dl <- filter(vl, key %in% keys); dc <- filter(vc, key %in% keys); dn <- filter(nlab, key %in% keys)
    ggplot(dl, aes(x = group, y = value, fill = group)) +
      geom_violin(alpha = 0.85, trim = FALSE, linewidth = 0.3, colour = "grey30") +
      geom_boxplot(width = 0.16, fill = "white", alpha = 0.9, outlier.shape = NA, linewidth = 0.3) +
      geom_text(data = dc, inherit.aes = FALSE,
                aes(x = 1.5, y = ymax + 0.12 * rng, label = ifelse(stars == "", "ns", stars)),
                fontface = "bold", size = 3.2) +
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
  block_m <- build_violin_block(c("M_whole_OS","M_whole_RFS","M_trt_OS","M_trt_RFS"), "Marisa (GSE39582)", TRUE)
  block_t <- build_violin_block(c("T_whole_OS","T_whole_RFS"), "TCGA-COAD", FALSE)
  fig_violin <- block_m + block_t +
    plot_layout(widths = c(4, 2), guides = "collect") +
    plot_annotation(
      title = if (identical(violin_modules, "all")) "Signature-score distributions by 3-year outcome"
              else "Score distributions for outcome-associated modules",
      subtitle = "Full cohort, true n shown; y-axis free within each dataset (ssGSEA scale differs GSE vs TCGA).  * = FDR-significant; ns = not.",
      theme = theme(plot.title = element_text(face = "bold", size = 14, hjust = 0),
                    plot.subtitle = element_text(size = 10, colour = "grey30", hjust = 0))) &
    theme(legend.position = "bottom")

  n_rows <- length(curated_levels); viol_h <- 1.5 * n_rows + 1.4
  .save_fig(fig_matrix, "Fig1A_outcome_effect_matrix", 9.5, 4.6, out_dir)
  .save_fig(fig_violin, "Fig1B_outcome_violins", 11.0, viol_h, out_dir)
  composite <- wrap_elements(fig_matrix) / wrap_elements(fig_violin) +
    plot_layout(heights = c(4.6, viol_h - 0.4)) + plot_annotation(tag_levels = "A") &
    theme(plot.tag = element_text(face = "bold", size = 16))
  .save_fig(composite, "Fig1_outcome_composite", 11.5, 4.6 + viol_h + 0.6, out_dir)
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
  END_LABEL <- c(OS = "3-yr survival", RFS = "3-yr recurrence")
  CONTRASTS$col_lab <- paste0(CONTRASTS$cohort, "\n(", END_LABEL[CONTRASTS$metric], ")")
  col_lab_map <- setNames(CONTRASTS$col_lab, CONTRASTS$key)

  metric_cols <- function(m) if (m == "OS") c(t = "OS3Y_delay", e = "OS3Y_event")
                             else            c(t = "RFS3Y_delay", e = "RFS3Y_event")
  clin  <- list(gse = .rd(run_dir, "GSE39582_clinical.csv"),
                tcga = .rd(run_dir, "TCGA_COAD_clinical.csv"))
  stats <- .rd(run_dir, "CRC_Statistical_Summary.csv")

  cohort_subset <- function(row) {
    d <- clin[[row$clin]]
    if (row$cohort == "Adjuvant-treated") d <- d[d$Chemo_adj == "Y" & !is.na(d$Chemo_adj), ]
    d
  }
  cox_stat <- function(row, score) {   # continuous-score Cox: inferential effect
    s <- stats[stats$Dataset == row$stats_ds & stats$Project == row$stats_cohort &
               stats$Test == "Cox" & stats$Metric == row$metric & stats$Score == score, ]
    if (!nrow(s) || !isTRUE(as.logical(s$Is_Testable[1])))
      return(list(hr_unit = NA_real_, fdr = NA_real_, sig = FALSE, testable = FALSE))
    list(hr_unit = as.numeric(s$HR[1]), fdr = as.numeric(s$FDR_P[1]),
         sig = isTRUE(as.logical(s$Is_Significant[1])), testable = TRUE)
  }
  logrank_stat <- function(row, score) {   # log-rank: shown on the curves
    s <- stats[stats$Dataset == row$stats_ds & stats$Project == row$stats_cohort &
               stats$Test == "KM" & stats$Metric == row$metric & stats$Score == score, ]
    if (!nrow(s)) return(list(fdr = NA_real_, sig = FALSE))
    list(fdr = as.numeric(s$FDR_P[1]), sig = isTRUE(as.logical(s$Is_Significant[1])))
  }

  curve_rows <- list(); cell_rows <- list()
  for (i in seq_len(nrow(CONTRASTS))) {
    row <- CONTRASTS[i, ]; d0 <- cohort_subset(row)
    mc <- metric_cols(row$metric); tcol <- mc[["t"]]; ecol <- mc[["e"]]
    for (sc in names(MODULES)) {
      if (!all(c(sc, tcol, ecol) %in% names(d0))) next
      d <- d0[!is.na(d0[[sc]]) & !is.na(d0[[tcol]]) & !is.na(d0[[ecol]]), ]
      score <- as.numeric(d[[sc]])
      d <- data.frame(time = pmin(as.numeric(d[[tcol]]), X_CAP), event = as.numeric(d[[ecol]]),
                      Group = factor(ifelse(score > median(score), "High", "Low"), levels = c("Low", "High")))
      tb <- table(d$Group)
      drawable <- length(tb) == 2 && all(tb >= MIN_KM_GROUP_N) && sum(d$event) >= MIN_EVENTS
      cx <- cox_stat(row, sc); lr <- logrank_stat(row, sc)
      hr_sd <- if (cx$testable && is.finite(cx$hr_unit)) exp(log(cx$hr_unit) * stats::sd(score)) else NA_real_
      cell_rows[[length(cell_rows) + 1]] <- tibble(
        key = row$key, module = MODULES[[sc]], score_col = sc, testable = cx$testable,
        hr = hr_sd, log2hr = log2(hr_sd), fdr = cx$fdr, sig = cx$sig,
        stars = .sig_stars(cx$fdr, "ns"), lr_fdr = lr$fdr, lr_sig = lr$sig)
      if (!drawable) next
      sf <- survfit(Surv(time, event) ~ Group, data = d)
      strata <- rep(sub(".*=", "", names(sf$strata)), sf$strata)
      cd <- rbind(data.frame(time = 0, surv = 1, ncens = 0, Group = levels(d$Group)),
                  data.frame(time = sf$time, surv = sf$surv, ncens = sf$n.censor, Group = strata))
      cd$key <- row$key; cd$module <- MODULES[[sc]]; cd$small <- nrow(d) < 500
      curve_rows[[length(curve_rows) + 1]] <- cd
    }
  }
  cells <- bind_rows(cell_rows) %>%
    mutate(module = factor(module, MODULE_LEVELS), key = factor(key, levels(CONTRASTS$key))) %>%
    left_join(select(CONTRASTS, key, dataset), by = "key")
  curves <- bind_rows(curve_rows) %>%
    mutate(module = factor(module, MODULE_LEVELS), key = factor(key, levels(CONTRASTS$key)),
           Group = factor(Group, levels = c("Low", "High")))

  # -- Fig 2A: Cox HR/SD matrix --
  lim <- min(max(abs(cells$log2hr[is.finite(cells$log2hr)]), na.rm = TRUE), 2)
  fig_matrix <- ggplot(cells, aes(x = key, y = fct_rev(module), fill = pmax(pmin(log2hr, lim), -lim))) +
    geom_tile(colour = "white", linewidth = 0.6) +
    geom_text(aes(label = ifelse(sig, stars, "")), size = 4.2, fontface = "bold",
              colour = "grey15", vjust = 0.78) +
    facet_grid(. ~ dataset, scales = "free_x", space = "free_x") +
    scale_fill_gradient2(low = "#3B6EA5", mid = "white", high = "#E64B35", midpoint = 0,
      limits = c(-lim, lim), name = expression(log[2]~"HR (per 1 SD of score)"),
      guide = guide_colourbar(barwidth = 9, barheight = 0.5, title.vjust = 1), na.value = "grey85") +
    scale_x_discrete(labels = col_lab_map) +
    labs(title = "Survival association of DTP-module scores across cohorts",
         subtitle = "Univariable continuous-score Cox PH, landmark 36 mo.  Colour = log2 HR per 1 SD (red = higher score worse).  * = FDR<0.05 (BH).",
         x = NULL, y = NULL) +
    .base_theme + theme(axis.text.x = element_text(size = 8, lineheight = 0.85),
                        panel.spacing.x = unit(6, "pt"))

  # -- Fig 2B: KM curve grid --
  km_levels <- if (identical(km_modules, "all")) MODULE_LEVELS else {
    sig_cols <- cells %>% filter(sig) %>% pull(score_col) %>% unique()
    keep <- MODULES[names(MODULES) %in% union(sig_cols, c("Up_ssGSEA", "Composite_ssGSEA"))]
    MODULE_LEVELS[MODULE_LEVELS %in% unname(keep)]
  }
  cv <- curves %>% filter(module %in% km_levels) %>% mutate(module = factor(module, km_levels))
  drawn <- distinct(cv, key, module)
  cc <- cells %>% filter(module %in% km_levels) %>% semi_join(drawn, by = c("key", "module")) %>%
    mutate(module = factor(module, km_levels),
           lab = paste0("log-rank FDR ", ifelse(is.na(lr_fdr), "NA",
                          ifelse(lr_fdr < 0.001, "p<0.001", sprintf("p=%.3f", lr_fdr))),
                        ifelse(lr_sig, " *", "")))

  build_km_block <- function(keys, title, show_y) {
    dc <- filter(cv, key %in% keys); da <- filter(cc, key %in% keys)
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
      labs(title = title, x = "Months since landmark", y = if (show_y) "Survival probability" else NULL) +
      .base_theme +
      theme(plot.title = element_text(face = "bold", size = 11, hjust = 0.5),
            strip.placement = "outside",
            strip.text.y.left = if (show_y) element_text(angle = 0, face = "bold", size = 8.5) else element_blank(),
            strip.text.x = element_text(size = 8, lineheight = 0.85),
            panel.spacing.x = unit(11, "pt"), panel.spacing.y = unit(5, "pt"),
            panel.grid.major.y = element_line(colour = "grey93", linewidth = 0.3))
  }
  block_m <- build_km_block(c("M_whole_OS","M_whole_RFS","M_trt_OS","M_trt_RFS"), "Marisa (GSE39582)", TRUE)
  block_t <- build_km_block(c("T_whole_OS","T_whole_RFS"), "TCGA-COAD", FALSE)
  fig_km <- block_m + block_t +
    plot_layout(widths = c(4, 2), guides = "collect") +
    plot_annotation(
      title = if (identical(km_modules, "all")) "Kaplan-Meier survival by signature score (median High/Low split)"
              else "Kaplan-Meier survival by score, outcome-associated modules",
      subtitle = "Median High/Low split; p = log-rank BH-FDR (* = significant). Continuous-Cox HR effect sizes are in the matrix panel. 3-yr landmark, x capped 36 mo.",
      theme = theme(plot.title = element_text(face = "bold", size = 14, hjust = 0),
                    plot.subtitle = element_text(size = 10, colour = "grey30", hjust = 0))) &
    theme(legend.position = "bottom")

  n_rows <- length(km_levels); km_h <- 2.0 * n_rows + 1.4
  .save_fig(fig_matrix, "Fig2A_km_hr_matrix", 9.5, 4.6, out_dir)
  .save_fig(fig_km,     "Fig2B_km_curves",    11.0, km_h, out_dir)
  composite <- wrap_elements(fig_matrix) / wrap_elements(fig_km) +
    plot_layout(heights = c(4.6, km_h - 0.4)) + plot_annotation(tag_levels = "A") &
    theme(plot.tag = element_text(face = "bold", size = 16))
  .save_fig(composite, "Fig2_km_composite", 11.5, 4.6 + km_h + 0.6, out_dir)
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
                 Metric  = c("OS","RFS","OS","RFS","OS","RFS"),
                 clab = c("Marisa\n(OS)","Marisa\n(RFS)","Marisa treated\n(OS)","Marisa treated\n(RFS)","TCGA\n(OS)","TCGA\n(RFS)"))
  col_lab_map <- setNames(COLS$clab, COLS$ckey)
  DS_PREFIX <- c("GSE39582" = "GSE", "GSE39582 (treated)" = "GSEt", "TCGA-COAD" = "TCGA")
  DIR_PAL <- c("Higher hazard (HR>1)" = "#E64B35", "Lower hazard (HR<1)" = "#3B6EA5")

  subg  <- .rd(run_dir, "Subgroup_Score_HRs.csv")
  inter <- .rd(run_dir, "Interaction_Cox_Summary.csv")
  clin  <- list(GSE39582 = .rd(run_dir, "GSE39582_clinical.csv"),
                "TCGA-COAD" = .rd(run_dir, "TCGA_COAD_clinical.csv"))
  # Treated-Marisa subset so sd_of() rescales the treated columns on the treated SD.
  clin[["GSE39582 (treated)"]] <- clin$GSE39582[which(clin$GSE39582$Chemo_adj == "Y"), , drop = FALSE]
  metric_evt <- c(OS = "OS3Y_event", RFS = "RFS3Y_event")
  sd_of <- function(ds, metric, score) {
    d <- clin[[ds]]; ev <- metric_evt[[metric]]
    if (!all(c(score, ev) %in% names(d))) return(NA_real_)
    v <- suppressWarnings(as.numeric(d[[score]])); ok <- !is.na(v) & !is.na(d[[ev]])
    stats::sd(v[ok])
  }

  subg <- subg %>%
    mutate(ckey = paste0(unname(DS_PREFIX[Dataset]), "_", Metric),
           sd_sc = mapply(sd_of, Dataset, Metric, Score),
           hr_sd = exp(log(HR) * sd_sc), lo_sd = exp(log(HR_lower) * sd_sc), hi_sd = exp(log(HR_upper) * sd_sc),
           Modifier = factor(Modifier, MOD_LEVELS), Level = factor(Level, LEVEL_ORDER),
           ckey = factor(ckey, COLS$ckey),
           offscale = hr_sd < FOREST_XLIM[1] | hr_sd > FOREST_XLIM[2],
           hr_d = squish(hr_sd, FOREST_XLIM), lo_d = squish(lo_sd, FOREST_XLIM), hi_d = squish(hi_sd, FOREST_XLIM))
  inter <- inter %>%
    mutate(ckey = paste0(unname(DS_PREFIX[Dataset]), "_", Metric),
           Modifier = factor(Modifier, MOD_LEVELS), ckey = factor(ckey, COLS$ckey),
           int_lab = paste0("interaction FDR ", ifelse(is.na(FDR_P), "NA",
                              ifelse(FDR_P < 0.001, "p<0.001", sprintf("p=%.3f", FDR_P))), .sig_stars(FDR_P)),
           int_sig = isTRUE(Is_Significant))

  # -- Fig 3A: interaction matrix --
  mat <- inter %>% filter(Score %in% names(SCORES)) %>%
    mutate(Score = factor(SCORES[Score], unname(SCORES)),
           Modifier = fct_rev(factor(MOD_LABEL[as.character(Modifier)], unname(MOD_LABEL))),
           mlog = pmin(-log10(FDR_P), 4))
  fig_matrix <- ggplot(mat, aes(x = ckey, y = Modifier, fill = mlog)) +
    geom_tile(colour = "white", linewidth = 0.6) +
    geom_text(aes(label = .sig_stars(FDR_P)), fontface = "bold", size = 4.2, vjust = 0.78) +
    facet_grid(. ~ Score) +
    scale_fill_gradient(low = "white", high = "#762A83", limits = c(0, 4),
      name = expression(-log[10]~"interaction FDR"),
      guide = guide_colourbar(barwidth = 9, barheight = 0.5, title.vjust = 1)) +
    scale_x_discrete(labels = col_lab_map) +
    labs(title = "Effect modification of the DTP signature (score-by-subgroup interaction)",
         subtitle = "Likelihood-ratio test of score*modifier vs score+modifier (landmark Cox).  * marks FDR-significant (BH, p<0.05).",
         x = NULL, y = NULL) +
    .base_theme + theme(axis.text.x = element_text(size = 8, lineheight = 0.85),
                        panel.spacing.x = unit(6, "pt"))

  # -- Fig 3B: per-score forests --
  build_forest <- function(score_col) {
    ds <- subg  %>% filter(Score == score_col)
    di <- inter %>% filter(Score == score_col)
    ds$dir <- ifelse(ds$hr_sd > 1, "Higher hazard (HR>1)", "Lower hazard (HR<1)")
    ggplot(ds, aes(x = hr_d, y = fct_rev(Level))) +
      geom_vline(xintercept = 1, linetype = "dashed", colour = "grey55") +
      geom_errorbarh(aes(xmin = lo_d, xmax = hi_d, colour = dir), height = 0.28, linewidth = 0.5) +
      geom_point(aes(colour = dir, shape = offscale), size = 2.1) +
      # Per-cohort n / events, drawn INSIDE each facet (cohort x endpoint) at the
      # left edge — the counts differ by cohort and by OS/RFS, so they cannot live
      # on the single shared y-axis label without mislabelling 3 of the 4 columns.
      geom_text(inherit.aes = FALSE,
                aes(x = FOREST_XLIM[1], y = fct_rev(Level),
                    label = sprintf("n=%d, e=%d", as.integer(N), as.integer(N_events))),
                hjust = 0, vjust = -0.5, size = 1.95, colour = "grey45") +
      geom_text(data = di, inherit.aes = FALSE, aes(x = FOREST_XLIM[1], y = Inf, label = int_lab),
                hjust = 0, vjust = 1.4, size = 2.5, colour = ifelse(di$int_sig, "#B2182B", "grey35")) +
      scale_y_discrete() +
      facet_grid(Modifier ~ ckey, scales = "free_y", space = "free_y", switch = "y",
                 labeller = labeller(ckey = col_lab_map, Modifier = function(x) MOD_LABEL[x])) +
      scale_x_log10(limits = FOREST_XLIM, breaks = c(0.1, 0.3, 1, 3, 10),
                    labels = c("0.1","0.3","1","3","10")) +
      scale_colour_manual(values = DIR_PAL, name = "Score effect", drop = FALSE) +
      scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 21),
                         labels = c(`FALSE` = "in range", `TRUE` = "off-scale (squished)"), name = NULL) +
      labs(title = paste0("Per-subgroup hazard ratio of ", SCORES[[score_col]], " score"),
           subtitle = "Continuous-score Cox HR per 1 SD within each subgroup (95% CI); x truncated to [0.1, 10]. n / e = subgroup size / events in that cohort.",
           x = "Hazard ratio per 1 SD (log scale)", y = NULL) +
      .base_theme +
      theme(strip.placement = "outside",
            strip.text.y.left = element_text(angle = 0, face = "bold", size = 8.5),
            strip.text.x = element_text(size = 8, lineheight = 0.85),
            axis.text.y = element_text(size = 7.5), panel.spacing = unit(4, "pt"),
            panel.grid.major.x = element_line(colour = "grey93", linewidth = 0.3))
  }
  .save_fig(fig_matrix, "Fig3A_interaction_matrix", 12.5, 3.4, out_dir)
  forests <- list()
  for (sc in names(SCORES)) {
    forests[[sc]] <- build_forest(sc)
    .save_fig(forests[[sc]], paste0("Fig3B_subgroup_forest_", sub("_ssGSEA", "", sc)), 14.0, 8.0, out_dir)
  }
  composite <- wrap_elements(fig_matrix) / wrap_elements(forests[[primary_score]]) +
    plot_layout(heights = c(3.4, 8.2)) + plot_annotation(tag_levels = "A") &
    theme(plot.tag = element_text(face = "bold", size = 16))
  .save_fig(composite, "Fig3_subgroup_composite", 14.0, 12.2, out_dir)
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
  metric_cols <- function(m) if (m == "OS") c(t = "OS3Y_delay", e = "OS3Y_event")
                             else            c(t = "RFS3Y_delay", e = "RFS3Y_event")
  # Within-subgroup score SD (non-NA score/time/event), computed once per member set.
  sd_of <- function(ds, project, score, metric) {
    d <- clin[[ds]]; mc <- metric_cols(metric)
    if (!all(c(score, mc[["t"]], mc[["e"]]) %in% names(d))) return(NA_real_)
    mem <- .subgroup_members(d, project, ds); d <- d[mem, , drop = FALSE]
    v <- suppressWarnings(as.numeric(d[[score]]))
    ok <- !is.na(v) & !is.na(d[[mc[["t"]]]]) & !is.na(d[[mc[["e"]]]])
    if (sum(ok) < 2) return(NA_real_)
    stats::sd(v[ok])
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
  proj$dtag     <- ifelse(proj$Dataset == "GSE39582", "GSE", "TCGA")
  proj$ylab     <- paste0(proj$dtag, " · ", proj$Subgroup)
  CAT_LEVELS <- intersect(c("CMS","PDS","Stage","MSI"), unique(proj$Category))

  END_LAB   <- c(OS = "3-yr OS", RFS = "3-yr RFS")
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
           tile_label = ifelse(!Is_Testable, "n/t",
                        ifelse(Cox_sig, .sig_stars(Cox_FDR_P, ""), "")))
  lim <- min(max(abs(pa$log2HR[is.finite(pa$log2HR)]), na.rm = TRUE), 2)
  fig_A <- ggplot(pa, aes(x = ckey, y = ylab, fill = pmax(pmin(log2HR, lim), -lim))) +
    geom_tile(colour = "white", linewidth = 0.6) +
    geom_text(aes(label = tile_label),
              colour = ifelse(pa$tile_label == "n/t", "grey45", "grey10"),
              fontface = "bold", size = 3.1, vjust = 0.78) +
    facet_grid(Category ~ ., scales = "free_y", space = "free_y", switch = "y") +
    scale_fill_gradient2(low = "#3B6EA5", mid = "white", high = "#E64B35", midpoint = 0,
      limits = c(-lim, lim), name = expression(log[2]~"HR (per 1 SD of score)"),
      guide = guide_colourbar(barwidth = 9, barheight = 0.5, title.vjust = 1), na.value = "grey85") +
    scale_x_discrete(labels = ck_lab) +
    labs(title = "DTP score vs 3-year outcome within molecular and clinical subgroups",
         subtitle = paste0("Univariable Cox HR per 1 SD (landmark 36 mo). * = FDR<0.05.\n",
                           "Grey \"n/t\" = not testable (< ", MIN_EVENTS, " events / ", MIN_COX_N,
                           " n), distinct from tested-but-null.\n",
                           "Most subgroups are underpowered: a non-significant or not-testable cell is NOT evidence of absence."),
         x = NULL, y = NULL) +
    .base_theme + theme(axis.text.x = element_text(size = 8, lineheight = 0.85),
                        axis.text.y = element_text(size = 8),
                        strip.placement = "outside",
                        strip.text.y.left = element_text(angle = 0, face = "bold", size = 9),
                        panel.spacing.y = unit(4, "pt"))

  # ============= Panel B: score-across-subtype (Kruskal-Wallis eps^2) =============
  pb <- sss %>% filter(Score %in% .CORE_SCORES, Subtype_Axis %in% c("CMS", "PDS")) %>%
    mutate(Score = factor(.CORE_LABEL[Score], rev(unname(.CORE_LABEL))),
           Axis  = factor(paste0(Subtype_Axis, " subtype"), c("CMS subtype", "PDS subtype")),
           eps2  = suppressWarnings(as.numeric(Eps2)),
           star  = .sig_stars(FDR_P, ""))
  fig_B <- ggplot(pb, aes(x = Dataset, y = Score, fill = eps2)) +
    geom_tile(colour = "white", linewidth = 0.6) +
    geom_text(aes(label = ifelse(star == "", sprintf("%.02f", eps2),
                                 paste0(sprintf("%.02f", eps2), " ", star))),
              size = 3.0, fontface = "bold", colour = "grey10") +
    facet_grid(. ~ Axis) +
    scale_fill_gradient(low = "white", high = "#762A83", limits = c(0, NA),
      name = expression(epsilon^2~"(effect size)"),
      guide = guide_colourbar(barwidth = 9, barheight = 0.5, title.vjust = 1)) +
    labs(title = "Does the DTP score itself differ across subtypes? (Kruskal-Wallis)",
         subtitle = "Epsilon-squared effect size across subtype levels; * = FDR<0.05 (BH).",
         x = NULL, y = NULL) +
    .base_theme + theme(axis.text.x = element_text(size = 9),
                        panel.spacing.x = unit(6, "pt"))

  # ---- Save panels + composite ----
  n_sub <- nlevels(cells$ylab); a_h <- max(4.5, 0.32 * n_sub + 2.6)
  .save_fig(fig_A, "Fig3_3A_subgroup_hr_matrix", 8.8, a_h, out_dir)
  .save_fig(fig_B, "Fig3_3B_score_across_subtype", 7.5, 3.2, out_dir)
  composite <- wrap_elements(fig_A) / wrap_elements(fig_B) +
    plot_layout(heights = c(a_h, 3.4)) + plot_annotation(tag_levels = "A") &
    theme(plot.tag = element_text(face = "bold", size = 16))
  .save_fig(composite, "Fig3_subtype_survival_composite", 9.5, a_h + 3.9, out_dir)

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
                 "4 (subtype survival)" = build_subtype_survival_figures)
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

# Auto-run only when executed as a script (skipped when source()'d by a module).
if (sys.nframe() == 0L) {
  .args <- commandArgs(trailingOnly = TRUE)
  build_crc_composites(if (length(.args) >= 1) .args[[1]] else .find_run_dir())
}
