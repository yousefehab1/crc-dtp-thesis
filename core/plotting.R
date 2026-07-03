# ==============================================================================
# core/plotting.R  —  ONE theme + ONE palette (#13), shared violin/KM/forest
# helpers, and three-way significance routing (#8).
# ==============================================================================

pub_theme <- ggplot2::theme_classic(base_size = 12) +
  ggplot2::theme(
    plot.title    = ggplot2::element_text(hjust = 0.5, face = "bold", size = 14),
    plot.subtitle = ggplot2::element_text(hjust = 0.5, face = "italic", size = 12),
    axis.title    = ggplot2::element_text(face = "bold"),
    axis.text     = ggplot2::element_text(colour = "black"),
    legend.position = "bottom",
    legend.title  = ggplot2::element_text(face = "bold"),
    strip.text    = ggplot2::element_text(face = "bold")
  )

# One palette covering every named state used across the three modules.
pub_palette <- c(
  "Dead_3yr"        = "#E64B35", "Alive_3yr"       = "#4DBBD5",
  "Recurred"        = "#E64B35", "Recurrence-Free" = "#4DBBD5",
  "High"            = "#E64B35", "Low"             = "#4DBBD5",
  "Normal"          = "#55A868", "Primary"         = "#4C72B0", "Metastasis" = "#C44E52"
)

# --- Significance lookup + folder routing (#8) --------------------------------
get_stat_info <- function(stats_df, dataset, proj, test, metric, score) {
  row <- stats_df %>%
    dplyr::filter(Dataset == dataset, Project == proj, Test == test,
                  Metric == metric, Score == score)
  if (nrow(row) == 0)
    return(list(sig = FALSE, fdr_p = NA_real_, effect_r = NA_real_, testable = FALSE))
  list(sig = isTRUE(any(row$Is_Significant)), fdr_p = row$FDR_P[1],
       effect_r = if ("Effect_r" %in% names(row)) row$Effect_r[1] else NA_real_,
       testable = isTRUE(any(row$Is_Testable)))
}

route_folder <- function(out_root, module, type, si) {
  sub <- if (!si$testable) "not_tested" else if (si$sig) "significant" else "non_significant"
  p <- file.path(out_root, module, paste0(type, "_plots"), sub)
  dir.create(p, recursive = TRUE, showWarnings = FALSE)
  p
}

# --- Violin (Wilcoxon) --------------------------------------------------------
generate_violin <- function(df, score_col, group_col, title, filename, si, out_root, module) {
  d <- df %>% dplyr::filter(!is.na(.data[[group_col]]) & .data[[group_col]] != "")
  if (!is.numeric(d[[score_col]])) return(invisible(NULL))
  tb <- table(droplevels(factor(d[[group_col]])))
  if (length(tb) < 2) return(invisible(NULL))

  out_path <- file.path(route_folder(out_root, module, "violin", si), paste0(filename, ".png"))

  # All patients are shown: geom_violin(scale = "width") makes the shapes
  # comparable regardless of group size, and the per-group n is annotated, so no
  # subsampling is needed for an honest comparison (and it would only add noise
  # to the density estimate).
  sm   <- d %>% dplyr::group_by(.data[[group_col]]) %>% dplyr::summarise(n = dplyr::n(), .groups = "drop")
  lbls <- setNames(paste0(sm[[group_col]], "\n(n=", sm$n, ")"), sm[[group_col]])
  ann  <- paste(Filter(nchar, c(
    if (!is.na(si$fdr_p))    sprintf("FDR p = %.3g", si$fdr_p) else "FDR p = N/A",
    if (!is.na(si$effect_r)) sprintf("r = %.3f", si$effect_r)  else "")),
    collapse = "  |  ")

  sig_subtitle <- if (!si$testable) "Not testable (insufficient data)"
                  else if (si$sig) "FDR Significant (p < 0.05)" else "Not FDR Significant"

  p <- ggplot2::ggplot(d, ggplot2::aes(x = .data[[group_col]], y = .data[[score_col]],
                                       fill = .data[[group_col]])) +
    ggplot2::geom_violin(alpha = 0.8, trim = FALSE, scale = "width") +
    ggplot2::geom_boxplot(width = 0.15, fill = "white", alpha = 0.9, outlier.shape = NA) +
    ggplot2::scale_fill_manual(values = pub_palette) +
    ggplot2::scale_x_discrete(labels = lbls) +
    ggplot2::annotate("text", x = 1.5, y = Inf, vjust = 2, label = ann,
                      size = 4, fontface = "italic", colour = "grey30") +
    ggplot2::labs(title = title,
                  subtitle = sig_subtitle,
                  x = NULL, y = paste0("ssGSEA: ", sub(SCORE_SUFFIX, "", score_col))) +
    pub_theme + ggplot2::theme(legend.position = "none")
  ggplot2::ggsave(out_path, p, width = 7, height = 6, dpi = 300, bg = "white")
}

# --- Score-across-subtype violin (Kruskal-Wallis; 2-4 groups) -----------------
# Distribution of one continuous score across molecular subtypes (CMS / PDS).
# `si` carries: testable, sig, fdr_p, eps2 (from get_kruskal_stats + apply_fdr).
generate_subtype_violin <- function(df, score_col, subtype_col, title, filename,
                                    si, out_root, module) {
  d <- df %>% dplyr::filter(!is.na(.data[[subtype_col]]) & .data[[subtype_col]] != "" &
                            !is.na(.data[[score_col]]))
  if (!is.numeric(d[[score_col]])) return(invisible(NULL))
  tb <- table(droplevels(factor(d[[subtype_col]])))
  if (length(tb) < 2) return(invisible(NULL))

  out_path <- file.path(route_folder(out_root, module, "subtype_violin", si),
                        paste0(filename, ".png"))

  # All patients are shown: geom_violin(scale = "width") keeps the subtype shapes
  # comparable despite unequal group sizes, and the per-group n is annotated.
  sm   <- d %>% dplyr::group_by(.data[[subtype_col]]) %>%
    dplyr::summarise(n = dplyr::n(), .groups = "drop")
  lbls <- setNames(paste0(sm[[subtype_col]], "\n(n=", sm$n, ")"), sm[[subtype_col]])
  ann  <- paste(Filter(nchar, c(
    if (!is.na(si$fdr_p)) sprintf("KW FDR p = %.3g", si$fdr_p) else "KW FDR p = N/A",
    if (!is.na(si$eps2)) sprintf("eps2 = %.3f", si$eps2) else "")),
    collapse = "  |  ")

  sig_subtitle <- if (!si$testable) "Not testable (insufficient data)"
                  else if (si$sig) "FDR Significant (p < 0.05)" else "Not FDR Significant"

  p <- ggplot2::ggplot(d, ggplot2::aes(x = .data[[subtype_col]], y = .data[[score_col]],
                                       fill = .data[[subtype_col]])) +
    ggplot2::geom_violin(alpha = 0.8, trim = FALSE, scale = "width") +
    ggplot2::geom_boxplot(width = 0.15, fill = "white", alpha = 0.9, outlier.shape = NA) +
    ggplot2::scale_fill_brewer(palette = "Set2") +
    ggplot2::scale_x_discrete(labels = lbls) +
    ggplot2::annotate("text", x = 1.5, y = Inf, vjust = 2, label = ann,
                      size = 4, fontface = "italic", colour = "grey30") +
    ggplot2::labs(title = title,
                  subtitle = sig_subtitle,
                  x = subtype_col, y = paste0("ssGSEA: ", sub(SCORE_SUFFIX, "", score_col))) +
    pub_theme + ggplot2::theme(legend.position = "none")
  ggplot2::ggsave(out_path, p, width = 7, height = 6, dpi = 300, bg = "white")
}

# --- Kaplan-Meier (median split for visualisation only; #7) -------------------
generate_km <- function(df, score_col, metric, title, filename, si, out_root, module) {
  mc <- metric_cols(metric); t_col <- mc[["t"]]; e_col <- mc[["e"]]
  d <- df[!is.na(df[[t_col]]) & !is.na(df[[e_col]]) & !is.na(df[[score_col]]), ]
  if (nrow(d) == 0) return(invisible(NULL))
  d$Group <- factor(ifelse(d[[score_col]] > median(d[[score_col]], na.rm = TRUE), "High", "Low"),
                    levels = c("Low", "High"))
  if (length(table(d$Group)) < 2 || any(table(d$Group) < MIN_KM_GROUP_N) ||
      sum(d[[e_col]]) < MIN_EVENTS) return(invisible(NULL))

  out_path <- file.path(route_folder(out_root, module, "km", si), paste0(filename, ".png"))
  ft <- paste0(title, "\n", if (!is.na(si$fdr_p)) sprintf("FDR p = %.3g", si$fdr_p) else "FDR p = N/A",
               if (si$sig) " (FDR Significant)" else "")
  tryCatch({
    fit <- eval(substitute(survival::survfit(Surv(T, E) ~ Group, data = d),
                           list(T = as.name(t_col), E = as.name(e_col))))
    pl  <- survminer::ggsurvplot(fit, data = d, pval = FALSE, risk.table = TRUE,
                                 risk.table.col = "strata", ggtheme = pub_theme, title = ft,
                                 xlab = sprintf("Time (months, truncated at %d)", OS_CUTPOINT),
                                 legend.labs = c("Low Score", "High Score"),
                                 palette = unname(c(pub_palette["Low"], pub_palette["High"])),
                                 censor = nrow(d) < 500)
    grDevices::png(out_path, width = 2000, height = 1800, res = 300)
    tryCatch(print(pl), error = function(e) message("  [KM print] ", e$message),
             finally = if (grDevices::dev.cur() > 1) grDevices::dev.off())
  }, error = function(e) message("  [KM] ", e$message))
}

# --- Forest plot: continuous-Cox HR across cohorts (#7) -----------------------
generate_forest <- function(stats_df, dataset, metric, score_name, out_root, module,
                            exclude_proj = NULL) {
  cd <- stats_df %>%
    dplyr::filter(Dataset == dataset, Test == "Cox", Metric == metric,
                  Score == score_name, Is_Testable, !is.na(HR))
  if (!is.null(exclude_proj)) cd <- cd %>% dplyr::filter(!Project %in% exclude_proj)
  if (nrow(cd) == 0) return(invisible(NULL))
  cd <- cd %>% dplyr::arrange(HR) %>%
    dplyr::mutate(Project = factor(Project, levels = Project),
                  Sig = ifelse(Is_Significant, "FDR Sig.", "NS"))

  p <- ggplot2::ggplot(cd, ggplot2::aes(x = HR, y = Project, colour = Sig)) +
    ggplot2::geom_vline(xintercept = 1, linetype = "dashed", colour = "grey50") +
    ggplot2::geom_point(size = 2.5) +
    ggplot2::geom_errorbarh(ggplot2::aes(xmin = HR_lower, xmax = HR_upper), height = 0.3) +
    ggplot2::scale_x_log10(name = "Hazard Ratio (log scale, 95% CI)") +
    ggplot2::scale_colour_manual(values = c("FDR Sig." = "#E64B35", "NS" = "grey55"), name = NULL) +
    ggplot2::labs(title = sprintf("Forest: %s  %s vs %s", dataset, sub(SCORE_SUFFIX, "", score_name), metric),
                  subtitle = sprintf("Univariable continuous-score Cox PH | n = %d cohort(s)", nrow(cd)),
                  y = NULL) +
    pub_theme + ggplot2::theme(legend.position = "top")

  dd <- file.path(out_root, module, "forest_plots"); dir.create(dd, recursive = TRUE, showWarnings = FALSE)
  # Dataset MUST be in the filename: crc_survival calls this per dataset, so
  # without it TCGA-COAD overwrites the GSE39582 forest of the same metric/score.
  ggplot2::ggsave(file.path(dd, sprintf("Forest_%s_%s_%s.png", dataset, metric, sub(SCORE_SUFFIX, "", score_name))),
                  p, width = 7, height = max(3.5, nrow(cd) * 0.32 + 1.5), dpi = 300, bg = "white")
}

# --- Subgroup forest: per-level score HR within one modifier ------------------
# hr_df: data.frame(Level, HR, HR_lower, HR_upper, N, N_events). Annotated with
# the effect-modification (interaction) FDR p for the modifier.
generate_subgroup_forest <- function(hr_df, dataset, metric, score_name, modifier,
                                     interaction_fdr, out_root, module) {
  cd <- hr_df[!is.na(hr_df$HR), , drop = FALSE]
  if (nrow(cd) == 0) return(invisible(NULL))
  cd$Label <- sprintf("%s (n=%d, e=%d)", cd$Level, cd$N, cd$N_events)
  cd <- cd[order(cd$HR), ]
  cd$Label <- factor(cd$Label, levels = cd$Label)

  ann <- if (!is.na(interaction_fdr))
    sprintf("Interaction FDR p = %.3g%s", interaction_fdr,
            if (interaction_fdr < 0.05) " (effect modification)" else "") else "Interaction p = N/A"

  p <- ggplot2::ggplot(cd, ggplot2::aes(x = HR, y = Label)) +
    ggplot2::geom_vline(xintercept = 1, linetype = "dashed", colour = "grey50") +
    ggplot2::geom_point(size = 2.5, colour = "#4C72B0") +
    ggplot2::geom_errorbarh(ggplot2::aes(xmin = HR_lower, xmax = HR_upper), height = 0.3) +
    ggplot2::scale_x_log10(name = "Score HR within subgroup (log scale, 95% CI)") +
    ggplot2::labs(title = sprintf("%s  %s vs %s  by %s",
                                  dataset, sub(SCORE_SUFFIX, "", score_name), metric, modifier),
                  subtitle = ann, y = NULL) +
    pub_theme + ggplot2::theme(legend.position = "none")

  dd <- file.path(out_root, module, "subgroup_forest_plots")
  dir.create(dd, recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(file.path(dd, sprintf("Subgroup_%s_%s_%s_by_%s.png",
                            dataset, metric, sub(SCORE_SUFFIX, "", score_name), modifier)),
                  p, width = 7, height = max(3, nrow(cd) * 0.4 + 1.5), dpi = 300, bg = "white")
}

# Convenience: run all violin + KM plots for a cohort block.
plot_survival_block <- function(cohorts, dataset, score_cols, stats_df, out_root, module,
                                group_os = "Surv_3yr", group_rfs = "Recurrence_3yr") {
  for (cn in names(cohorts)) {
    cd <- cohorts[[cn]]
    for (sc in score_cols) {
      generate_violin(cd, sc, group_os,
                      paste(cn, sub(SCORE_SUFFIX, "", sc), "vs OS"),
                      paste0(cn, "_OS_", sc),
                      get_stat_info(stats_df, dataset, cn, "Wilcoxon", "OS", sc), out_root, module)
      generate_violin(cd, sc, group_rfs,
                      paste(cn, sub(SCORE_SUFFIX, "", sc), "vs RFS"),
                      paste0(cn, "_RFS_", sc),
                      get_stat_info(stats_df, dataset, cn, "Wilcoxon", "RFS", sc), out_root, module)
      generate_km(cd, sc, "OS",
                  paste(cn, sub(SCORE_SUFFIX, "", sc), "OS"),
                  paste0(cn, "_OS_", sc),
                  get_stat_info(stats_df, dataset, cn, "KM", "OS", sc), out_root, module)
      generate_km(cd, sc, "RFS",
                  paste(cn, sub(SCORE_SUFFIX, "", sc), "RFS"),
                  paste0(cn, "_RFS_", sc),
                  get_stat_info(stats_df, dataset, cn, "KM", "RFS", sc), out_root, module)
    }
  }
}
