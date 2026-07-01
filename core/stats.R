# ==============================================================================
# core/stats.R  —  Shared statistics: continuous-score Cox (#7), effect size +
# testability (#8), unified gating thresholds (#9), provisional FDR (#10),
# and a high-level survival runner so both modules treat stats identically.
# ==============================================================================

metric_cols <- function(metric)
  {if (metric == "OS") c(t = "OS3Y_delay", e = "OS3Y_event")
  else                c(t = "RFS3Y_delay", e = "RFS3Y_event")
}
# --- Wilcoxon + rank-biserial effect size (#8). Gate: MIN_GROUP_N (#9). -------
# r sign is relative to the FIRST factor level (Dead_3yr / Recurred).
get_wilcox_stats <- function(df, score_col, group_col) {
  d <- df %>% dplyr::filter(!is.na(.data[[group_col]]) & .data[[group_col]] != "")
  if (!is.numeric(d[[score_col]])) return(list(p = NA, r = NA, n = NA))
  tbl <- table(droplevels(factor(d[[group_col]])))
  if (length(tbl) < 2 || any(tbl < MIN_GROUP_N)) return(list(p = NA, r = NA, n = NA))
  lvls <- names(tbl)
  g1 <- d[[score_col]][d[[group_col]] == lvls[1]]
  g2 <- d[[score_col]][d[[group_col]] == lvls[2]]
  wt <- suppressWarnings(wilcox.test(g1, g2, exact = FALSE))
  r  <- as.numeric(1 - 2 * wt$statistic / (length(g1) * length(g2)))
  list(p = wt$p.value, r = r, n = length(g1) + length(g2))
}

# --- Kruskal-Wallis + epsilon-squared effect for score-across-subtype (>=2 --
# groups). Gate: MIN_GROUP_N per retained level. Used to compare a continuous
# score across molecular subtypes (CMS1-4, PDS1-3), unlike the 2-group Wilcoxon.
get_kruskal_stats <- function(df, score_col, group_col) {
  d <- df %>% dplyr::filter(!is.na(.data[[group_col]]) & .data[[group_col]] != "" &
                            !is.na(.data[[score_col]]))
  if (!is.numeric(d[[score_col]])) return(list(p = NA, eps2 = NA, n = NA, k = NA))
  tbl  <- table(droplevels(factor(d[[group_col]])))
  keep <- names(tbl)[tbl >= MIN_GROUP_N]
  d    <- d[d[[group_col]] %in% keep, , drop = FALSE]
  if (length(keep) < 2) return(list(p = NA, eps2 = NA, n = nrow(d), k = length(keep)))
  kt <- tryCatch(kruskal.test(d[[score_col]], factor(d[[group_col]])),
                 error = function(e) NULL)
  if (is.null(kt)) return(list(p = NA, eps2 = NA, n = nrow(d), k = length(keep)))
  n    <- nrow(d)
  eps2 <- as.numeric(kt$statistic) / ((n^2 - 1) / (n + 1))   # epsilon-squared
  list(p = kt$p.value, eps2 = eps2, n = n, k = length(keep))
}

# --- Log-rank p on median-split KM. Gate: MIN_KM_GROUP_N + MIN_EVENTS (#9). ---
get_km_pval <- function(df, score_col, t_col, e_col) {
  d <- df[!is.na(df[[t_col]]) & !is.na(df[[e_col]]) & !is.na(df[[score_col]]), ]
  if (nrow(d) == 0) return(NA)
  d$Group <- factor(ifelse(d[[score_col]] > median(d[[score_col]], na.rm = TRUE), "High", "Low"),
                    levels = c("Low", "High"))
  if (length(table(d$Group)) < 2 || any(table(d$Group) < MIN_KM_GROUP_N) ||
      sum(d[[e_col]]) < MIN_EVENTS) return(NA)
  tryCatch({
    sd <- survival::survdiff(as.formula(paste0("Surv(", t_col, ",", e_col, ") ~ Group")), data = d)
    1 - pchisq(sd$chisq, length(sd$n) - 1)
  }, error = function(e) NA)
}

# --- Univariable Cox on the CONTINUOUS score (#7). Gate: MIN_COX_N+MIN_EVENTS -
get_cox_stats <- function(df, score_col, t_col, e_col) {
  d <- df[!is.na(df[[t_col]]) & !is.na(df[[e_col]]) & !is.na(df[[score_col]]), ]
  if (nrow(d) < MIN_COX_N || sum(d[[e_col]]) < MIN_EVENTS) return(NULL)
  tryCatch({
    fit <- survival::coxph(as.formula(paste0("Surv(", t_col, ",", e_col, ") ~ ", score_col)), data = d)
    s   <- summary(fit)
    list(HR = exp(coef(fit)[[1]]),
         HR_lower = s$conf.int[1, "lower .95"], HR_upper = s$conf.int[1, "upper .95"],
         P = s$coefficients[1, "Pr(>|z|)"], C_index = unname(s$concordance["C"]),
         N = nrow(d), N_events = sum(d[[e_col]]))
  }, error = function(e) NULL)
}

# --- One tidy stats row -------------------------------------------------------
.stat_row <- function(dataset, proj, test, metric, score,
                      raw_p = NA, effect_r = NA, hr = NA, hr_l = NA, hr_u = NA,
                      c_index = NA, n = NA, n_events = NA, is_testable = !is.na(raw_p)) {
  data.frame(Dataset = dataset, Project = proj, Test = test, Metric = metric, Score = score,
             Raw_P = raw_p, Effect_r = effect_r, HR = hr, HR_lower = hr_l, HR_upper = hr_u,
             C_index = c_index, N = n, N_events = n_events, Is_Testable = is_testable,
             stringsAsFactors = FALSE)
}

# Run all three tests for one score x metric; returns up to 3 rows.
run_all_stats <- function(df, score_col, metric, dataset, proj, score_name, group_col) {
  mc <- metric_cols(metric)
  ws <- get_wilcox_stats(df, score_col, group_col)
  km <- get_km_pval(df, score_col, mc[["t"]], mc[["e"]])
  cs <- get_cox_stats(df, score_col, mc[["t"]], mc[["e"]])
  rows <- list(
    .stat_row(dataset, proj, "Wilcoxon", metric, score_name,
              raw_p = ws$p, effect_r = ws$r, n = ws$n, is_testable = !is.na(ws$p)),
    .stat_row(dataset, proj, "KM", metric, score_name, raw_p = km, is_testable = !is.na(km)),
    if (is.null(cs))
      .stat_row(dataset, proj, "Cox", metric, score_name, is_testable = FALSE)
    else
      .stat_row(dataset, proj, "Cox", metric, score_name,
                raw_p = cs$P, hr = cs$HR, hr_l = cs$HR_lower, hr_u = cs$HR_upper,
                c_index = cs$C_index, n = cs$N, n_events = cs$N_events)
  )
  do.call(rbind, rows)
}

# ------------------------------------------------------------------------------
# FDR correction (#10) — PROVISIONAL.
# ------------------------------------------------------------------------------
# !!! The cross-module unified FDR strategy is an OPEN DECISION. Until it is
# !!! settled, each module calls apply_fdr() with the `by` key that REPLICATES
# !!! its original behaviour, so no result silently changes:
# !!!   CRC arm:        by = c("Dataset","Project","Test","Metric")
# !!!   pan-cancer arm: by = c("Family","Test")
# !!! Changing the agreed scheme later is a one-line edit at each call site.
# ------------------------------------------------------------------------------
apply_fdr <- function(stats_df, by) {
  stats_df %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(by))) %>%
    dplyr::mutate(FDR_P = p.adjust(Raw_P, method = "BH")) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(Is_Significant = !is.na(FDR_P) & FDR_P < 0.05)
}

# ------------------------------------------------------------------------------
# High-level runner: stats for a named list of cohorts x score columns.
# Returns the raw (un-FDR'd) stats rows; the caller applies apply_fdr().
# ------------------------------------------------------------------------------
run_survival_block <- function(cohorts, dataset, score_cols,
                               group_os = "Surv_3yr", group_rfs = "Recurrence_3yr") {
  rows <- list()
  for (cn in names(cohorts)) {
    cd <- cohorts[[cn]]
    for (sc in score_cols) {
      rows[[length(rows) + 1]] <- run_all_stats(cd, sc, "OS",  dataset, cn, sc, group_os)
      rows[[length(rows) + 1]] <- run_all_stats(cd, sc, "RFS", dataset, cn, sc, group_rfs)
    }
  }
  do.call(rbind, rows)
}
