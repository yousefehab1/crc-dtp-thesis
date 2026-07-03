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

# ------------------------------------------------------------------------------
# Confounder + effect-modification Cox (CRC only). All models use the 3-year
# landmarked endpoints and the CONTINUOUS score, consistent with get_cox_stats().
# ------------------------------------------------------------------------------

# Harmonise the four modifiers into clean, unified factor columns so one code
# path serves both cohorts. cohort is "GSE39582" or "TCGA-COAD".
harmonize_crc_modifiers <- function(df, cohort) {
  stage_raw <- if (cohort == "GSE39582") suppressWarnings(as.numeric(df$TNM_stage))
               else df$Stage
  df$Stage_bin <- if (cohort == "GSE39582")
    factor(dplyr::case_when(stage_raw %in% 0:2 ~ "Early",
                            stage_raw %in% 3:4 ~ "Late", TRUE ~ NA_character_),
           levels = c("Early", "Late"))
  else
    factor(dplyr::case_when(df$Stage %in% c("I", "II")   ~ "Early",
                            df$Stage %in% c("III", "IV")  ~ "Late", TRUE ~ NA_character_),
           levels = c("Early", "Late"))

  msi_raw <- if (cohort == "GSE39582") df$MMR else df[["paper_MSI_status"]]
  # MSI-L is grouped with MSS, not MSI: MSI-L lacks the hypermutator/immune phenotype that
  # defines the distinct MSI-high (dMMR) entity and behaves like MSS (Rantanen 2023 pooled
  # MSS/MSI-low vs MSI-high). This also makes GSE (dMMR/pMMR only) and TCGA consistent.
  df$MSI_group <- factor(dplyr::case_when(
    msi_raw %in% c("pMMR", "MSS", "MSI-L") ~ "MSS",
    msi_raw %in% c("dMMR", "MSI-H")        ~ "MSI",
    TRUE                                   ~ NA_character_), levels = c("MSS", "MSI"))

  if ("CMS" %in% colnames(df)) {
    cms <- df$CMS; cms[!cms %in% paste0("CMS", 1:4)] <- NA
    df$CMS <- stats::relevel(factor(cms), ref = "CMS2")   # canonical majority ref
  }
  if ("PDS" %in% colnames(df)) {
    pds <- df$PDS; pds[!pds %in% paste0("PDS", 1:3)] <- NA   # drop "Mixed" + NA
    df$PDS <- droplevels(factor(pds))
    if (nlevels(df$PDS) > 0) {
      ref <- names(which.max(table(df$PDS)))
      df$PDS <- stats::relevel(df$PDS, ref = ref)
    }
  }
  df
}

# Second-row p-value from an anova.coxph LRT table, robust to the p-column name.
.anova_lrt_p <- function(reduced, full) {
  av   <- stats::anova(reduced, full)
  pcol <- grep("Chi", colnames(av), value = TRUE)
  if (length(pcol) == 0) return(NA_real_)
  av[[tail(pcol, 1)]][2]
}

# Events-per-variable gate: TRUE when the model is admissible.
.cox_epv_ok <- function(n, n_events, n_params)
  n >= MIN_COX_N && n_events >= MIN_EVENTS && n_events >= MIN_EPV * n_params

# Drop rows missing any modelled column, then drop unused factor levels.
.cox_complete <- function(df, cols) {
  d <- df[stats::complete.cases(df[, cols, drop = FALSE]), , drop = FALSE]
  droplevels(d)
}

# Adjusted Cox: score term after adjusting for `covars` (character vector).
# Returns the SCORE-term HR/CI/p plus a model LRT (does score add signal beyond
# covars) and a PH check (cox.zph for the score term). NULL if inadmissible.
get_cox_adjusted <- function(df, score_col, t_col, e_col, covars) {
  cols <- c(score_col, t_col, e_col, covars)
  d <- .cox_complete(df, cols)
  covars <- covars[vapply(covars, function(c) length(unique(d[[c]])) > 1, logical(1))]
  n_params <- length(covars) + 1L +
    sum(vapply(covars, function(c) if (is.factor(d[[c]])) nlevels(d[[c]]) - 2L else 0L, integer(1)))
  n_ev <- sum(d[[e_col]])
  if (!.cox_epv_ok(nrow(d), n_ev, n_params)) return(NULL)
  rhs  <- paste(c(score_col, covars), collapse = " + ")
  base <- paste(covars, collapse = " + ")
  tryCatch({
    full <- survival::coxph(as.formula(sprintf("Surv(%s,%s) ~ %s", t_col, e_col, rhs)),
                            data = d, ties = "efron")
    s    <- summary(full)
    lrt  <- if (length(covars) == 0) unname(s$logtest[["pvalue"]]) else {
      red <- survival::coxph(as.formula(sprintf("Surv(%s,%s) ~ %s", t_col, e_col, base)), data = d)
      .anova_lrt_p(red, full)
    }
    ph <- tryCatch(survival::cox.zph(full)$table[score_col, "p"], error = function(e) NA_real_)
    adj_hr <- exp(coef(full)[[score_col]])
    # unadjusted score fit on the SAME complete-case subset -> clean Δ estimate
    unadj  <- survival::coxph(as.formula(sprintf("Surv(%s,%s) ~ %s", t_col, e_col, score_col)), data = d)
    unadj_hr <- exp(coef(unadj)[[1]])
    # Change-in-estimate for confounding. Delta_logHR is the stable absolute change in the
    # score's log-HR after adjustment; the percent version is floored so it does not explode
    # when the unadjusted log-HR is ~0 (HR near 1) — read the absolute value in that regime.
    delta_loghr <- log(adj_hr) - log(unadj_hr)
    delta_pct   <- 100 * delta_loghr / max(abs(log(unadj_hr)), 0.05)
    list(HR = adj_hr,
         HR_lower = s$conf.int[score_col, "lower .95"],
         HR_upper = s$conf.int[score_col, "upper .95"],
         P = s$coefficients[score_col, "Pr(>|z|)"], C_index = unname(s$concordance["C"]),
         LRT_P = lrt, PH_P = ph, Unadj_HR = unadj_hr,
         Delta_logHR = delta_loghr, Delta_logHR_pct = delta_pct,
         N = nrow(d), N_events = n_ev)
  }, error = function(e) NULL)
}

# Interaction Cox: LRT of `score * modifier` vs `score + modifier` (effect
# modification), plus the per-level simple score HRs. NULL if inadmissible.
get_cox_interaction <- function(df, score_col, t_col, e_col, modifier) {
  cols <- c(score_col, t_col, e_col, modifier)
  d <- .cox_complete(df, cols)
  if (!is.factor(d[[modifier]])) d[[modifier]] <- factor(d[[modifier]])
  d <- droplevels(d)
  k <- nlevels(d[[modifier]])
  if (k < 2) return(NULL)
  n_params <- 2L * k - 1L                       # score + (k-1) main + (k-1) interaction
  n_ev <- sum(d[[e_col]])
  if (!.cox_epv_ok(nrow(d), n_ev, n_params)) return(NULL)
  tryCatch({
    red  <- survival::coxph(as.formula(sprintf("Surv(%s,%s) ~ %s + %s", t_col, e_col, score_col, modifier)), data = d)
    full <- survival::coxph(as.formula(sprintf("Surv(%s,%s) ~ %s * %s", t_col, e_col, score_col, modifier)), data = d)
    lrt  <- .anova_lrt_p(red, full)
    # per-level simple score HRs (univariable score within each modifier level)
    per_level <- lapply(levels(d[[modifier]]), function(lv) {
      cs <- get_cox_stats(d[d[[modifier]] == lv, , drop = FALSE], score_col, t_col, e_col)
      if (is.null(cs)) return(NULL)
      data.frame(Level = lv, HR = cs$HR, HR_lower = cs$HR_lower, HR_upper = cs$HR_upper,
                 N = cs$N, N_events = cs$N_events, stringsAsFactors = FALSE)
    })
    list(Interaction_P = lrt, K = k, N = nrow(d), N_events = n_ev,
         per_level = do.call(rbind, per_level))
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
