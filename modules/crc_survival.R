# ==============================================================================
# modules/crc_survival.R  —  CRC ssGSEA & survival (GSE39582 + TCGA-COAD).
# Was Final.R. Subtyping (CMS/CRIS/PDS) has been pulled out to
# subtyping/crc_subtyping.R (#11) and is NOT run here.
#
# Fixes applied: TPM for TCGA (#2a/#2b), log-space collapse (#3), shared CDR
# endpoint=PFI (#4), Dead/Alive labels (#5), censoring-aware endpoints (#6),
# continuous-score Cox (#7), effect size + 3-way routing (#8), gating (#9),
# primary-only (#17). FDR is PROVISIONAL (#10): by=Dataset×Project×Test×Metric.
# ==============================================================================

run_crc_survival <- function(out_root, panel, crc_signatures = NULL) {
  module <- "crc_survival"
  message("\n==================== MODULE: CRC SURVIVAL ====================")
  dir.create(file.path(out_root, module), recursive = TRUE, showWarnings = FALSE)
  gs         <- build_gene_sets(panel, which = crc_signatures)   # NULL = whole panel
  all_stats  <- list()

  # ------------------------------------------------------------------------
  # PART A — GSE39582 (Affymetrix microarray)
  # ------------------------------------------------------------------------
  message("\n--- GSE39582 ---")
  eset    <- cache_rds("GSE39582_eset", function() { g <- GEOquery::getGEO("GSE39582", GSEMatrix = TRUE); g[[1]] })
  mat_sym <- prep_microarray_symbols(Biobase::exprs(eset))

  scores_gse <- run_ssgsea(mat_sym, gs)
  scores_gse <- add_composite(scores_gse)
  scores_gse$Sample_ID <- rownames(scores_gse)
  score_cols <- setdiff(colnames(scores_gse), "Sample_ID")

  clinical_df <- as.data.frame(Biobase::pData(eset))

  # GEO characteristics_ch1.N are positional — verify field order is consistent
  # across all samples before trusting the hardcoded index map (kept from Final.R).
  characteristics_cols_used <- c(
    "characteristics_ch1.2",  "characteristics_ch1.3",  "characteristics_ch1.4",
    "characteristics_ch1.9",  "characteristics_ch1.11", "characteristics_ch1.12",
    "characteristics_ch1.13", "characteristics_ch1.14", "characteristics_ch1.15",
    "characteristics_ch1.18", "characteristics_ch1.22", "characteristics_ch1.26")
  consistency_issues <- character(0)
  for (col in characteristics_cols_used) {
    vals   <- clinical_df[[col]]; vals <- vals[!is.na(vals)]
    labels <- unique(trimws(sub(":.*$", "", vals)))
    if (length(labels) > 1)
      consistency_issues <- c(consistency_issues,
                              paste0(col, " has ", length(labels), " labels: ", paste(labels, collapse = ", ")))
  }
  if (length(consistency_issues) > 0)
    stop("GSE39582 characteristics_ch1.N field order inconsistent:\n",
         paste(consistency_issues, collapse = "\n"))

  parse_numeric_field <- function(raw, field_name) {
    stripped <- stringr::str_remove(raw, "^.*: ")
    out      <- suppressWarnings(as.numeric(stripped))
    newly_na <- is.na(out) & !is.na(stripped)
    if (any(newly_na))
      message(sprintf("  %s: %d value(s) unparseable; e.g. %s", field_name, sum(newly_na),
                      paste(unique(stripped[newly_na])[seq_len(min(5, sum(newly_na)))], collapse = " | ")))
    out
  }

  clinical <- clinical_df %>%
    tibble::rownames_to_column(var = "Sample_ID") %>%
    dplyr::distinct(Sample_ID, .keep_all = TRUE) %>%
    dplyr::select(Sample_ID,
                  Age_raw = characteristics_ch1.3,  Sex_raw = characteristics_ch1.2,
                  Chemo_adj_raw = characteristics_ch1.9, TNM_stage_raw = characteristics_ch1.4,
                  MMR_raw = characteristics_ch1.15,
                  RFS_event_raw = characteristics_ch1.11, RFS_delay_raw = characteristics_ch1.12,
                  OS_event_raw = characteristics_ch1.13,  OS_delay_raw = characteristics_ch1.14,
                  TP53_raw = characteristics_ch1.18, KRAS_raw = characteristics_ch1.22, BRAF_raw = characteristics_ch1.26) %>%
    dplyr::mutate(
      Age       = parse_numeric_field(Age_raw, "Age"),
      RFS_event = parse_numeric_field(RFS_event_raw, "RFS_event"),
      RFS_delay = parse_numeric_field(RFS_delay_raw, "RFS_delay"),
      OS_event  = parse_numeric_field(OS_event_raw, "OS_event"),
      OS_delay  = parse_numeric_field(OS_delay_raw, "OS_delay"),
      Sex       = stringr::str_trim(stringr::str_remove(Sex_raw, "^.*: ")),
      Chemo_adj = stringr::str_trim(stringr::str_remove(Chemo_adj_raw, "^.*: ")),
      TNM_stage = stringr::str_trim(stringr::str_remove(TNM_stage_raw, "^.*: ")),
      MMR       = stringr::str_trim(stringr::str_remove(MMR_raw, "^.*: "))
    ) %>%
    dplyr::select(-dplyr::ends_with("_raw")) %>%
    dplyr::left_join(scores_gse, by = "Sample_ID") %>%
    # feed shared censoring-aware derivation (#5/#6)
    dplyr::mutate(Recurrence_event = RFS_event, Recurrence_delay = RFS_delay) %>%
    derive_endpoints(os_event = "OS_event", os_time = "OS_delay",
                     rfs_event = "Recurrence_event", rfs_time = "Recurrence_delay")

  # --- Molecular subtyping (CMS + PDS), symbol-keyed, on the same matrix ------
  emat_gse_sym <- build_symbol_matrix(Biobase::exprs(eset), "microarray")
  cms_gse <- cache_rds("cms_gse", function() call_cms(emat_gse_sym))
  pds_gse <- cache_rds("pds_gse", function() call_pds(emat_gse_sym))
  clinical <- clinical %>%
    attach_subtypes(cms_gse, clinical_key = "Sample_ID") %>%
    attach_subtypes(pds_gse, clinical_key = "Sample_ID")
  write.csv(csv_safe(dplyr::select(clinical, Sample_ID, dplyr::any_of(c("CMS", "PDS")))),
            file.path(out_root, module, "GSE39582_subtypes.csv"), row.names = FALSE)

  write.csv(csv_safe(clinical), file.path(out_root, module, "GSE39582_clinical.csv"), row.names = FALSE)

  gse_violin_cohorts <- list(
    "GSE_All_Patients"          = clinical,
    "GSE_Treated"               = clinical %>% dplyr::filter(Chemo_adj == "Y"),
    "GSE_Stage3_4_Treated"      = clinical %>% dplyr::filter(TNM_stage %in% c(3, 4), Chemo_adj == "Y"),
    "GSE_Stage3_Treated"        = clinical %>% dplyr::filter(TNM_stage %in% c(3),    Chemo_adj == "Y"),
    "GSE_Stage3_4_Treated_MSS"  = clinical %>% dplyr::filter(TNM_stage %in% c(3, 4), Chemo_adj == "Y", MMR == "pMMR")
  )
  gse_km_cohorts <- list(
    "GSE_All_Treated"           = clinical %>% dplyr::filter(Chemo_adj == "Y"),
    "GSE_All_MSS"               = clinical %>% dplyr::filter(MMR == "pMMR"),
    "GSE_Stage2_Untreated_MSS"  = clinical %>% dplyr::filter(MMR == "pMMR", TNM_stage == 2, Chemo_adj == "N"),
    "GSE_Stage3_Treated_MSS"    = clinical %>% dplyr::filter(MMR == "pMMR", TNM_stage == 3, Chemo_adj == "Y"),
    "GSE_Stage34_Treated_MSS"   = clinical %>% dplyr::filter(MMR == "pMMR", TNM_stage %in% c(3, 4), Chemo_adj == "Y")
  )
  gse_cohorts <- c(gse_violin_cohorts, gse_km_cohorts,
                   subtype_cohorts(clinical, "CMS", "GSE"),
                   subtype_cohorts(clinical, "PDS", "GSE"))
  all_stats[["gse"]] <- run_survival_block(gse_cohorts, "GSE39582", score_cols)

  # ------------------------------------------------------------------------
  # PART B — TCGA-COAD (STAR TPM, primary-only)
  # ------------------------------------------------------------------------
  message("\n--- TCGA-COAD ---")
  query <- TCGAbiolinks::GDCquery(project = "TCGA-COAD",
                                  data.category = "Transcriptome Profiling",
                                  data.type     = "Gene Expression Quantification",
                                  workflow.type = "STAR - Counts",
                                  sample.type   = "Primary Tumor")
  # Two-level cache: SE (large) -> TPM matrix + colData (small).
  # colData is cached separately because MSI status is read from it downstream.
  mat_tcga <- cache_rds("TCGA_COAD_tpm", function() {
    se <- cache_rds("TCGA_COAD_se", function() {
      TCGAbiolinks::GDCdownload(query); TCGAbiolinks::GDCprepare(query)
    })
    m  <- prep_tcga_tpm(se)
    cd <- as.data.frame(SummarizedExperiment::colData(se))
    saveRDS(cd, file.path(CACHE_DIR, "TCGA_COAD_coldata.rds"))
    rm(se); gc()
    m
  })
  coad_coldata <- cache_rds("TCGA_COAD_coldata", function() {
    stop("colData cache missing — delete TCGA_COAD_tpm.rds and rerun to rebuild.")
  })

  scores_tcga <- run_ssgsea(mat_tcga, gs) %>%
    add_composite() %>%
    tibble::rownames_to_column("Full_Barcode") %>%
    dplyr::arrange(Full_Barcode) %>%             # deterministic vial selection
    dplyr::mutate(ID = tcga_patient_id(Full_Barcode)) %>%
    dplyr::distinct(ID, .keep_all = TRUE) %>%
    dplyr::select(-Full_Barcode)

  cdr <- load_tcga_cdr() %>% dplyr::filter(Project_ID == "TCGA-COAD")

  # Treatment status is CRC-specific (pan-cancer doesn't use it).
  has_value <- function(col, val) vapply(col, function(x) val %in% x, logical(1))
  treat <- TCGAbiolinks::GDCquery_clinic(project = "TCGA-COAD", type = "clinical") %>%
    dplyr::mutate(
      ID = trimws(toupper(gsub("\\.", "-", submitter_id))),
      Treatment_Status = dplyr::case_when(
        has_value(treatments_pharmaceutical_treatment_or_therapy, "yes") |
          has_value(treatments_radiation_treatment_or_therapy, "yes") ~ "Treated",
        has_value(treatments_pharmaceutical_treatment_or_therapy, "no") &
          has_value(treatments_radiation_treatment_or_therapy, "no")  ~ "Not Treated",
        TRUE ~ "Unknown/Not Reported")) %>%
    dplyr::select(ID, Treatment_Status) %>% dplyr::distinct(ID, .keep_all = TRUE)

  msi <- coad_coldata %>%
    dplyr::mutate(ID = substr(barcode, 1, 12)) %>%
    dplyr::arrange(barcode) %>% dplyr::distinct(ID, .keep_all = TRUE) %>%
    dplyr::select(ID, dplyr::any_of("paper_MSI_status"))

  clinical_tcga <- scores_tcga %>%
    dplyr::inner_join(cdr %>% dplyr::select(ID = Patient_ID, Stage, OS_event, OS_months, DFS_event, DFS_months),
                      by = "ID") %>%
    dplyr::left_join(treat, by = "ID") %>%
    dplyr::left_join(msi,   by = "ID") %>%
    dplyr::mutate(OS_delay = OS_months, Recurrence_event = DFS_event, Recurrence_delay = DFS_months) %>%
    derive_endpoints(os_event = "OS_event", os_time = "OS_delay",
                     rfs_event = "Recurrence_event", rfs_time = "Recurrence_delay")

  # --- Molecular subtyping (CMS + PDS), symbol-keyed, on the same SE ----------
  se_coad <- cache_rds("TCGA_COAD_se", function()
    stop("SE cache missing — delete TCGA_COAD_tpm.rds and rerun to rebuild."))
  emat_tcga_sym <- build_symbol_matrix(se_coad, "tcga")
  rm(se_coad); gc()
  cms_tcga <- cache_rds("cms_tcga", function() call_cms(emat_tcga_sym))
  pds_tcga <- cache_rds("pds_tcga", function() call_pds(emat_tcga_sym))
  clinical_tcga <- clinical_tcga %>%
    attach_subtypes(cms_tcga, clinical_key = "ID", key_fun = tcga_patient_id) %>%
    attach_subtypes(pds_tcga, clinical_key = "ID", key_fun = tcga_patient_id)
  write.csv(csv_safe(dplyr::select(clinical_tcga, ID, dplyr::any_of(c("CMS", "PDS")))),
            file.path(out_root, module, "TCGA_COAD_subtypes.csv"), row.names = FALSE)

  write.csv(csv_safe(clinical_tcga), file.path(out_root, module, "TCGA_COAD_clinical.csv"), row.names = FALSE)

  has_msi <- "paper_MSI_status" %in% colnames(clinical_tcga)
  tcga_cohorts <- list(
    "TCGA_All_Patients"     = clinical_tcga,
    "TCGA_Untreated"        = clinical_tcga %>% dplyr::filter(Treatment_Status == "Not Treated"),
    "TCGA_Treated"          = clinical_tcga %>% dplyr::filter(Treatment_Status == "Treated"),
    "TCGA_Stage1_Untreated" = clinical_tcga %>% dplyr::filter(Stage == "I",  Treatment_Status == "Not Treated"),
    "TCGA_Stage2_Untreated" = clinical_tcga %>% dplyr::filter(Stage == "II", Treatment_Status == "Not Treated")
  )
  if (has_msi) {
    tcga_cohorts[["TCGA_All_MSS"]] <- clinical_tcga %>% dplyr::filter(paper_MSI_status == "MSS")
    tcga_cohorts[["TCGA_All_MSI"]] <- clinical_tcga %>% dplyr::filter(paper_MSI_status %in% c("MSI-H", "MSI-L"))
  }
  tcga_cohorts <- c(tcga_cohorts,
                    subtype_cohorts(clinical_tcga, "CMS", "TCGA"),
                    subtype_cohorts(clinical_tcga, "PDS", "TCGA"))
  all_stats[["tcga"]] <- run_survival_block(tcga_cohorts, "TCGA-COAD", score_cols)

  # ------------------------------------------------------------------------
  # FDR (PROVISIONAL, #10) + plots + forest
  # ------------------------------------------------------------------------
  stats_df <- apply_fdr(do.call(rbind, all_stats),
                        by = c("Dataset", "Project", "Test", "Metric"))   # replicates Final.R
  write.csv(stats_df, file.path(out_root, module, "CRC_Statistical_Summary.csv"), row.names = FALSE)

  plot_survival_block(gse_cohorts,  "GSE39582",  score_cols, stats_df, out_root, module)
  plot_survival_block(tcga_cohorts, "TCGA-COAD", score_cols, stats_df, out_root, module)
  for (ds in c("GSE39582", "TCGA-COAD"))
    for (sc in score_cols)
      for (m in c("OS", "RFS"))
        generate_forest(stats_df, ds, m, sc, out_root, module)

  # ------------------------------------------------------------------------
  # Score-across-subtype characterisation (Kruskal-Wallis, DTP score by CMS/PDS)
  # ------------------------------------------------------------------------
  subtype_data <- list("GSE39582" = clinical, "TCGA-COAD" = clinical_tcga)
  sub_rows <- list()
  for (ds in names(subtype_data)) {
    cd <- subtype_data[[ds]]
    for (ax in c("CMS", "PDS")) {
      if (!ax %in% colnames(cd)) next
      for (sc in score_cols) {
        ks <- get_kruskal_stats(cd, sc, ax)
        sub_rows[[length(sub_rows) + 1]] <- data.frame(
          Dataset = ds, Subtype_Axis = ax, Score = sc,
          Raw_P = ks$p, Eps2 = ks$eps2, N = ks$n, K = ks$k,
          Is_Testable = !is.na(ks$p), stringsAsFactors = FALSE)
      }
    }
  }
  subtype_stats <- do.call(rbind, sub_rows)
  if (!is.null(subtype_stats)) {
    subtype_stats <- subtype_stats %>%
      dplyr::group_by(Dataset, Subtype_Axis) %>%
      dplyr::mutate(FDR_P = p.adjust(Raw_P, method = "BH")) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(Is_Significant = !is.na(FDR_P) & FDR_P < 0.05)
    write.csv(subtype_stats, file.path(out_root, module, "Subtype_Score_Stats.csv"), row.names = FALSE)
    for (i in seq_len(nrow(subtype_stats))) {
      r  <- subtype_stats[i, ]
      sc_lbl <- sub(SCORE_SUFFIX, "", r$Score)
      si <- list(testable = r$Is_Testable, sig = r$Is_Significant, fdr_p = r$FDR_P, eps2 = r$Eps2)
      generate_subtype_violin(subtype_data[[r$Dataset]], r$Score, r$Subtype_Axis,
        title    = sprintf("%s  %s by %s", r$Dataset, sc_lbl, r$Subtype_Axis),
        filename = sprintf("%s_%s_by_%s", r$Dataset, sc_lbl, r$Subtype_Axis),
        si = si, out_root = out_root, module = module)
    }
  }

  # ------------------------------------------------------------------------
  # Confounding + effect-modification of the DTP signature (CRC only)
  # Does CMS / PDS / Stage / MSI explain away (confound) or modify the score's
  # OS/RFS association? Adjusted Cox + interaction Cox on the core DTP scores.
  # ------------------------------------------------------------------------
  # Treated-Marisa is carried as an extra "cohort" (a GSE subset) so the subgroup
  # forest gets Treated-Marisa OS/RFS columns matching Figs 1-2. FDR is grouped by
  # Dataset, so these rows form their own family and leave the whole-cohort cells intact.
  cox_data    <- list("GSE39582" = harmonize_crc_modifiers(clinical, "GSE39582"),
                      "GSE39582 (treated)" = harmonize_crc_modifiers(dplyr::filter(clinical, Chemo_adj == "Y"), "GSE39582"),
                      "TCGA-COAD" = harmonize_crc_modifiers(clinical_tcga, "TCGA-COAD"))
  core_scores <- intersect(CORE_DTP_SCORES, score_cols)
  adj_specs   <- list(Clinicopath = c("Stage_bin", "MSI_group"),
                      CMS_adjusted = "CMS", PDS_adjusted = "PDS")

  adj_rows <- list(); int_rows <- list(); level_rows <- list()
  for (ds in names(cox_data)) {
    cd <- cox_data[[ds]]
    for (metric in c("OS", "RFS")) {
      mc <- metric_cols(metric); tcol <- mc[["t"]]; ecol <- mc[["e"]]
      for (sc in core_scores) {
        for (mlab in names(adj_specs)) {
          covs <- intersect(adj_specs[[mlab]], colnames(cd))
          a <- if (length(covs) == 0) NULL else get_cox_adjusted(cd, sc, tcol, ecol, covs)
          adj_rows[[length(adj_rows) + 1]] <- data.frame(
            Dataset = ds, Metric = metric, Score = sc, Model = mlab,
            HR = if (is.null(a)) NA else a$HR,
            HR_lower = if (is.null(a)) NA else a$HR_lower, HR_upper = if (is.null(a)) NA else a$HR_upper,
            Raw_P = if (is.null(a)) NA else a$P, C_index = if (is.null(a)) NA else a$C_index,
            LRT_score_P = if (is.null(a)) NA else a$LRT_P, PH_P = if (is.null(a)) NA else a$PH_P,
            Unadj_HR = if (is.null(a)) NA else a$Unadj_HR,
            Delta_logHR = if (is.null(a)) NA else a$Delta_logHR,
            Delta_logHR_pct = if (is.null(a)) NA else a$Delta_logHR_pct,
            N = if (is.null(a)) NA else a$N, N_events = if (is.null(a)) NA else a$N_events,
            Is_Testable = !is.null(a), stringsAsFactors = FALSE)
        }
        for (modf in intersect(CRC_MODIFIERS, colnames(cd))) {
          it <- get_cox_interaction(cd, sc, tcol, ecol, modf)
          int_rows[[length(int_rows) + 1]] <- data.frame(
            Dataset = ds, Metric = metric, Score = sc, Modifier = modf,
            Raw_P = if (is.null(it)) NA else it$Interaction_P, K = if (is.null(it)) NA else it$K,
            N = if (is.null(it)) NA else it$N, N_events = if (is.null(it)) NA else it$N_events,
            Is_Testable = !is.null(it), stringsAsFactors = FALSE)
          if (!is.null(it) && !is.null(it$per_level)) {
            pl <- it$per_level
            pl$Dataset <- ds; pl$Metric <- metric; pl$Score <- sc; pl$Modifier <- modf
            level_rows[[length(level_rows) + 1]] <- pl
          }
        }
      }
    }
  }

  adj_df <- apply_fdr(do.call(rbind, adj_rows), by = c("Dataset", "Metric", "Model"))
  write.csv(adj_df, file.path(out_root, module, "Adjusted_Cox_Summary.csv"), row.names = FALSE)
  int_df <- apply_fdr(do.call(rbind, int_rows), by = c("Dataset", "Metric", "Modifier"))
  write.csv(int_df, file.path(out_root, module, "Interaction_Cox_Summary.csv"), row.names = FALSE)

  level_df <- if (length(level_rows)) do.call(rbind, level_rows) else NULL
  if (!is.null(level_df)) {
    write.csv(level_df, file.path(out_root, module, "Subgroup_Score_HRs.csv"), row.names = FALSE)
    for (i in seq_len(nrow(int_df))) {
      r  <- int_df[i, ]
      hr <- level_df[level_df$Dataset == r$Dataset & level_df$Metric == r$Metric &
                     level_df$Score == r$Score & level_df$Modifier == r$Modifier, , drop = FALSE]
      if (nrow(hr) > 0)
        generate_subgroup_forest(hr, r$Dataset, r$Metric, r$Score, r$Modifier,
                                 r$FDR_P, out_root, module)
    }
  }

  # ---- Publication composite figures (groups 1-3) from this run's CSVs -------
  # Built from the tabular outputs just written by core/composite_figures.R
  # (sourced in main.R); guarded so a figure error never fails the analysis, but
  # any failure is warned so it surfaces in the end-of-run summary.
  if (exists("build_crc_composites")) {
    tryCatch(build_crc_composites(file.path(out_root, module)),
             error = function(e) {
               msg <- paste0("[composite figures] builder aborted: ", conditionMessage(e))
               message("  ", msg); warning(msg, call. = FALSE)
             })
  } else {
    warning("[composite figures] build_crc_composites() not found — is core/composite_figures.R sourced in main.R?",
            call. = FALSE)
  }

  message("CRC survival module complete.")
  invisible(stats_df)
}
