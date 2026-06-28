# ==============================================================================
# modules/pancan_survival.R  —  Pan-cancer DTP ssGSEA & survival across TCGA.
# Was Panncan.R. Now shares the core (loader/scoring/stats/plotting).
#
# Fixes applied: shares the single sig.csv loader (#1, selects Up/Down), log2
# ssGSEA + dropped back-transform (#2b), log-space collapse (#3), Dead/Alive
# labels (#5), censoring-aware endpoints (#6), continuous Cox (#7), effect size
# + 3-way routing (#8), gating (#9), primary-only — metastatic fallback removed
# (#17). FDR PROVISIONAL (#10): by=Family×Test (replicates Panncan.R).
# ==============================================================================

# Returns TRUE for errors that may resolve on retry (network, server-side).
# vctrs_error_subscript_oob and similar parsing errors are NOT transient.
.is_transient_gdc_error <- function(e) {
  cls <- class(e)
  if (any(grepl("vctrs_error|subscript|parse|type|coerce|match", cls, ignore.case = TRUE)))
    return(FALSE)
  msg <- tolower(paste(conditionMessage(e), collapse = " "))
  any(grepl("timeout|connection|network|curl|ssl|503|502|500|temporarily", msg))
}

# Build a minimal SummarizedExperiment directly from the downloaded STAR-Counts
# TSV files, bypassing GDCprepare entirely.  Used as a fallback when
# GDCprepare's clinical-annotation step crashes with vctrs_error_subscript_oob.
# The resulting SE has exactly what prep_tcga_tpm() needs:
#   - assay "tpm_unstrand"
#   - rowData: gene_id, gene_name
#   - colData: sample_type  (inferred from the TCGA barcode sample-code field)
.read_star_tpm_as_se <- function(query) {
  manifest <- as.data.frame(TCGAbiolinks::getResults(query))
  proj     <- manifest$project[1]
  gdc_dir  <- file.path("GDCdata", proj,
                        "Transcriptome_Profiling",
                        "Gene_Expression_Quantification")

  # Locate one TSV per file-UUID directory
  tsv_paths <- vapply(manifest$id, function(fid) {
    d <- file.path(gdc_dir, fid)
    f <- list.files(d, pattern = "\\.tsv$", full.names = TRUE)
    if (length(f) == 0) stop("No TSV found in: ", d)
    f[1]
  }, character(1))

  # Gene annotation from the first file (identical across all samples)
  gene_ann <- readr::read_tsv(tsv_paths[1], comment = "#", show_col_types = FALSE) %>%
    dplyr::filter(!grepl("^N_", gene_id)) %>%
    dplyr::select(gene_id, gene_name)

  # TPM matrix: genes × samples  (vapply → nrow(gene_ann) × n_samples)
  tpm_mat <- vapply(tsv_paths, function(f) {
    df <- readr::read_tsv(f, comment = "#", show_col_types = FALSE,
                          col_select = c("gene_id", "tpm_unstranded")) %>%
      dplyr::filter(!grepl("^N_", gene_id))
    df$tpm_unstranded[match(gene_ann$gene_id, df$gene_id)]
  }, numeric(nrow(gene_ann)))

  # Column names: prefer sample barcode, fall back to case/patient ID.
  # Guard: stop with a clear message if none of the expected columns exist.
  barcode_col <- intersect(c("sample.submitter_id", "cases.submitter_id", "cases"),
                           colnames(manifest))[1]
  if (is.na(barcode_col))
    stop("Cannot find a barcode column in GDC manifest. Available: ",
         paste(colnames(manifest), collapse = ", "))
  barcodes <- manifest[[barcode_col]]
  if (anyNA(barcodes) || length(barcodes) != ncol(tpm_mat))
    stop(sprintf("Barcode vector length (%d) != matrix columns (%d); NAs: %d",
                 length(barcodes), ncol(tpm_mat), sum(is.na(barcodes))))
  colnames(tpm_mat) <- barcodes

  # sample_type: use manifest column if present, else infer from barcode
  # TCGA barcode field 4 (1-based, 14-15 chars): "01"=Primary, "06"=Metastatic, "11"=Normal
  sample_type <- if ("sample_type" %in% colnames(manifest)) {
    manifest$sample_type
  } else {
    codes <- substr(barcodes, 14, 15)
    dplyr::case_when(
      codes == "01" ~ "Primary Tumor",
      codes == "06" ~ "Metastatic",
      codes == "11" ~ "Solid Tissue Normal",
      TRUE          ~ paste0("Sample_", codes)
    )
  }

  SummarizedExperiment::SummarizedExperiment(
    assays  = list(tpm_unstrand = tpm_mat),
    rowData = S4Vectors::DataFrame(gene_id   = gene_ann$gene_id,
                                   gene_name = gene_ann$gene_name,
                                   row.names = gene_ann$gene_id),
    colData = S4Vectors::DataFrame(sample_type = sample_type,
                                   row.names   = barcodes)
  )
}

# GDCprepare with direct-read fallback.
# If GDCprepare fails with vctrs_error_subscript_oob (TCGAbiolinks clinical
# annotation bug), we fall back to reading the already-downloaded TSV files
# directly.  All other errors are re-raised.
.gdc_prepare_safe <- function(query) {
  se <- tryCatch(TCGAbiolinks::GDCprepare(query), error = function(e) e)
  if (!inherits(se, "error")) return(se)

  is_clinical_err <- inherits(se, "vctrs_error_subscript_oob") ||
    any(grepl("vctrs_error_subscript", class(se), ignore.case = TRUE))
  if (!is_clinical_err) stop(se)

  message("   -> clinical annotation failed; reading TSV files directly")
  .read_star_tpm_as_se(query)
}

# Download + prepare one TCGA project with up to 3 retries on transient
# network/server errors.  Cached via cache_rds so it only runs on cache miss.
.gdc_with_retry <- function(query, proj, retries = 3, wait_sec = 15) {
  cache_rds(paste0("se_", proj), function() {
    last_err <- NULL
    for (attempt in seq_len(retries)) {
      result <- tryCatch({
        TCGAbiolinks::GDCdownload(query, method = "api", files.per.chunk = 6)
        .gdc_prepare_safe(query)
      }, error = function(e) e)
      if (!inherits(result, "error")) return(result)
      last_err <- result
      msg <- if (nzchar(trimws(last_err$message))) last_err$message else class(last_err)[1]
      # Don't retry non-transient errors (data/parsing failures won't resolve)
      if (!.is_transient_gdc_error(last_err)) {
        message(sprintf("   [!] non-transient error, skipping retries: %s", msg))
        break
      }
      if (attempt < retries) {
        message(sprintf("   [retry %d/%d] %s — waiting %ds", attempt, retries, msg, wait_sec))
        Sys.sleep(wait_sec)
      }
    }
    stop(last_err)
  })
}

# Max samples per GDCprepare call.  Cohorts larger than this are split into
# chunks so that peak RAM never exceeds ~chunk/total × full-cohort memory.
# 300 keeps BRCA (1230 samples) within ~4 GB per chunk on a 16 GB machine.
.PANCAN_CHUNK_SIZE <- 300L

# Build a TPM matrix for proj by chunking valid_patients into groups of
# chunk_size, running GDCprepare separately on each, and cbind-ing results.
# Files are already downloaded; only the in-memory SE construction is chunked.
.prep_tpm_chunked <- function(proj, valid_patients, chunk_size = .PANCAN_CHUNK_SIZE) {
  chunks <- split(valid_patients, ceiling(seq_along(valid_patients) / chunk_size))
  message(sprintf("   -> %d samples split into %d chunks of ~%d (memory optimisation)",
                  length(valid_patients), length(chunks), chunk_size))
  tpm_chunks <- lapply(seq_along(chunks), function(k) {
    bc <- chunks[[k]]
    message(sprintf("   -> chunk %d/%d  (%d samples)", k, length(chunks), length(bc)))
    q <- TCGAbiolinks::GDCquery(project       = proj,
                                data.category = "Transcriptome Profiling",
                                data.type     = "Gene Expression Quantification",
                                workflow.type = "STAR - Counts",
                                barcode       = bc)
    se <- .gdc_prepare_safe(q)
    m  <- prep_tcga_tpm(se)
    rm(se); gc()
    m
  })
  genes <- rownames(tpm_chunks[[1]])
  do.call(cbind, lapply(tpm_chunks, function(m) m[genes, , drop = FALSE]))
}

run_pancan_survival <- function(out_root, panel,
                                up_name = "Up", down_name = "Down") {
  module <- "pancan_survival"
  message("\n================== MODULE: PAN-CANCER SURVIVAL ==================")
  check_updown_overlap(panel, up_name, down_name)
  gs <- build_gene_sets(panel, which = c(up_name, down_name))   # only Up/Down (#1)
  up_genes <- panel[[up_name]]; down_genes <- panel[[down_name]]

  cohort_dir <- file.path(out_root, module, "cohort_csvs")
  dir.create(cohort_dir, recursive = TRUE, showWarnings = FALSE)

  # ---- Phase 1: CDR + project inventory --------------------------------------
  cdr_clinical <- load_tcga_cdr()
  cohort_summary <- cdr_clinical %>%
    dplyr::group_by(Project_ID) %>%
    dplyr::summarise(CDR_Patients = dplyr::n(),
                     OS_Available  = sum(!is.na(OS_months)  & !is.na(OS_event)),
                     DFS_Available = sum(!is.na(DFS_months) & !is.na(DFS_event)),
                     .groups = "drop") %>%
    dplyr::arrange(dplyr::desc(CDR_Patients))
  write.csv(cohort_summary, file.path(out_root, module, "cohort_summary.csv"), row.names = FALSE)

  valid_projects <- cohort_summary %>%
    dplyr::filter(!Project_ID %in% EXCLUDE_PROJECTS, CDR_Patients >= MIN_COHORT_N) %>%
    dplyr::pull(Project_ID)

  excluded_df <- cohort_summary %>%
    dplyr::mutate(
      Excluded = Project_ID %in% EXCLUDE_PROJECTS | CDR_Patients < MIN_COHORT_N,
      Reason   = dplyr::case_when(
        Project_ID %in% EXCLUDE_PROJECTS ~ paste0("Manually excluded (non-solid / DTP N/A): ",
                                                   paste(EXCLUDE_PROJECTS, collapse = ", ")),
        CDR_Patients < MIN_COHORT_N      ~ sprintf("Too few CDR patients (%d < threshold %d)",
                                                    CDR_Patients, MIN_COHORT_N),
        TRUE ~ NA_character_
      )
    ) %>%
    dplyr::filter(Excluded) %>%
    dplyr::select(Project_ID, CDR_Patients, Reason)

  if (nrow(excluded_df) > 0) {
    message(sprintf("\n--- Pre-queue exclusions (%d cancer type(s)) ---", nrow(excluded_df)))
    for (i in seq_len(nrow(excluded_df)))
      message(sprintf("  [EXCLUDED] %-15s  n=%4d  |  %s",
                      excluded_df$Project_ID[i], excluded_df$CDR_Patients[i], excluded_df$Reason[i]))
    message("---")
  }
  message(sprintf("%d projects queued.", length(valid_projects)))

  # ---- Phase 2: per-cohort RNA-seq -> ssGSEA -> clinical (checkpointed) -------
  cohort_outcomes <- vector("list", length(valid_projects))
  names(cohort_outcomes) <- valid_projects

  for (i in seq_along(valid_projects)) {
    proj     <- valid_projects[i]
    out_file <- file.path(cohort_dir, sprintf("clin_%s.csv", proj))
    message(sprintf("\n[%d/%d] %s", i, length(valid_projects), proj))

    if (file.exists(out_file)) {
      message("   -> already processed, skipping.")
      cohort_outcomes[[proj]] <- list(status = "Success (cached)", n = NA_integer_, reason = NA_character_)
      next
    }

    tryCatch({
      cdr_proj <- cdr_clinical %>% dplyr::filter(Project_ID == proj)
      valid_patients <- unique(cdr_proj$Patient_ID)
      if (length(valid_patients) == 0) {
        message("   -> no CDR patients.")
        cohort_outcomes[[proj]] <- list(status = "Skipped", n = 0L, reason = "No CDR patients")
        return(invisible(NULL))
      }

      query <- TCGAbiolinks::GDCquery(project = proj,
                                      data.category = "Transcriptome Profiling",
                                      data.type     = "Gene Expression Quantification",
                                      workflow.type = "STAR - Counts",
                                      barcode       = valid_patients)

      # Two-level cache: build TPM matrix once, reload the small RDS on reruns.
      # Large cohorts (> .PANCAN_CHUNK_SIZE) are processed in chunks so that
      # GDCprepare never loads more than chunk_size samples at once, keeping
      # peak RAM proportional to chunk_size rather than the full cohort.
      mat <- cache_rds(paste0("tpm_", proj), function() {
        if (length(valid_patients) > .PANCAN_CHUNK_SIZE) {
          .prep_tpm_chunked(proj, valid_patients)
        } else {
          se <- .gdc_with_retry(query, proj)
          m  <- prep_tcga_tpm(se)
          rm(se); gc()
          m
        }
      })

      pct_up <- gene_set_coverage(mat, up_genes); pct_down <- gene_set_coverage(mat, down_genes)
      message(sprintf("   -> coverage: Up=%.1f%%, Down=%.1f%%", pct_up, pct_down))
      if (pct_up < MIN_GENE_COV || pct_down < MIN_GENE_COV) {
        message("   -> coverage below threshold.")
        cohort_outcomes[[proj]] <- list(
          status = "Skipped", n = NA_integer_,
          reason = sprintf("Gene coverage too low (Up=%.1f%%, Down=%.1f%%; threshold=%.0f%%)",
                           pct_up, pct_down, MIN_GENE_COV))
        next
      }

      scores <- run_ssgsea(mat, gs)
      scores <- add_composite(scores, up_name, down_name)
      scores <- scores %>%
        tibble::rownames_to_column("Full_Barcode") %>%
        dplyr::arrange(Full_Barcode) %>%
        dplyr::mutate(Patient_ID = tcga_patient_id(Full_Barcode)) %>%
        dplyr::distinct(Patient_ID, .keep_all = TRUE) %>%
        dplyr::select(-Full_Barcode)

      clinical <- scores %>%
        dplyr::inner_join(cdr_proj %>% dplyr::select(Patient_ID, OS_event, OS_months, DFS_event, DFS_months),
                          by = "Patient_ID") %>%
        dplyr::mutate(Project_ID = proj)

      if (nrow(clinical) < MIN_COHORT_N) {
        message(sprintf("   -> %d matched < MIN_COHORT_N=%d.", nrow(clinical), MIN_COHORT_N))
        cohort_outcomes[[proj]] <- list(
          status = "Skipped", n = nrow(clinical),
          reason = sprintf("Only %d patients after CDR join (threshold=%d)", nrow(clinical), MIN_COHORT_N))
        next
      }
      write.csv(csv_safe(clinical), out_file, row.names = FALSE)
      message(sprintf("   -> saved %d patients.", nrow(clinical)))
      cohort_outcomes[[proj]] <- list(status = "Success", n = nrow(clinical), reason = NA_character_)
      rm(mat, scores, clinical); gc()
    }, error = function(e) {
      err_msg <- if (nzchar(trimws(e$message))) e$message else class(e)[1]
      message("   [!] failed: ", err_msg)
      cohort_outcomes[[proj]] <<- list(status = "Failed", n = NA_integer_, reason = err_msg)
    })
  }

  # ---- Phase 2 summary --------------------------------------------------------
  summary_df <- dplyr::bind_rows(lapply(names(cohort_outcomes), function(p) {
    o <- cohort_outcomes[[p]]
    data.frame(Project  = p,
               Status   = o$status,
               N        = if (is.null(o$n))      NA_integer_ else o$n,
               Reason   = if (is.null(o$reason)) NA_character_ else o$reason,
               stringsAsFactors = FALSE)
  }))
  write.csv(summary_df, file.path(out_root, module, "cohort_run_summary.csv"), row.names = FALSE)
  write.csv(excluded_df, file.path(out_root, module, "cohort_prequeu_exclusions.csv"), row.names = FALSE)

  n_ok   <- sum(grepl("^Success", summary_df$Status))
  n_skip <- sum(summary_df$Status == "Skipped")
  n_fail <- sum(summary_df$Status == "Failed")
  message(sprintf(
    "\n========== PAN-CANCER COHORT SUMMARY ==========\n  Succeeded : %d\n  Skipped   : %d\n  Failed    : %d\n  Pre-queue excluded: %d\n  Total queued: %d",
    n_ok, n_skip, n_fail, nrow(excluded_df), nrow(summary_df)))
  if (n_skip + n_fail > 0) {
    not_ok <- summary_df[summary_df$Status %in% c("Skipped", "Failed"), , drop = FALSE]
    max_proj <- max(nchar(not_ok$Project))
    message("")
    for (j in seq_len(nrow(not_ok))) {
      r <- not_ok[j, ]
      message(sprintf("  [%-*s]  %-8s  %s", max_proj, r$Project, r$Status, r$Reason))
    }
  }
  message("================================================\n")

  # ---- Phase 3: master compile, endpoints, batch correction ------------------
  message("\n--- compiling master ---")
  csv_files <- list.files(cohort_dir, pattern = "^clin_.*\\.csv$", full.names = TRUE)
  if (length(csv_files) == 0) stop("No per-cohort CSVs produced.")
  master <- dplyr::bind_rows(lapply(csv_files, \(f) readr::read_csv(f, col_types = readr::cols(.default = "c"))))

  up_col <- paste0(up_name, SCORE_SUFFIX); down_col <- paste0(down_name, SCORE_SUFFIX)
  comp_col <- paste0("Composite", SCORE_SUFFIX)
  master <- master %>%
    dplyr::mutate(dplyr::across(c(dplyr::all_of(c(up_col, down_col, comp_col)),
                                  OS_event, OS_months, DFS_event, DFS_months), as.numeric)) %>%
    dplyr::mutate(OS_delay = OS_months, Recurrence_event = DFS_event, Recurrence_delay = DFS_months) %>%
    derive_endpoints(os_event = "OS_event", os_time = "OS_delay",
                     rfs_event = "Recurrence_event", rfs_time = "Recurrence_delay")
  write.csv(master, file.path(out_root, module, "MASTER_Pan_Cancer_Clinical.csv"), row.names = FALSE)

  # event-rate sanity check per cohort
  event_rate <- master %>% dplyr::filter(!is.na(Surv_3yr)) %>%
    dplyr::group_by(Project_ID) %>%
    dplyr::summarise(N = dplyr::n(), Pct_Dead = round(mean(Surv_3yr == "Dead_3yr") * 100, 1), .groups = "drop") %>%
    dplyr::mutate(Flag = dplyr::case_when(Pct_Dead < 5  ~ "WARN <5% dead",
                                          Pct_Dead > 95 ~ "WARN >95% dead", TRUE ~ "OK")) %>%
    dplyr::arrange(Pct_Dead)
  write.csv(event_rate, file.path(out_root, module, "OS_Cutpoint_EventRates.csv"), row.names = FALSE)

  # batch correction across cancer types (PAN-CANCER ONLY — kept module-scoped)
  score_matrix <- t(as.matrix(master[, c(up_col, down_col, comp_col)]))
  corrected    <- limma::removeBatchEffect(score_matrix, batch = master$Project_ID)
  stopifnot(ncol(corrected) == nrow(master))
  corr_df <- as.data.frame(t(corrected))
  colnames(corr_df) <- paste0(c(up_name, down_name, "Composite"), SCORE_SUFFIX, "_Corrected")
  master <- dplyr::bind_cols(master, corr_df)
  message("Batch correction applied (corrected scores: pan-cancer aggregate only).")

  corr_cols   <- colnames(corr_df)
  uncorr_cols <- c(up_col, down_col, comp_col)

  # ---- Phase 4: statistics ---------------------------------------------------
  rows <- list()
  # Pan-cancer aggregate uses CORRECTED scores
  for (sc in corr_cols) {
    rows[[length(rows) + 1]] <- run_all_stats(master, sc, "OS",  "PanCancer", "Pan-Cancer", sc, "Surv_3yr")
    rows[[length(rows) + 1]] <- run_all_stats(master, sc, "RFS", "PanCancer", "Pan-Cancer", sc, "Recurrence_3yr")
  }
  # Per-cohort uses UNCORRECTED scores
  for (proj in unique(master$Project_ID)) {
    cd <- master %>% dplyr::filter(Project_ID == proj)
    for (sc in uncorr_cols) {
      rows[[length(rows) + 1]] <- run_all_stats(cd, sc, "OS",  "PanCancer", proj, sc, "Surv_3yr")
      rows[[length(rows) + 1]] <- run_all_stats(cd, sc, "RFS", "PanCancer", proj, sc, "Recurrence_3yr")
    }
  }
  stats_df <- do.call(rbind, rows) %>%
    dplyr::mutate(Family = ifelse(Project == "Pan-Cancer", "PanCancer", "PerCohort"))
  stats_df <- apply_fdr(stats_df, by = c("Family", "Test"))      # replicates Panncan.R (#10)
  write.csv(stats_df, file.path(out_root, module, "FDR_Stats_Summary.csv"), row.names = FALSE)

  cox_summary <- stats_df %>% dplyr::filter(Test == "Cox", !is.na(HR)) %>%
    dplyr::select(Project, Metric, Score, N, N_events, HR, HR_lower, HR_upper,
                  C_index, Raw_P, FDR_P, Is_Significant) %>%
    dplyr::arrange(Metric, Score, FDR_P)
  write.csv(cox_summary, file.path(out_root, module, "Cox_HR_Summary.csv"), row.names = FALSE)

  # ---- Phase 5: plots --------------------------------------------------------
  # Pan-cancer aggregate (corrected)
  agg <- setNames(list(master), "Pan-Cancer")
  plot_survival_block(agg, "PanCancer", corr_cols, stats_df, out_root, module)
  # Per-cohort (uncorrected)
  per <- lapply(unique(master$Project_ID), function(p) master %>% dplyr::filter(Project_ID == p))
  names(per) <- unique(master$Project_ID)
  plot_survival_block(per, "PanCancer", uncorr_cols, stats_df, out_root, module)
  # Forest plots across cohorts (per uncorrected score x metric), excluding aggregate
  for (sc in uncorr_cols)
    for (m in c("OS", "RFS"))
      generate_forest(stats_df, "PanCancer", m, sc, out_root, module, exclude_proj = "Pan-Cancer")

  message("Pan-cancer survival module complete.")
  invisible(stats_df)
}
