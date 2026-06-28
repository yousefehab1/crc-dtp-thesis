# ==============================================================================
# core/clinical.R  —  CDR loader (#4) + ONE censoring-aware endpoint
# derivation (#5 labels, #6 logic) used by BOTH survival modules.
# ==============================================================================

# Load TCGA Clinical Data Resource (Liu et al. 2018). Returns one row per
# patient with OS + the selected recurrence endpoint (shared CDR_RFS_TYPE).
load_tcga_cdr <- function(path = CDR_FILE, rfs_type = CDR_RFS_TYPE) {
  stopifnot(file.exists(path))
  if (!rfs_type %in% c("OS", "DSS", "DFI", "PFI"))
    stop("CDR_RFS_TYPE must be one of OS/DSS/DFI/PFI")

  rfs_time_col <- paste0(rfs_type, ".time")
  cdr <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE,
                  na.strings = c("", "NA", "#N/A",
                                 "[Not Available]", "[Not Applicable]", "[Not Evaluated]"))
  colnames(cdr) <- trimws(gsub("^\uFEFF", "", colnames(cdr)))

  required <- c("bcr_patient_barcode", "type", "OS", "OS.time", rfs_type, rfs_time_col)
  miss <- setdiff(required, colnames(cdr))
  if (length(miss)) stop("CDR missing columns: ", paste(miss, collapse = ", "))

  has_stage <- "ajcc_pathologic_tumor_stage" %in% colnames(cdr)

  neg_os  <- sum(cdr$OS.time         < 0, na.rm = TRUE)
  neg_rfs <- sum(cdr[[rfs_time_col]] < 0, na.rm = TRUE)
  if (neg_os  > 0) warning(neg_os,  " patients with negative OS.time",       call. = FALSE)
  if (neg_rfs > 0) warning(neg_rfs, " patients with negative ", rfs_time_col, call. = FALSE)
  message(sprintf("CDR: %d patients, %d cancer types. Recurrence endpoint = %s (%d non-null).",
                  nrow(cdr), dplyr::n_distinct(cdr$type), rfs_type, sum(!is.na(cdr[[rfs_type]]))))

  cdr %>%
    dplyr::mutate(
      Patient_ID = toupper(trimws(bcr_patient_barcode)),
      Project_ID = paste0("TCGA-", toupper(trimws(type))),
      Stage      = if (has_stage)
        gsub("Stage |[ABC]$| ", "", ajcc_pathologic_tumor_stage) else NA_character_,
      OS_event   = as.integer(OS),
      OS_months  = as.numeric(OS.time) / DAYS_PER_MONTH,
      DFS_event  = as.integer(.data[[rfs_type]]),
      DFS_months = as.numeric(.data[[rfs_time_col]]) / DAYS_PER_MONTH
    ) %>%
    dplyr::filter(!is.na(Patient_ID), Patient_ID != "") %>%
    dplyr::distinct(Patient_ID, .keep_all = TRUE)
}

# Censoring-aware 3-year endpoints (#6) with corrected labels (#5).
# Expects numeric columns: <os_event>, <os_time>, <rfs_event>, <rfs_time> (months).
# Adds: Surv_3yr, Recurrence_3yr (factors) and OS3Y_/RFS3Y_ delay+event.
derive_endpoints <- function(df,
                             os_event = "OS_event",  os_time = "OS_delay",
                             rfs_event = "Recurrence_event", rfs_time = "Recurrence_delay",
                             os_cut = OS_CUTPOINT, rfs_cut = RFS_CUTPOINT) {
  df %>%
    dplyr::mutate(
      .oe = as.numeric(.data[[os_event]]),  .ot = as.numeric(.data[[os_time]]),
      .re = as.numeric(.data[[rfs_event]]), .rt = as.numeric(.data[[rfs_time]]),

      Surv_3yr = factor(dplyr::case_when(
        !is.na(.ot) & !is.na(.oe) & .ot <= os_cut & .oe == 1 ~ "Dead_3yr",
        !is.na(.ot) & !is.na(.oe)                            ~ "Alive_3yr",
        TRUE                                                 ~ NA_character_),
        levels = c("Dead_3yr", "Alive_3yr")),

      Recurrence_3yr = factor(dplyr::case_when(
        .re == 1 & .rt <= rfs_cut ~ "Recurred",
        .re == 1 & .rt >  rfs_cut ~ "Recurrence-Free",
        .re == 0 & .rt >= rfs_cut ~ "Recurrence-Free",
        TRUE                      ~ NA_character_),
        levels = c("Recurred", "Recurrence-Free")),

      OS3Y_event = dplyr::case_when(is.na(.oe) | is.na(.ot) ~ NA_real_,
                                    .ot > os_cut | .oe == 0 ~ 0, TRUE ~ 1),
      OS3Y_delay = ifelse(is.na(.ot), NA_real_, ifelse(.ot >= os_cut, os_cut, .ot)),

      RFS3Y_event = dplyr::case_when(is.na(.re) | is.na(.rt) ~ NA_real_,
                                     .rt > rfs_cut | .re == 0 ~ 0, TRUE ~ 1),
      RFS3Y_delay = ifelse(is.na(.rt), NA_real_, ifelse(.rt >= rfs_cut, rfs_cut, .rt))
    ) %>%
    dplyr::select(-.oe, -.ot, -.re, -.rt)
}
