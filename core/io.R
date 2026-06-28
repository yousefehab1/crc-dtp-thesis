# ==============================================================================
# core/io.R  —  Run scaffolding, caching, reproducibility (#14, #15, #16).
# ==============================================================================

# Start a run: set the global seed, create a timestamped output root so nothing
# is ever silently overwritten (the pan-cancer convention, applied everywhere).
init_run <- function(prefix = "run") {
  set.seed(GLOBAL_SEED)
  dir.create(CACHE_DIR, showWarnings = FALSE)
  root <- paste0(prefix, "_", format(Sys.time(), "%Y%m%d_%H%M"))
  dir.create(root, recursive = TRUE, showWarnings = FALSE)
  message("Run output root: ", root)
  root
}

# Capture exact package versions at the end of a run (#15).
finalize_run <- function(root) {
  writeLines(capture.output(sessionInfo()), file.path(root, "sessionInfo.txt"))
  message("sessionInfo saved to ", file.path(root, "sessionInfo.txt"))
}

# Generic on-disk cache (#16). `fun` is only evaluated on a cache miss, so
# slow GEO/GDC downloads happen once and resume across re-runs.
cache_rds <- function(key, fun) {
  f <- file.path(CACHE_DIR, paste0(gsub("[^A-Za-z0-9_]+", "_", key), ".rds"))
  if (file.exists(f)) {
    message("[cache] ", key, " (load)")
    return(readRDS(f))
  }
  message("[cache] ", key, " (build)")
  v <- fun()
  saveRDS(v, f)
  v
}

# Make a data frame CSV-safe (flatten list columns, strip newlines).
csv_safe <- function(df) {
  df %>%
    dplyr::mutate(dplyr::across(dplyr::where(is.list),
                                ~ sapply(., \(x) paste(unlist(x), collapse = " | ")))) %>%
    dplyr::mutate(dplyr::across(dplyr::where(is.character), ~ gsub("\n|\r", " ", .)))
}
