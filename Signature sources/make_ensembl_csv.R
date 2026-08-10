#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# Turn gene_lists_symbols.csv into gene_lists_ensembl.csv:
# same column layout, Ensembl gene IDs instead of HGNC symbols.
#
# Usage
#   Rscript make_ensembl_csv.R
#       -> maps against the MSigDB Ensembl CHIP (downloads it; no packages needed)
#
#   Rscript make_ensembl_csv.R --gtf gencode.v44.annotation.gtf.gz
#       -> maps against your own annotation. USE THIS ONE if you have the GTF
#          your count matrix was quantified against: an Ensembl ID is only
#          meaningful against a specific release, and matching your own
#          annotation removes release-mismatch dropout.
#
#   Rscript make_ensembl_csv.R --source orgdb
#       -> uses Bioconductor org.Hs.eg.db instead (must already be installed)
#
# Other options
#   --csv  <file>   input  (default gene_lists_symbols.csv)
#   --out  <file>   output (default gene_lists_ensembl.csv)
#   --chip <url|file>
#
# Outputs
#   gene_lists_ensembl.csv          same headers, ENSG IDs, blank-padded
#   gene_lists_mapping_report.tsv   every symbol + status: OK | MULTI | UNMAPPED
#
# Base R only unless you choose --source orgdb.
# ---------------------------------------------------------------------------

DEFAULT_CHIP <- paste0(
  "https://data.broadinstitute.org/gsea-msigdb/msigdb/annotations/human/",
  "Human_Ensembl_Gene_ID_MSigDB.v2026.1.Hs.chip"
)

## ---- argument parsing -----------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
getopt <- function(flag, default = NULL) {
  i <- match(flag, args)
  if (is.na(i) || i == length(args)) default else args[[i + 1L]]
}
csv_in  <- getopt("--csv",    "gene_lists_symbols.csv")
csv_out <- getopt("--out",    "gene_lists_ensembl.csv")
chip    <- getopt("--chip",   DEFAULT_CHIP)
gtf     <- getopt("--gtf",    NULL)
source_ <- getopt("--source", NULL)

if (!file.exists(csv_in)) stop("Input not found: ", csv_in)

## ---- mapping sources ------------------------------------------------------
# Each returns a named list: symbol -> character vector of ENSG IDs.

split_map <- function(symbols, ids) {
  keep <- !is.na(symbols) & !is.na(ids) & nzchar(symbols) & nzchar(ids) &
          symbols != "---" & startsWith(ids, "ENSG")
  symbols <- symbols[keep]; ids <- ids[keep]
  lapply(split(ids, symbols), unique)
}

load_chip <- function(path) {
  local_path <- path
  if (grepl("^https?://", path)) {
    local_path <- tempfile(fileext = ".chip")
    message("Downloading CHIP ...")
    utils::download.file(path, local_path, quiet = TRUE, mode = "wb")
  }
  d <- utils::read.delim(local_path, stringsAsFactors = FALSE,
                         quote = "", comment.char = "")
  if (!all(c("Probe.Set.ID", "Gene.Symbol") %in% names(d)))
    stop("Not a CHIP file: ", path)
  # Probe Set ID is the Ensembl gene ID; non-ENSG rows are patch/scaffold entries
  split_map(d$Gene.Symbol, d$Probe.Set.ID)
}

load_gtf <- function(path) {
  con <- if (grepl("\\.gz$", path)) gzfile(path, "rt") else file(path, "rt")
  on.exit(close(con))
  lines <- readLines(con, warn = FALSE)
  lines <- lines[!startsWith(lines, "#")]
  f3 <- vapply(strsplit(lines, "\t", fixed = TRUE),
               function(x) if (length(x) >= 3L) x[[3]] else "", character(1))
  lines <- lines[f3 == "gene"]
  if (!length(lines)) stop("No 'gene' features found in ", path)
  attr9 <- vapply(strsplit(lines, "\t", fixed = TRUE), `[[`, character(1), 9L)
  gid   <- sub('.*gene_id "([^"]+)".*',   "\\1", attr9)
  gname <- sub('.*gene_name "([^"]+)".*', "\\1", attr9)
  gid   <- sub("\\..*$", "", gid)          # strip version suffix
  ok    <- gid != attr9 & gname != attr9   # sub() returns input when no match
  split_map(gname[ok], gid[ok])
}

load_orgdb <- function() {
  if (!requireNamespace("org.Hs.eg.db", quietly = TRUE))
    stop("org.Hs.eg.db is not installed. Use the default CHIP source, or:\n",
         '  BiocManager::install("org.Hs.eg.db")')
  m <- AnnotationDbi::select(org.Hs.eg.db::org.Hs.eg.db,
                             keys    = AnnotationDbi::keys(org.Hs.eg.db::org.Hs.eg.db, "SYMBOL"),
                             keytype = "SYMBOL",
                             columns = "ENSEMBL")
  split_map(m$SYMBOL, m$ENSEMBL)
}

mapping <- if (!is.null(gtf)) {
  message("Mapping source: ", gtf); load_gtf(gtf)
} else if (identical(source_, "orgdb")) {
  message("Mapping source: org.Hs.eg.db"); load_orgdb()
} else {
  message("Mapping source: ", chip); load_chip(chip)
}
message("  ", format(length(mapping), big.mark = ","), " symbols available")

## ---- map each column ------------------------------------------------------
dat <- utils::read.csv(csv_in, stringsAsFactors = FALSE, check.names = FALSE,
                       colClasses = "character")
set_names <- names(dat)

out_cols <- vector("list", length(set_names))
report   <- list()

for (j in seq_along(set_names)) {
  genes <- dat[[j]]
  genes <- genes[!is.na(genes) & nzchar(trimws(genes))]
  ids <- character(0)
  for (g in genes) {
    hit <- mapping[[g]]
    status <- if (is.null(hit) || !length(hit)) "UNMAPPED"
              else if (length(hit) > 1L)        "MULTI"
              else                              "OK"
    report[[length(report) + 1L]] <- data.frame(
      gene_set        = set_names[[j]],
      gene_symbol     = g,
      ensembl_gene_id = if (is.null(hit)) "" else paste(hit, collapse = ";"),
      status          = status,
      stringsAsFactors = FALSE
    )
    # keep the first ID; MULTI cases are listed in full in the report
    if (!is.null(hit) && length(hit) && !(hit[[1]] %in% ids)) ids <- c(ids, hit[[1]])
  }
  out_cols[[j]] <- ids
  message(sprintf("  %s: %d symbols -> %d ENSG", set_names[[j]], length(genes), length(ids)))
}

## ---- write ----------------------------------------------------------------
depth <- max(c(0L, vapply(out_cols, length, integer(1))))
padded <- lapply(out_cols, function(x) c(x, rep("", depth - length(x))))
res <- as.data.frame(padded, stringsAsFactors = FALSE)
names(res) <- set_names
utils::write.csv(res, csv_out, row.names = FALSE, quote = FALSE)

rep_df <- do.call(rbind, report)
utils::write.table(rep_df, "gene_lists_mapping_report.tsv", sep = "\t",
                   row.names = FALSE, quote = FALSE)

n <- table(factor(rep_df$status, levels = c("OK", "MULTI", "UNMAPPED")))
message(sprintf("\nOK %d   MULTI %d   UNMAPPED %d", n[["OK"]], n[["MULTI"]], n[["UNMAPPED"]]))
if (n[["UNMAPPED"]] > 0L) {
  message("Unmapped - check these before scoring:")
  u <- rep_df[rep_df$status == "UNMAPPED", c("gene_set", "gene_symbol")]
  for (i in seq_len(nrow(u))) message("  ", u$gene_set[i], "\t", u$gene_symbol[i])
}
message("\nWrote ", csv_out, " and gene_lists_mapping_report.tsv")
