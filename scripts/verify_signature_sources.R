# verify_signature_sources.R
# Objectively identify the published source of each reference gene set in the
# panel by overlapping it against MSigDB (msigdbr). Run in the pipeline's R
# environment, where msigdbr is already a dependency:
#     Rscript verify_signature_sources.R
# If msigdbr is missing:  install.packages("msigdbr")
#
# For each reference set it prints (a) the top MSigDB sets by Jaccard and by
# raw intersection, and (b) overlaps with named candidate signatures matched by
# keyword, so the source can be confirmed by gene-list overlap rather than by
# name alone.

suppressMessages(library(msigdbr))

p <- read.csv("Data/Fatemeh.csv", check.names = FALSE, stringsAsFactors = FALSE)
names(p) <- trimws(gsub("﻿", "", names(p)))
sets  <- c("Myc", "CSC", "regStemCell", "revCSC", "Foetal", "IBD")
panel <- lapply(sets, function(s) {
  v <- trimws(p[[s]]); v <- v[v != "" & tolower(v) != "na"]; unique(v)
})
names(panel) <- sets

m <- msigdbr(species = "Homo sapiens")
gs_col  <- if ("gs_name" %in% names(m)) "gs_name" else grep("name",  names(m), value = TRUE)[1]
sym_col <- if ("gene_symbol" %in% names(m)) "gene_symbol" else grep("symbol", names(m), value = TRUE)[1]
cat_col <- intersect(c("gs_cat", "gs_collection", "gs_collection_name"), names(m))[1]
if (!is.na(cat_col)) m <- m[m[[cat_col]] %in% c("C2", "H", "C8"), ]  # curated + hallmark + cell-type (C8)

bysets <- lapply(split(m[[sym_col]], m[[gs_col]]), unique)
cat(sprintf("Scanned %d MSigDB C2+H sets.\n\n", length(bysets)))

# keyword patterns to surface named candidate signatures directly
kw <- c(Myc = "MYC", CSC = "MERLOS|INTESTINAL_STEM|STEM_CELL|EPHB2",
        regStemCell = "REGENERAT|FETAL|STEM_CELL", revCSC = "AYYAZ|REVIVAL|CLU",
        Foetal = "FETAL|MUSTATA|YAP|ONCOFETAL", IBD = "OSM|ONCOSTATIN|INFLAMMAT|WEST")

for (s in sets) {
  g <- panel[[s]]; n <- length(g)
  inter <- vapply(bysets, function(x) length(intersect(x, g)), integer(1))
  jac   <- vapply(bysets, function(x) length(intersect(x, g)) / length(union(x, g)), numeric(1))
  cat(sprintf("### %s  (panel n = %d)\n", s, n))
  cat("  top by Jaccard:\n")
  for (i in order(jac, decreasing = TRUE)[1:6])
    cat(sprintf("    %-58s  inter=%3d/%d  jac=%.3f\n", substr(names(bysets)[i], 1, 58), inter[i], n, jac[i]))
  cat("  top by intersection:\n")
  for (i in order(inter, decreasing = TRUE)[1:6])
    cat(sprintf("    %-58s  inter=%3d/%d  jac=%.3f\n", substr(names(bysets)[i], 1, 58), inter[i], n, jac[i]))
  hit <- grep(kw[[s]], names(bysets))
  if (length(hit)) {
    hit <- hit[order(inter[hit], decreasing = TRUE)][1:min(6, length(hit))]
    cat(sprintf("  keyword-matched candidates (/%s/):\n", kw[[s]]))
    for (i in hit)
      cat(sprintf("    %-58s  inter=%3d/%d  jac=%.3f\n", substr(names(bysets)[i], 1, 58), inter[i], n, jac[i]))
  }
  cat("\n")
}
