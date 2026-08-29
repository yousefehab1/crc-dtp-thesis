# ==============================================================================
# modules/mets_de.R  —  GSE50760 paired DE, GSEA, PCA, single-sample ROC.
# Was Mets.R. This is a DIFFERENT analysis class (tissue-stage differential
# expression, no survival), so it shares INFRASTRUCTURE only: the single
# sig.csv loader (#1), the ssGSEA wrapper + naming (#12), the shared theme/
# palette (#13), timestamped output (#14), global seed (#15), GEO cache (#16).
# Its analytical core (limma paired DE, liver-purity covariate, GSEA, PCA, ROC)
# is preserved unchanged. FPKM duplicate SUMMING is kept on purpose (#3 — these
# are additive Cufflinks sub-features, not redundant measurements).
# ==============================================================================

run_mets_de <- function(out_root, panel, core_signature = "UP") {
  module <- "mets_de"
  message("\n==================== MODULE: METS DE ====================")
  base <- file.path(out_root, module)
  for (d in c("results/PurityAdjusted", "results/SingleSampleDiagnostics",
              "plots/GSEA/PrimaryVsMetastasis",
              "plots/GSEA/PrimaryVsMetastasis_purityAdjusted",
              "plots/GSEA/PrimaryVsNormal",
              "plots/PCA", "plots/LiverValidation", "plots/SingleSampleDiagnostics"))
    dir.create(file.path(base, d), recursive = TRUE, showWarnings = FALSE)

  # ---- Download (cached) + parse identity ------------------------------------
  meta <- cache_rds("GSE50760_meta", function() GEOquery::getGEO("GSE50760", GSEMatrix = TRUE))
  pdata  <- Biobase::pData(meta[[1]]); titles <- as.character(pdata$title); gsm_ids <- rownames(pdata)

  tissue <- dplyr::case_when(
    grepl("normal",  titles, ignore.case = TRUE) ~ "Normal",
    grepl("primary", titles, ignore.case = TRUE) ~ "Primary",
    grepl("metasta", titles, ignore.case = TRUE) ~ "Metastasis", TRUE ~ NA_character_)
  patient_code <- sub(".*?(AMC_[0-9]+).*", "\\1", titles)
  patient_code[!grepl("^AMC_[0-9]+$", patient_code)] <- NA_character_
  sample_meta <- data.frame(GSM = gsm_ids, Tissue = tissue, Patient = patient_code, stringsAsFactors = FALSE)
  if (anyNA(sample_meta$Tissue) || anyNA(sample_meta$Patient)) stop("Could not parse Tissue/Patient from titles.")
  if (!all(table(sample_meta$Patient, sample_meta$Tissue) == 1)) stop("Not a clean 18x3 design.")
  message(sprintf("Parsed %d samples / %d patients / %d tissues.",
                  nrow(sample_meta), dplyr::n_distinct(sample_meta$Patient), dplyr::n_distinct(sample_meta$Tissue)))

  # ---- Build FPKM matrix (sum duplicate symbols; #3) -------------------------
  supp <- cache_rds("GSE50760_supp_dir", function() {
    sf  <- GEOquery::getGEOSuppFiles("GSE50760", baseDir = CACHE_DIR)
    ex  <- file.path(CACHE_DIR, "GSE50760_extracted"); dir.create(ex, showWarnings = FALSE)
    untar(rownames(sf)[1], exdir = ex); ex
  })
  all_files <- list.files(supp, full.names = TRUE)
  files_to_read <- character()
  for (id in sample_meta$GSM) {
    m <- all_files[grepl(paste0(id, "([._-]|$)"), all_files)]
    if (length(m)) files_to_read[id] <- m[1]
  }
  if (length(setdiff(sample_meta$GSM, names(files_to_read))) > 0)
    stop("Missing FPKM file(s): ", paste(setdiff(sample_meta$GSM, names(files_to_read)), collapse = ", "))

  count_list <- list()
  for (id in names(files_to_read)) {
    td <- read.delim(files_to_read[id], header = TRUE, stringsAsFactors = FALSE)[, c(1, 2)]
    colnames(td) <- c("gene_id", id); td$gene_id <- as.character(td$gene_id)
    td[[id]] <- suppressWarnings(as.numeric(td[[id]]))
    count_list[[id]] <- td %>% dplyr::group_by(gene_id) %>%
      dplyr::summarise(!!id := sum(.data[[id]], na.rm = TRUE), .groups = "drop") %>% as.data.frame()
  }
  fpkm_df <- Reduce(function(x, y) merge(x, y, by = "gene_id", all = TRUE), count_list)
  fpkm <- as.matrix(fpkm_df[, setdiff(colnames(fpkm_df), "gene_id")]); rownames(fpkm) <- fpkm_df$gene_id
  fpkm <- fpkm[, sample_meta$GSM]
  if (anyNA(fpkm)) fpkm <- fpkm[rowSums(is.na(fpkm)) == 0, ]
  log2_matrix <- log2(fpkm + 1)                       # log2 convention, consistent with core

  # Harmonise row IDs to ID_TYPE. The actual ID type in the FPKM files is
  # auto-detected (Ensembl or symbol) and converted as needed.
  log2_matrix <- harmonise_matrix_ids(log2_matrix)

  # Drop undetected / invariant genes before any modelling. An all-zero (or
  # near-constant) FPKM row yields limma t = 0, and in the GSEA ranked list all
  # such genes are ordered purely by the random tie-breaker in .run_gsea() — a
  # large contiguous block of noise-ranked genes that dilutes the enrichment-
  # score normalisation and shows up as an empty band mid-plot. Require
  # detectable expression (FPKM > 0, i.e. log2 > 0) in at least ~10 % of samples.
  min_det     <- max(3L, floor(0.10 * ncol(log2_matrix)))
  detected    <- rowSums(log2_matrix > 0) >= min_det
  message(sprintf("  Expression filter: %d/%d genes retained (detected in >= %d samples).",
                  sum(detected), length(detected), min_det))
  log2_matrix <- log2_matrix[detected, , drop = FALSE]

  # ---- Paired DE (limma ~ Patient + Tissue) ----------------------------------
  cd <- sample_meta; rownames(cd) <- cd$GSM; cd <- cd[colnames(log2_matrix), ]
  cd$Tissue  <- factor(cd$Tissue, levels = c("Normal", "Primary", "Metastasis"))
  cd$Patient <- factor(cd$Patient)
  cd$Tissue_for_DE <- relevel(cd$Tissue, ref = "Primary")
  design  <- model.matrix(~ Patient + Tissue_for_DE, data = cd)
  fit     <- limma::eBayes(limma::lmFit(log2_matrix, design))
  coef_nm <- "Tissue_for_DEMetastasis"
  res <- limma::topTable(fit, coef = coef_nm, number = Inf, sort.by = "none") %>%
    tibble::rownames_to_column("gene") %>%
    dplyr::rename(stat = t, padj = adj.P.Val, pvalue = P.Value) %>% dplyr::filter(!is.na(stat))
  rownames(res) <- res$gene

  keep2 <- cd$Tissue %in% c("Primary", "Metastasis")
  log2_2way <- log2_matrix[, keep2]; cd2 <- droplevels(cd[keep2, ])

  # ---- Liver-purity covariate + negative-control validation ------------------
  # Cache the gene-set lookup (msigdbr loads ~100k rows just for one set).
  # ID_TYPE-independent: msigdbr always returns symbols; convert_gene_ids() below
  # handles ID_TYPE. If msigdbr is upgraded to a newer MSigDB release, delete
  # cache/msigdbr_hsiao_liver_specific_genes.rds manually to pick up the update
  # (same manual-invalidation convention as every other cache_rds() call).
  liver_genes_sym <- cache_rds("msigdbr_hsiao_liver_specific_genes", function() {
    msigdbr::msigdbr(species = "Homo sapiens") %>%
      dplyr::filter(gs_name == "HSIAO_LIVER_SPECIFIC_GENES") %>%
      dplyr::pull(gene_symbol) %>% unique()
  })
  # Convert liver gene symbols to ID_TYPE so they match the expression matrix.
  liver_genes <- if (ID_TYPE == "ensembl")
    convert_gene_ids(liver_genes_sym, from = "SYMBOL", to = "ENSEMBL")
  else
    liver_genes_sym
  liver_present <- intersect(liver_genes, rownames(log2_matrix))
  liver_gs <- GSEABase::GeneSetCollection(GSEABase::GeneSet(liver_present, setName = "Liver"))
  liver_scores <- GSVA::gsva(GSVA::ssgseaParam(exprData = log2_matrix, geneSets = liver_gs))
  liver_df <- as.data.frame(t(liver_scores)) %>% tibble::rownames_to_column("GSM") %>%
    dplyr::left_join(cd %>% dplyr::select(GSM, Tissue, Patient), by = "GSM") %>% dplyr::arrange(Patient, Tissue)
  write.csv(liver_df, file.path(base, "results/PurityAdjusted/liver_validation_scores_3way.csv"), row.names = FALSE)
  LiverScore <- as.numeric(liver_scores["Liver", colnames(log2_2way)]); cd2$LiverScore <- LiverScore

  paired_wilcox_p <- function(df, value_col, group_col, id_col, a, b) {
    w <- df %>% dplyr::filter(.data[[group_col]] %in% c(a, b)) %>%
      dplyr::select(dplyr::all_of(c(id_col, group_col, value_col))) %>%
      tidyr::pivot_wider(names_from = dplyr::all_of(group_col), values_from = dplyr::all_of(value_col))
    w <- w[stats::complete.cases(w[, c(a, b)]), ]
    stats::wilcox.test(w[[a]], w[[b]], paired = TRUE)$p.value
  }
  comp <- list(c("Normal", "Primary"), c("Primary", "Metastasis"), c("Normal", "Metastasis"))
  pvl  <- sapply(comp, function(x) paired_wilcox_p(liver_df, "Liver", "Tissue", "Patient", x[1], x[2]))
  p_traj <- ggplot2::ggplot(liver_df, ggplot2::aes(Tissue, Liver, group = Patient)) +
    ggplot2::geom_line(colour = "gray60", alpha = 0.6) +
    ggplot2::geom_point(ggplot2::aes(colour = Tissue), size = 3) +
    ggplot2::scale_colour_manual(values = pub_palette) + pub_theme +
    ggplot2::labs(title = NULL, x = NULL, y = "Liver score") +
    ggplot2::theme(legend.position = "none")
  ggplot2::ggsave(file.path(base, "plots/LiverValidation/Liver_Trajectory_3Way.png"), p_traj, width = 7, height = 6, dpi = 300)
  message(sprintf("Liver validation — N/P p=%.3g | P/M p=%.3g | N/M p=%.3g", pvl[1], pvl[2], pvl[3]))

  # ---- Purity sensitivity model + corrected matrix ---------------------------
  design_adj <- model.matrix(~ Patient + Tissue_for_DE + LiverScore, data = cd2)
  res_adj    <- limma::topTable(limma::eBayes(limma::lmFit(log2_2way, design_adj)), coef = coef_nm,
                                number = Inf, sort.by = "none")
  comparison_de <- res %>% dplyr::select(gene, logFC, stat, pvalue, padj) %>%
    dplyr::rename(logFC_unadj = logFC, t_unadj = stat, pval_unadj = pvalue, padj_unadj = padj) %>%
    dplyr::left_join(res_adj %>% tibble::rownames_to_column("gene") %>%
                       dplyr::select(gene, logFC, t, P.Value, adj.P.Val) %>%
                       dplyr::rename(logFC_adj = logFC, t_adj = t, pval_adj = P.Value, padj_adj = adj.P.Val),
                     by = "gene") %>%
    dplyr::mutate(Is_Liver_Gene = gene %in% liver_genes,
                  Status = dplyr::case_when(
                    padj_unadj < 0.05 & padj_adj >= 0.05 ~ "Lost significance after adjustment",
                    padj_unadj >= 0.05 & padj_adj < 0.05 ~ "Gained significance after adjustment",
                    padj_unadj < 0.05 & padj_adj < 0.05  ~ "Robust to adjustment",
                    TRUE ~ "Not significant either way")) %>% dplyr::arrange(padj_adj)
  write.csv(comparison_de, file.path(base, "results/PurityAdjusted/DE_unadjusted_vs_purityAdjusted.csv"), row.names = FALSE)
  preserve <- model.matrix(~ Patient + Tissue, data = cd2)
  log2_2way_corrected <- limma::removeBatchEffect(log2_2way, covariates = LiverScore, design = preserve)

  # ---- GSEA (clusterProfiler) on the sig.csv panel ---------------------------
  # Two comparisons share the same gene-set definitions; only the ranked list
  # differs.  custom_t2g is built once and reused.
  custom_t2g <- dplyr::bind_rows(lapply(names(panel),
                   function(nm) data.frame(Term = nm, Gene = panel[[nm]])))

  # Helper: rank a topTable result, run GSEA, save summary CSV, return result.
  .run_gsea <- function(de_res, csv_path) {
    rnk <- de_res$stat; names(rnk) <- de_res$gene
    rnk <- sort(rnk + runif(length(rnk), 1e-9, 1e-8), decreasing = TRUE)
    g <- clusterProfiler::GSEA(geneList = rnk, TERM2GENE = custom_t2g,
                                pvalueCutoff = 1, eps = 0,
                                minGSSize = 5, maxGSSize = 500,
                                BPPARAM = BiocParallel::SerialParam())
    write.csv(as.data.frame(g) %>% dplyr::arrange(p.adjust), csv_path, row.names = FALSE)
    g
  }

  # Helper: save individual + combined enrichment plots for one GSEA result.
  .save_gsea_plots <- function(gsea_obj, out_dir, label_left, label_right,
                               model = NULL) {
    tbl <- as.data.frame(gsea_obj) %>% dplyr::arrange(p.adjust)
    if (nrow(tbl) == 0) { message("  No GSEA results to plot."); return(invisible(NULL)) }

    all_ids <- tbl$ID
    # Model tag (e.g. "purity-adjusted") appended to every title so plots from
    # the adjusted and unadjusted rankings of the same contrast are never
    # confusable once separated from their directory.
    mtag <- if (!is.null(model)) paste0("  [", model, "]") else ""

    # Dotplot overview. Direction is a factor so the facets follow the
    # label_left/label_right convention rather than alphabetical order; each
    # direction is its own row (free + free-space y) so no empty rows appear.
    # p.adjust drives point size only — colouring by it too was redundant.
    tbl_dir <- tbl %>%
      dplyr::mutate(Direction = factor(ifelse(NES > 0, label_left, label_right),
                                       levels = c(label_left, label_right)))
    p_dot <- tbl_dir %>%
      dplyr::group_by(Direction) %>%
      dplyr::slice_min(p.adjust, n = 10, with_ties = FALSE) %>%
      dplyr::ungroup() %>%
      ggplot2::ggplot(ggplot2::aes(NES, reorder(Description, NES),
                                   size = -log10(p.adjust))) +
      ggplot2::geom_point(colour = "#3C5488") +
      ggplot2::facet_grid(Direction ~ ., scales = "free_y", space = "free_y") +
      ggplot2::scale_size_continuous(name = expression(-log[10] ~ italic(p)[adj])) +
      pub_theme +
      ggplot2::labs(x = "NES", y = NULL,
                    title = sprintf("GSEA overview  |  %s vs %s%s",
                                    label_left, label_right, mtag))
    ggplot2::ggsave(file.path(out_dir, "GSEA_Dotplot_Overview.png"),
                    p_dot, width = 10, height = 7, dpi = 300)

    # Individual enrichment plot per signature. Sanitised filenames are made
    # unique so two signatures differing only in punctuation (e.g. "DTP.Up" vs
    # "DTP-Up") cannot silently overwrite each other's PNG. tbl is arranged by
    # p.adjust, so row i corresponds to all_ids[i] / safe_ids[i].
    safe_ids <- make.unique(gsub("[[:punct:]]+", "_", all_ids), sep = "_")
    for (i in seq_along(all_ids)) {
      row <- tbl[i, ]
      le  <- if (!is.na(row$core_enrichment) && nzchar(row$core_enrichment))
               strsplit(row$core_enrichment, "/")[[1]] else character(0)
      ttl <- sprintf("%s (%s)%s\nNES=%.2f | p.adj=%.2e | leading edge=%d\n%s (left)  |  %s (right)",
                     row$Description,
                     ifelse(row$p.adjust < 0.05, "sig", "ns"), mtag,
                     row$NES, row$p.adjust, length(le),
                     label_left, label_right)
      # Single-signature plot uses the canonical enrichplot::gseaplot2() (the
      # standard clusterProfiler/enrichplot GSEA figure). Pole direction is
      # folded into the title since gseaplot2 has no end annotations. The
      # multi-curve combined + composite plots below keep plot_gsea_enrichment().
      p_ind <- enrichplot::gseaplot2(gsea_obj, geneSetID = row$ID,
                                     title = ttl, color = "#E64B35")
      ggplot2::ggsave(
        file.path(out_dir, sprintf("GSEA_%s.png", safe_ids[i])),
        p_ind, width = 8, height = 7, dpi = 300)
    }

    # Combined plot — all signatures overlaid in one panel
    p_all <- plot_gsea_enrichment(gsea_obj,
                                   gene_set_ids = all_ids,
                                   title = sprintf("All signatures  |  %s vs %s%s",
                                                   label_left, label_right, mtag),
                                   label_left  = label_left,
                                   label_right = label_right)
    ggplot2::ggsave(file.path(out_dir, "GSEA_AllSignatures_Combined.png"),
                    p_all,
                    width  = 10,
                    height = 5 + length(all_ids) * 0.4,
                    dpi    = 300)
    invisible(NULL)
  }

  # --- Comparison 1: Primary vs Metastasis (liver-UNADJUSTED) -----------------
  # Kept for transparency; the ADJUSTED ranking below is the primary figure.
  message("  Running GSEA: Primary vs Metastasis (liver-unadjusted)")
  gsea_pm <- .run_gsea(res,
                        file.path(base, "results/GSEA_PrimaryVsMetastasis_Summary.csv"))
  saveRDS(gsea_pm, file.path(base, "results/gsea_PrimaryVsMetastasis.rds"))  # transparency
  .save_gsea_plots(gsea_pm,
                   out_dir     = file.path(base, "plots/GSEA/PrimaryVsMetastasis"),
                   label_left  = "Up in Metastasis",
                   label_right = "Up in Primary",
                   model       = "liver-unadjusted")

  # --- Comparison 1b: Primary vs Metastasis, PURITY-ADJUSTED (PRIMARY figure) --
  # Metastasis samples carry heavy hepatic signal (liver-validation figure), so
  # the metastasis-direction enrichment is read off the liver-adjusted ranking
  # (res_adj: ~ Patient + Tissue_for_DE + LiverScore). res_adj is a raw topTable
  # (rownames = gene, stat column = t), so shape it like res before ranking.
  message("  Running GSEA: Primary vs Metastasis (purity-adjusted for liver score)")
  res_adj_ranked <- res_adj %>%
    tibble::rownames_to_column("gene") %>%
    dplyr::rename(stat = t, padj = adj.P.Val, pvalue = P.Value) %>%
    dplyr::filter(!is.na(stat))
  gsea_pm_adj <- .run_gsea(res_adj_ranked,
                            file.path(base, "results/GSEA_PrimaryVsMetastasis_purityAdjusted_Summary.csv"))
  saveRDS(gsea_pm_adj, file.path(base, "results/gsea_PrimaryVsMetastasis_adjusted.rds"))  # for the composite
  # Enrichment plots for the PRIMARY (purity-adjusted) ranking. The plots in
  # plots/GSEA/PrimaryVsMetastasis come from the UNADJUSTED ranking (kept only
  # for transparency); this adjusted set is what the thesis composite is built
  # from, so it gets its own directory and a "[purity-adjusted]" title tag.
  .save_gsea_plots(gsea_pm_adj,
                   out_dir     = file.path(base, "plots/GSEA/PrimaryVsMetastasis_purityAdjusted"),
                   label_left  = "Up in Metastasis",
                   label_right = "Up in Primary",
                   model       = "purity-adjusted")

  # Verification: adjusted NES / FDR for every signature.
  .adf <- as.data.frame(gsea_pm_adj) %>% dplyr::arrange(dplyr::desc(NES))
  message("  Primary vs Metastasis (purity-adjusted) — NES / FDR per signature:")
  for (i in seq_len(nrow(.adf)))
    message(sprintf("    %-16s NES=%6.3f  FDR=%.4f%s", .adf$ID[i], .adf$NES[i],
                    .adf$p.adjust[i], ifelse(.adf$p.adjust[i] < 0.05, " *", "")))

  # Supplementary: what liver adjustment changes for the P/M enrichment.
  .pm_u <- as.data.frame(gsea_pm)     %>% dplyr::select(ID, NES, p.adjust)
  .pm_a <- as.data.frame(gsea_pm_adj) %>% dplyr::select(ID, NES, p.adjust)
  pvm_cmp <- dplyr::full_join(
      dplyr::rename(.pm_u, NES_unadjusted = NES, FDR_unadjusted = p.adjust),
      dplyr::rename(.pm_a, NES_adjusted   = NES, FDR_adjusted   = p.adjust), by = "ID") %>%
    dplyr::mutate(dNES = NES_adjusted - NES_unadjusted,
                  Sig_unadjusted = FDR_unadjusted < 0.05,
                  Sig_adjusted   = FDR_adjusted   < 0.05) %>%
    dplyr::arrange(FDR_adjusted)
  write.csv(pvm_cmp, file.path(base, "results/GSEA_PvM_adjusted_vs_unadjusted.csv"), row.names = FALSE)

  # Small two-column NES heatmap: unadjusted vs adjusted (fill = NES, +FDR stars).
  .hm <- dplyr::bind_rows(
      dplyr::transmute(pvm_cmp, ID, Model = "Unadjusted", NES = NES_unadjusted, FDR = FDR_unadjusted),
      dplyr::transmute(pvm_cmp, ID, Model = "Adjusted",   NES = NES_adjusted,   FDR = FDR_adjusted)) %>%
    dplyr::mutate(Model = factor(Model, c("Unadjusted", "Adjusted")),
                  Signature = factor(ID, levels = rev(pvm_cmp$ID[order(pvm_cmp$NES_adjusted)])),
                  lab = paste0(sprintf("%.2f", NES), .sig_stars(FDR)))
  .lim <- max(abs(.hm$NES), na.rm = TRUE)
  p_cmp <- ggplot2::ggplot(.hm, ggplot2::aes(Model, Signature, fill = NES)) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.6) +
    ggplot2::geom_text(ggplot2::aes(label = lab), size = 3, fontface = "bold", colour = "grey10") +
    ggplot2::scale_fill_gradient2(low = "#3B6EA5", mid = "white", high = "#E64B35",
      midpoint = 0, limits = c(-.lim, .lim), name = "NES") +
    ggplot2::labs(x = "Primary vs Metastasis ranking", y = NULL,
                  caption = "NES (+FDR stars): which enrichments survive liver adjustment.") +
    pub_theme +
    ggplot2::theme(panel.grid = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank(),
                   plot.caption = ggplot2::element_text(hjust = 0, size = 8, colour = "grey30"))
  ggplot2::ggsave(file.path(base, "plots/GSEA/Fig_Mets_GSEA_PvM_adj_vs_unadj.png"),
                  p_cmp, width = 5.5, height = 5, dpi = 300, bg = "white")

  # --- Comparison 2: Primary vs Normal ----------------------------------------
  # Extracted from the same 3-way limma fit; coefficient Tissue_for_DENormal
  # gives Normal − Primary contrast by construction (Tissue_for_DE releveled to
  # ref = "Primary"), so negate the ranking statistic so positive = up in Primary.
  message("  Running GSEA: Primary vs Normal")
  res_pn <- limma::topTable(fit, coef = "Tissue_for_DENormal",
                             number = Inf, sort.by = "none") %>%
    tibble::rownames_to_column("gene") %>%
    dplyr::rename(stat = t, padj = adj.P.Val, pvalue = P.Value) %>%
    dplyr::filter(!is.na(stat)) %>%
    # Coefficient is Normal - Primary; negate so positive = up in PRIMARY.
    dplyr::mutate(stat = -stat, logFC = -logFC)
  gsea_pn <- .run_gsea(res_pn,
                        file.path(base, "results/GSEA_PrimaryVsNormal_Summary.csv"))
  saveRDS(gsea_pn, file.path(base, "results/gsea_PrimaryVsNormal.rds"))  # for the composite
  .save_gsea_plots(gsea_pn,
                   out_dir     = file.path(base, "plots/GSEA/PrimaryVsNormal"),
                   label_left  = "Up in Primary",
                   label_right = "Up in Normal")

  # The GSEA composite (build_mets_gsea_composite) and the liver-purity composite
  # are built separately by modules/composites.R (run_composites), decoupled from
  # this analysis run — see main.R. The per-signature enrichment plots above use
  # plot_gsea_enrichment() from core/gsea_plots.R and remain part of this run.

  # ---- PCA (3-way) -----------------------------------------------------------
  save_pca <- function(mat, meta, filename, title) {
    mv  <- mat[apply(mat, 1, var) > 0, ]
    pca <- prcomp(t(mv), scale. = TRUE); pv <- round(100 * pca$sdev^2 / sum(pca$sdev^2))
    pd  <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2], condition = meta$Tissue)
    p <- ggplot2::ggplot(pd, ggplot2::aes(PC1, PC2, colour = condition)) +
      ggplot2::geom_point(size = 4, alpha = 0.8) +
      ggplot2::stat_ellipse(ggplot2::aes(fill = condition), geom = "polygon", alpha = 0.1, show.legend = FALSE) +
      ggplot2::scale_colour_manual(values = pub_palette) + ggplot2::scale_fill_manual(values = pub_palette) +
      pub_theme + ggplot2::labs(title = title, x = paste0("PC1: ", pv[1], "%"), y = paste0("PC2: ", pv[2], "%"))
    ggplot2::ggsave(filename, p, width = 8, height = 6, dpi = 300)
  }
  save_pca(log2_matrix, cd, file.path(base, "plots/PCA/PCA_Global_3Way.png"), "Global PCA (3-Way)")
  if (core_signature %in% names(panel)) {
    cg <- intersect(panel[[core_signature]], rownames(log2_matrix))
    save_pca(log2_matrix[cg, ], cd, file.path(base, "plots/PCA/PCA_CoreSig_3Way.png"),
             sprintf("PCA: '%s' Signature (3-Way)", core_signature))
  }

  # ---- Single-sample diagnostics: ROC + paired violin (shared ssGSEA wrapper) -
  gs_collection <- build_gene_sets(panel)
  scores <- run_ssgsea(log2_2way_corrected, gs_collection)        # "<Sig>_ssGSEA" (#12)
  scores <- add_composite(scores)                                  # Composite = Up - Down
  scores_df <- scores %>% tibble::rownames_to_column("GSM") %>%
    dplyr::left_join(cd2 %>% dplyr::select(GSM, Tissue, Patient), by = "GSM") %>%
    dplyr::arrange(Patient, Tissue)
  scores_df$Tissue <- factor(scores_df$Tissue, levels = c("Primary", "Metastasis"))
  write.csv(scores_df, file.path(base, "results/SingleSampleDiagnostics/signature_scores_paired.csv"), row.names = FALSE)

  # Include Composite alongside individual signatures for ROC and Wilcoxon.
  sig_names <- c(names(panel),
                 if (paste0("Composite", SCORE_SUFFIX) %in% colnames(scores_df)) "Composite")

  roc_summary <- list()
  for (gs in sig_names) {
    col <- paste0(gs, SCORE_SUFFIX); if (!col %in% colnames(scores_df)) next
    v <- scores_df[[col]]
    mp <- mean(v[scores_df$Tissue == "Primary"]); mm <- mean(v[scores_df$Tissue == "Metastasis"])
    if (mm >= mp) { pos <- "Metastasis"; neg <- "Primary" } else { pos <- "Primary"; neg <- "Metastasis" }
    pred <- ROCR::prediction(v, scores_df$Tissue, label.ordering = c(neg, pos))
    auc  <- ROCR::performance(pred, "auc")@y.values[[1]]
    roc_summary[[gs]] <- data.frame(Signature = gs, AUC = auc, Predicts_For = pos)
    grDevices::png(file.path(base, "plots/SingleSampleDiagnostics", paste0(gs, "_ROC.png")), width = 1200, height = 1200, res = 200)
    plot(ROCR::performance(pred, "tpr", "fpr"), colorize = TRUE,
         main = sprintf("%s\nPrediction for %s: AUC=%.2f", gs, pos, auc))
    grDevices::dev.off()
  }
  write.csv(dplyr::bind_rows(roc_summary) %>% dplyr::arrange(dplyr::desc(AUC)),
            file.path(base, "results/SingleSampleDiagnostics/ROC_AUC_summary.csv"), row.names = FALSE)

  raw_p <- sapply(sig_names, function(gs) {
    col <- paste0(gs, SCORE_SUFFIX); if (!col %in% colnames(scores_df)) return(NA_real_)
    paired_wilcox_p(scores_df, col, "Tissue", "Patient", "Primary", "Metastasis")
  })
  adj_p <- p.adjust(raw_p, method = "BH")
  write.csv(data.frame(Signature = sig_names, p_raw = raw_p, p_adj_BH = adj_p) %>% dplyr::arrange(p_adj_BH),
            file.path(base, "results/SingleSampleDiagnostics/ssGSEA_paired_pvalues_BHadjusted.csv"), row.names = FALSE)

  message("Mets DE module complete.")
  invisible(list(res = res, comparison_de = comparison_de,
                 gsea_pm = as.data.frame(gsea_pm),
                 gsea_pm_adj = as.data.frame(gsea_pm_adj),
                 gsea_pn = as.data.frame(gsea_pn)))
}
