# ==============================================================================
# core/gsea_plots.R  —  Custom GSEA enrichment plots.
#
# Provides plot_gsea_enrichment(), which replaces enrichplot::gseaplot2 with:
#   • A dashed ES = 0 reference line on the running-score panel.
#   • Direction annotations ("Up in X" / "Up in Y") at both ends of the rank
#     axis so the reader knows what the two poles mean.
#   • Multi-signature overlay: pass a vector of IDs to show all running-score
#     curves in one plot (used for the "all signatures" combined plot).
#   • Clean 3-panel layout via patchwork: running ES / hit marks / ranked stat.
#
# Depends on: pub_theme (core/plotting.R), patchwork, scales.
# ==============================================================================

# Extract running-ES data from a clusterProfiler gseaResult object.
# Uses the enrichplot internal helper gsInfo; falls back gracefully.
.gsea_running_data <- function(gsea_obj, gene_set_id) {
  tryCatch(
    enrichplot:::gsInfo(gsea_obj, gene_set_id),
    error = function(e)
      stop("Cannot extract running-ES data for gene set '", gene_set_id,
           "'. Verify the ID exists in the GSEA result.\nDetails: ", e$message)
  )
}

# ---------------------------------------------------------------------------
# plot_gsea_enrichment()
#
# gsea_obj      — clusterProfiler gseaResult object
# gene_set_ids  — character vector; one ID = single curve, many = overlay
# title         — plot title string
# label_left    — annotation placed at the LEFT of the rank axis
#                 (= genes ranked highest = most up in the POSITIVE direction
#                  of the contrast used to build the ranked list)
# label_right   — annotation placed at the RIGHT
# colors        — optional named/unnamed colour vector; auto-assigned if NULL
#
# Returns a patchwork object (invisible). Use ggplot2::ggsave() on the return
# value to save.
# ---------------------------------------------------------------------------
plot_gsea_enrichment <- function(gsea_obj,
                                  gene_set_ids,
                                  title      = "",
                                  label_left  = "Up in Comparison",
                                  label_right = "Up in Reference",
                                  colors      = NULL) {
  if (!requireNamespace("patchwork", quietly = TRUE))
    stop("patchwork is required for plot_gsea_enrichment(). ",
         "Install with: install.packages('patchwork')")

  n_sets <- length(gene_set_ids)

  # Assign colours
  if (is.null(colors)) {
    pal <- c("#E64B35", "#4DBBD5", "#00A087", "#3C5488",
             "#F39B7F", "#8491B4", "#91D1C2", "#DC0000")
    colors <- head(rep(pal, length.out = n_sets), n_sets)
  }
  colors <- setNames(rep_len(colors, n_sets), gene_set_ids)

  # Collect per-gene-set running-ES data frames
  es_list <- lapply(gene_set_ids, function(id) {
    df <- .gsea_running_data(gsea_obj, id)
    df$GeneSet <- id
    df
  })
  es_all <- dplyr::bind_rows(es_list)
  N      <- max(es_all$x)
  y_rng  <- range(es_all$runningScore, na.rm = TRUE)
  y_pad  <- diff(range(y_rng, 0)) * 0.13   # offset so labels don't overlap ES=0

  # ---- Panel 1: Running enrichment score --------------------------------
  p_es <- ggplot2::ggplot(es_all,
           ggplot2::aes(x = x, y = runningScore, colour = GeneSet)) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed",
                        colour = "grey45", linewidth = 0.55) +
    ggplot2::geom_line(linewidth = 0.9) +
    # Direction annotations at top of the axis
    ggplot2::annotate("text",
                      x     = N * 0.02,
                      y     = max(y_rng[2], 0) - y_pad,
                      label = label_left,
                      hjust = 0, size = 3, colour = "grey30", fontface = "italic") +
    ggplot2::annotate("text",
                      x     = N * 0.98,
                      y     = max(y_rng[2], 0) - y_pad,
                      label = label_right,
                      hjust = 1, size = 3, colour = "grey30", fontface = "italic") +
    ggplot2::scale_colour_manual(values = colors) +
    ggplot2::scale_x_continuous(expand = c(0, 0)) +
    ggplot2::labs(title  = title,
                  x      = NULL,
                  y      = "Enrichment Score",
                  colour = NULL) +
    pub_theme +
    ggplot2::theme(
      axis.text.x     = ggplot2::element_blank(),
      axis.ticks.x    = ggplot2::element_blank(),
      legend.position = if (n_sets > 1) "bottom" else "none"
    )

  # ---- Panel 2: Gene hit marks (one facet row per gene set) --------------
  hits <- dplyr::filter(es_all, .data$position == 1)
  p_hits <- ggplot2::ggplot(hits,
             ggplot2::aes(x = x, colour = GeneSet)) +
    ggplot2::geom_segment(ggplot2::aes(xend = x, y = 0, yend = 1),
                          linewidth = 0.35, alpha = 0.75) +
    ggplot2::facet_wrap(~ GeneSet, ncol = 1, strip.position = "left") +
    ggplot2::scale_colour_manual(values = colors, guide = "none") +
    ggplot2::scale_x_continuous(expand = c(0, 0)) +
    ggplot2::scale_y_continuous(breaks = NULL, expand = c(0, 0)) +
    ggplot2::labs(x = NULL, y = NULL) +
    pub_theme +
    ggplot2::theme(
      axis.text.x        = ggplot2::element_blank(),
      axis.ticks.x       = ggplot2::element_blank(),
      panel.spacing      = ggplot2::unit(0.5, "mm"),
      strip.text.y.left  = ggplot2::element_text(angle = 0, size = 7, face = "bold")
    )

  # ---- Panel 3: Ranked metric bar (same gene list for all sets) ----------
  stat_df <- es_list[[1]] %>%
    dplyr::mutate(fill = dplyr::if_else(.data$geneList >= 0, "pos", "neg"))
  p_stat <- ggplot2::ggplot(stat_df,
             ggplot2::aes(x = x, y = geneList, fill = fill)) +
    ggplot2::geom_col(width = 1) +
    ggplot2::scale_fill_manual(
      values = c(pos = "#E64B35", neg = "#4DBBD5"), guide = "none") +
    ggplot2::scale_x_continuous(expand = c(0, 0)) +
    ggplot2::labs(x = "Gene Rank", y = "Ranked\nMetric") +
    pub_theme

  # ---- Assemble with patchwork -------------------------------------------
  p_es / p_hits / p_stat +
    patchwork::plot_layout(heights = c(5, n_sets * 0.65 + 0.3, 2))
}

# ---------------------------------------------------------------------------
# .mets_nes_heatmap()  —  one-column NES heatmap (fill = NES, label = NES +
# FDR stars) over all signatures for a single GSEA contrast. Signatures are
# ordered top->bottom by `sig_order`; fill limits `lim` are shared across the
# two contrasts so the colours are directly comparable.
# ---------------------------------------------------------------------------
.mets_stars <- function(p) ifelse(is.na(p), "",
  ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", ""))))

.mets_nes_heatmap <- function(gsea_df, contrast_lab, pos_label, sig_order, lim) {
  d <- gsea_df
  d$Signature <- factor(d$ID, levels = rev(sig_order))
  d$Contrast  <- contrast_lab
  d$lab <- paste0(sprintf("%.2f", d$NES), .mets_stars(d$p.adjust))
  ggplot2::ggplot(d, ggplot2::aes(x = Contrast, y = Signature, fill = NES)) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.6) +
    ggplot2::geom_text(ggplot2::aes(label = lab), size = 3, fontface = "bold",
                       colour = "grey10") +
    ggplot2::scale_fill_gradient2(low = "#3B6EA5", mid = "white", high = "#E64B35",
      midpoint = 0, limits = c(-lim, lim), name = "NES") +
    ggplot2::labs(title = sprintf("NES + FDR (+NES = up in %s)", pos_label),
                  x = NULL, y = NULL) +
    pub_theme +
    ggplot2::theme(legend.position = "right",
                   panel.grid  = ggplot2::element_blank(),
                   axis.ticks  = ggplot2::element_blank(),
                   plot.title  = ggplot2::element_text(size = 11, face = "bold", hjust = 0))
}

# ---------------------------------------------------------------------------
# build_mets_gsea_composite()
#
# Four-panel composite for the mets GSEA analysis (thesis mets figure):
#   A  Normal-vs-Primary     all-signature enrichment plot
#   B  Normal-vs-Primary     NES heatmap (fill = NES, label = NES + FDR stars)
#   C  Primary-vs-Metastasis all-signature enrichment plot
#   D  Primary-vs-Metastasis NES heatmap
#
# Rebuildable from a completed run: reads the two saved gseaResult objects
# (results/gsea_{NormalVsPrimary,PrimaryVsMetastasis}.rds written by
# run_mets_de). NES / FDR for the heatmaps come from as.data.frame(obj).
# ---------------------------------------------------------------------------
build_mets_gsea_composite <- function(mets_dir,
                                      out_dir = file.path(mets_dir, "plots", "GSEA")) {
  np_rds <- file.path(mets_dir, "results", "gsea_NormalVsPrimary.rds")
  pm_rds <- file.path(mets_dir, "results", "gsea_PrimaryVsMetastasis.rds")
  if (!file.exists(np_rds) || !file.exists(pm_rds))
    stop("build_mets_gsea_composite(): missing gsea RDS in ",
         file.path(mets_dir, "results"), " (run run_mets_de first).")
  gsea_np <- readRDS(np_rds); gsea_pm <- readRDS(pm_rds)

  df_np <- as.data.frame(gsea_np); df_pm <- as.data.frame(gsea_pm)
  if (!nrow(df_np) || !nrow(df_pm))
    stop("build_mets_gsea_composite(): empty GSEA result(s).")

  # Shared signature order (top->bottom) = descending Primary-vs-Metastasis NES,
  # the headline contrast; any signature absent there is appended.
  ord       <- df_pm$ID[order(-df_pm$NES)]
  sig_order <- c(ord, setdiff(union(df_np$ID, df_pm$ID), ord))
  glim      <- max(abs(c(df_np$NES, df_pm$NES)), na.rm = TRUE)

  A <- plot_gsea_enrichment(gsea_np, gene_set_ids = sig_order,
         title = "Normal vs Primary: all signatures",
         label_left = "Up in Normal", label_right = "Up in Primary")
  B <- .mets_nes_heatmap(df_np, "N / P", pos_label = "Normal",     sig_order, glim)
  C <- plot_gsea_enrichment(gsea_pm, gene_set_ids = sig_order,
         title = "Primary vs Metastasis: all signatures",
         label_left = "Up in Metastasis", label_right = "Up in Primary")
  D <- .mets_nes_heatmap(df_pm, "P / M", pos_label = "Metastasis", sig_order, glim)

  composite <- patchwork::wrap_plots(
      patchwork::wrap_elements(A), patchwork::wrap_elements(B),
      patchwork::wrap_elements(C), patchwork::wrap_elements(D),
      ncol = 2, widths = c(3, 1)) +
    patchwork::plot_annotation(
      tag_levels = "A",
      title = "Mets GSEA: signature enrichment across the Normal to Primary to Metastasis axis",
      theme = ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", size = 15))) &
    ggplot2::theme(plot.tag = ggplot2::element_text(face = "bold", size = 16))

  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  n_sig <- length(sig_order); h <- max(14, 2 * (n_sig * 0.55 + 4))
  ggplot2::ggsave(file.path(out_dir, "Fig_Mets_GSEA_composite.pdf"), composite,
                  width = 16, height = h, limitsize = FALSE, device = grDevices::pdf)
  ggplot2::ggsave(file.path(out_dir, "Fig_Mets_GSEA_composite.png"), composite,
                  width = 16, height = h, dpi = 300, bg = "white", limitsize = FALSE)
  message("  wrote Fig_Mets_GSEA_composite (", n_sig, " signatures, ",
          round(h, 1), " in tall)")
  invisible(composite)
}
