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
