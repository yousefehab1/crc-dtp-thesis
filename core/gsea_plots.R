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
# Uses the enrichplot internal helper gsInfo(). The ID is validated up front so
# a genuinely missing gene set and a broken/renamed gsInfo() (it is unexported,
# so a Bioconductor upgrade can remove it) produce different, accurate errors.
.gsea_running_data <- function(gsea_obj, gene_set_id) {
  if (!requireNamespace("enrichplot", quietly = TRUE))
    stop("enrichplot is required for plot_gsea_enrichment() (uses its internal ",
         "gsInfo() helper). Install with: BiocManager::install('enrichplot')")
  if (!gene_set_id %in% gsea_obj@result$ID)
    stop("Gene set '", gene_set_id, "' is not present in the GSEA result. ",
         "Available IDs: ",
         paste(utils::head(gsea_obj@result$ID, 20L), collapse = ", "))
  tryCatch(
    enrichplot:::gsInfo(gsea_obj, gene_set_id),
    error = function(e)
      stop("enrichplot:::gsInfo() failed for gene set '", gene_set_id,
           "' — the installed enrichplot version may have changed or removed ",
           "this unexported helper.\nDetails: ", conditionMessage(e),
           call. = FALSE)
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

  # Assign colours. Keyed to signature identity (alphabetical rank), NOT to the
  # order the caller happened to pass, so a given signature keeps its colour
  # across figures (the individual, combined, and composite plots pass the same
  # set in different orders). When there are more signatures than palette
  # swatches the palette is extended by a generated qualitative ramp, so colours
  # never silently recycle two signatures onto the same colour.
  if (is.null(colors)) {
    pal  <- c("#E64B35", "#4DBBD5", "#00A087", "#3C5488",
              "#F39B7F", "#8491B4", "#91D1C2", "#DC0000")
    base <- if (n_sets <= length(pal)) pal[seq_len(n_sets)]
            else grDevices::hcl.colors(n_sets, palette = "Dark 3")
    colors <- setNames(base[rank(gene_set_ids, ties.method = "first")], gene_set_ids)
  } else if (!is.null(names(colors)) && all(gene_set_ids %in% names(colors))) {
    colors <- colors[gene_set_ids]               # honour caller-supplied names
  } else {
    colors <- setNames(rep_len(colors, n_sets), gene_set_ids)
  }

  # Collect per-gene-set running-ES data frames
  es_list <- lapply(gene_set_ids, function(id) {
    df <- .gsea_running_data(gsea_obj, id)
    df$GeneSet <- id
    df
  })
  es_all <- dplyr::bind_rows(es_list)
  # Preserve the caller's order for the faceted hit panel and the legend — a
  # character column would otherwise facet/sort alphabetically, desyncing the
  # hit rows from the NES heatmaps drawn beside them in the composite.
  es_all$GeneSet <- factor(es_all$GeneSet, levels = gene_set_ids)
  N      <- max(es_all$x)
  # Pole labels sit in headroom ABOVE the curve envelope, so they never overprint
  # the running-score line (for a negative-NES set the curve peaks below zero and
  # the old data-max offset dropped the label straight onto the curve).
  y_hi   <- max(es_all$runningScore, 0, na.rm = TRUE)
  y_lo   <- min(es_all$runningScore, 0, na.rm = TRUE)
  y_span <- y_hi - y_lo
  y_lab  <- y_hi + y_span * 0.08

  # ---- Panel 1: Running enrichment score --------------------------------
  p_es <- ggplot2::ggplot(es_all,
           ggplot2::aes(x = x, y = runningScore, colour = GeneSet)) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed",
                        colour = "grey45", linewidth = 0.55) +
    ggplot2::geom_line(linewidth = 0.9) +
    # Direction annotations in the headroom above the curve
    ggplot2::annotate("text",
                      x     = N * 0.02,
                      y     = y_lab,
                      label = label_left,
                      hjust = 0, size = 3, colour = "grey30", fontface = "italic") +
    ggplot2::annotate("text",
                      x     = N * 0.98,
                      y     = y_lab,
                      label = label_right,
                      hjust = 1, size = 3, colour = "grey30", fontface = "italic") +
    ggplot2::expand_limits(y = y_lab + y_span * 0.04) +
    ggplot2::scale_colour_manual(values = colors) +
    ggplot2::scale_x_continuous(expand = c(0, 0), limits = c(1, N)) +
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
    ggplot2::scale_x_continuous(expand = c(0, 0), limits = c(1, N)) +
    ggplot2::scale_y_continuous(breaks = NULL, expand = c(0, 0)) +
    ggplot2::labs(x = NULL, y = NULL) +
    pub_theme +
    ggplot2::theme(
      axis.text.x        = ggplot2::element_blank(),
      axis.ticks.x       = ggplot2::element_blank(),
      panel.spacing      = ggplot2::unit(0.5, "mm"),
      strip.text.y.left  = ggplot2::element_text(angle = 0, size = 7, face = "bold")
    )

  # ---- Panel 3: Ranked metric (same gene list for all sets) --------------
  # Two filled areas (positive red / negative blue) rather than one geom_col bar
  # per gene: identical appearance, two polygons instead of ~N rectangles, and
  # no sub-pixel dropout where |metric| rounds below one device pixel.
  stat_df <- es_list[[1]]
  p_stat <- ggplot2::ggplot(stat_df, ggplot2::aes(x = x)) +
    ggplot2::geom_area(ggplot2::aes(y = pmax(.data$geneList, 0)), fill = "#E64B35") +
    ggplot2::geom_area(ggplot2::aes(y = pmin(.data$geneList, 0)), fill = "#4DBBD5") +
    ggplot2::scale_x_continuous(expand = c(0, 0), limits = c(1, N)) +
    ggplot2::labs(x = "Gene Rank", y = "Ranked\nMetric") +
    pub_theme

  # ---- Assemble with patchwork -------------------------------------------
  # guides = "collect" gathers the single colour legend to one place instead of
  # leaving it wedged between the ES curve and the hit-mark rows it annotates.
  assembled <- p_es / p_hits / p_stat +
    patchwork::plot_layout(heights = c(5, n_sets * 0.65 + 0.3, 2),
                           guides = "collect")
  if (n_sets > 1)
    assembled <- assembled & ggplot2::theme(legend.position = "bottom")
  assembled
}

# ---------------------------------------------------------------------------
# NOTE: .mets_nes_heatmap() and build_mets_gsea_composite() moved to
# modules/composites.R (the consolidated composite module). plot_gsea_enrichment()
# above stays here — it is a run-time primitive used by modules/mets_de.R.
# ---------------------------------------------------------------------------
