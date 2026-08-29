#!/usr/bin/env Rscript
# ==============================================================================
# figure_paired_violins_3way.R
# ------------------------------------------------------------------------------
# Single-plot 3-way paired signature score violin figure on one set of axes.
# Uses GSE50760 paired data across 3 tissue axes: Normal, Primary, Metastasis.
# Generates and updates signature_scores_paired.csv (18 patients x 3 tissues).
# All signatures are plotted on ONE single pair of axes with a shared Y-scale.
# ==============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(stats)
})

# Source core infrastructure helpers
source("core/config.R")
source("core/io.R")
source("core/expression.R")
source("core/id_conversion.R")
source("core/signatures.R")
source("core/scoring.R")
source("core/plotting.R")

message("\n==================== SINGLE-AXIS 3-WAY PAIRED VIOLIN PLOT ====================")

# 1. Load metadata & expression data for GSE50760 ------------------------------
message("--> Loading GSE50760 metadata and expression data...")
meta <- cache_rds("GSE50760_meta", function() GEOquery::getGEO("GSE50760", GSEMatrix = TRUE))
pdata  <- Biobase::pData(meta[[1]])
titles <- as.character(pdata$title)
gsm_ids <- rownames(pdata)

tissue <- dplyr::case_when(
  grepl("normal",  titles, ignore.case = TRUE) ~ "Normal",
  grepl("primary", titles, ignore.case = TRUE) ~ "Primary",
  grepl("metasta", titles, ignore.case = TRUE) ~ "Metastasis",
  TRUE ~ NA_character_
)
patient_code <- sub(".*?(AMC_[0-9]+).*", "\\1", titles)
patient_code[!grepl("^AMC_[0-9]+$", patient_code)] <- NA_character_
sample_meta <- data.frame(GSM = gsm_ids, Tissue = tissue, Patient = patient_code, stringsAsFactors = FALSE)

supp <- file.path(CACHE_DIR, "GSE50760_extracted")
all_files <- list.files(supp, full.names = TRUE)
files_to_read <- character()
for (id in sample_meta$GSM) {
  m <- all_files[grepl(paste0(id, "([._-]|$)"), all_files)]
  if (length(m)) files_to_read[id] <- m[1]
}

count_list <- list()
for (id in names(files_to_read)) {
  td <- read.delim(files_to_read[id], header = TRUE, stringsAsFactors = FALSE)[, c(1, 2)]
  colnames(td) <- c("gene_id", id)
  td$gene_id <- as.character(td$gene_id)
  td[[id]] <- suppressWarnings(as.numeric(td[[id]]))
  count_list[[id]] <- td %>%
    dplyr::group_by(gene_id) %>%
    dplyr::summarise(!!id := sum(.data[[id]], na.rm = TRUE), .groups = "drop") %>%
    as.data.frame()
}
fpkm_df <- Reduce(function(x, y) merge(x, y, by = "gene_id", all = TRUE), count_list)
fpkm <- as.matrix(fpkm_df[, setdiff(colnames(fpkm_df), "gene_id")])
rownames(fpkm) <- fpkm_df$gene_id
fpkm <- fpkm[, sample_meta$GSM]
if (anyNA(fpkm)) fpkm <- fpkm[rowSums(is.na(fpkm)) == 0, ]
log2_matrix <- log2(fpkm + 1)
log2_matrix <- harmonise_matrix_ids(log2_matrix)

# 2. Compute 3-way ssGSEA scores & save signature_scores_paired.csv -----------
message("--> Computing 3-way ssGSEA scores across Normal, Primary, Metastasis...")
panel <- load_signature_panel()
gs_collection <- build_gene_sets(panel)
scores <- run_ssgsea(log2_matrix, gs_collection)
scores <- add_composite(scores)

scores_df <- scores %>%
  tibble::rownames_to_column("GSM") %>%
  dplyr::left_join(sample_meta, by = "GSM") %>%
  dplyr::arrange(Patient, Tissue)
scores_df$Tissue <- factor(scores_df$Tissue, levels = c("Normal", "Primary", "Metastasis"))

# Write signature_scores_paired.csv
out_csv_dir <- "results/SingleSampleDiagnostics"
dir.create(out_csv_dir, recursive = TRUE, showWarnings = FALSE)
out_csv_path <- file.path(out_csv_dir, "signature_scores_paired.csv")
write.csv(scores_df, out_csv_path, row.names = FALSE)
message("--> Saved 3-way signature scores to: ", out_csv_path)

# Also check for newest run directory and save there if present
run_dirs <- list.dirs(".", recursive = FALSE)
run_dirs <- run_dirs[grepl("^\\./CRC_DTP_[0-9]{8}_[0-9]{4}$", run_dirs)]
if (length(run_dirs) > 0) {
  latest_run <- sort(run_dirs, decreasing = TRUE)[1]
  run_csv_path <- file.path(latest_run, "mets_de/results/SingleSampleDiagnostics/signature_scores_paired.csv")
  dir.create(dirname(run_csv_path), recursive = TRUE, showWarnings = FALSE)
  write.csv(scores_df, run_csv_path, row.names = FALSE)
  message("--> Saved copy to active run folder: ", run_csv_path)
}

# 3. Prepare dataset for single-axis plotting -----------------------------------
score_cols <- grep("_ssGSEA$", colnames(scores_df), value = TRUE)

long_df <- scores_df %>%
  tidyr::pivot_longer(
    cols = dplyr::all_of(score_cols),
    names_to = "Signature_Col",
    values_to = "Score"
  ) %>%
  dplyr::mutate(
    Signature = sub("_ssGSEA$", "", Signature_Col),
    Tissue    = factor(Tissue, levels = c("Normal", "Primary", "Metastasis"))
  )

sig_label_map <- c(
  "Up"        = "DTP Up",
  "Down"      = "DTP Down",
  "Composite" = "DTP Composite",
  "Fetal"     = "Foetal Intestinal",
  "RSC"       = "Regenerative Stem Cell",
  "revSC"     = "Revival revCSC",
  "IBD"       = "IBD Module",
  "CSC"       = "Columnar Stem Cell",
  "MYC"       = "MYC Module"
)
long_df$Sig_Label <- ifelse(long_df$Signature %in% names(sig_label_map),
                            sig_label_map[long_df$Signature], long_df$Signature)

sig_order <- c("DTP Up", "DTP Down", "DTP Composite", "Foetal Intestinal",
               "Regenerative Stem Cell", "Revival revCSC", "IBD Module",
               "Columnar Stem Cell", "MYC Module")
present_sigs <- intersect(sig_order, unique(long_df$Sig_Label))
remaining_sigs <- setdiff(unique(long_df$Sig_Label), present_sigs)
all_sig_levels <- c(present_sigs, remaining_sigs)
long_df$Sig_Label <- factor(long_df$Sig_Label, levels = all_sig_levels)

# Signature-specific color palette (each signature keeps its uniform color across all 3 tissues)
sig_colors <- c(
  "DTP Up"                 = "#2B5C8F",
  "DTP Down"               = "#3690C0",
  "DTP Composite"          = "#02818A",
  "Foetal Intestinal"      = "#67A9CF",
  "Regenerative Stem Cell" = "#D95F02",
  "Revival revCSC"         = "#7570B3",
  "IBD Module"             = "#A6761D",
  "Columnar Stem Cell"     = "#E7298A",
  "MYC Module"             = "#66A61E"
)

# Numeric coordinate mapping for single-axis layout with extra wide signature spacing
# Signature index i = 1..N spaced out by 1.6 units on X-axis
# Tissue offset: Normal = -0.35, Primary = 0.00, Metastasis = +0.35
sig_center_spacing <- 1.6
tissue_offsets <- c("Normal" = -0.35, "Primary" = 0.00, "Metastasis" = 0.35)
tissue_shapes  <- c("Normal" = 21, "Primary" = 24, "Metastasis" = 22) # Circle, Triangle, Square

long_df <- long_df %>%
  dplyr::mutate(
    Sig_Index = as.numeric(Sig_Label),
    Sig_Center = Sig_Index * sig_center_spacing,
    X_Coord    = Sig_Center + tissue_offsets[as.character(Tissue)]
  )

# Compute global Y range & extra wide expanded breaks
y_min_global <- min(long_df$Score, na.rm = TRUE)
y_max_global <- max(long_df$Score, na.rm = TRUE)
y_rng_global <- y_max_global - y_min_global

y_min_expanded <- floor((y_min_global - 0.28 * y_rng_global) * 5) / 5
y_max_expanded <- ceiling((y_max_global + 0.28 * y_rng_global) * 5) / 5

# Compute group mean trajectories per signature and tissue for cohort overlay
mean_df <- long_df %>%
  dplyr::group_by(Sig_Label, Tissue, Sig_Center) %>%
  dplyr::summarise(
    Mean_Score = mean(Score, na.rm = TRUE),
    X_Coord    = unique(X_Coord),
    .groups    = "drop"
  )

# Add sub-labels for Normal (N), Primary (P), Metastasis (M)
all_x_breaks <- c(
  as.vector(sapply(seq_along(all_sig_levels) * sig_center_spacing, function(xc) xc + tissue_offsets)),
  seq_along(all_sig_levels) * sig_center_spacing
)
all_x_labels <- c(
  rep(c("N", "P", "M"), length(all_sig_levels)),
  as.character(all_sig_levels)
)

message("--> Building single plot with extra wide expanded X/Y axis scale & statistical annotations...")

dividers <- seq_along(all_sig_levels)[-length(all_sig_levels)] * sig_center_spacing + (sig_center_spacing / 2)

p_single <- ggplot(long_df, aes(x = X_Coord, y = Score, group = interaction(Patient, Sig_Label))) +
  # Vertical background dividers separating signature blocks
  geom_vline(xintercept = dividers, color = "grey85", linetype = "dashed", linewidth = 0.55) +
  # Individual paired patient lines across Normal -> Primary -> Metastasis
  geom_line(color = "grey72", alpha = 0.40, linewidth = 0.45) +
  # Bold cohort mean trajectory line per signature
  geom_line(data = mean_df, aes(x = X_Coord, y = Mean_Score, group = Sig_Label, color = Sig_Label),
            linewidth = 1.15, inherit.aes = FALSE, alpha = 0.95) +
  # Violins per signature/tissue position
  geom_violin(aes(x = X_Coord, y = Score, fill = Sig_Label, group = interaction(Tissue, Sig_Label)),
              alpha = 0.25, color = NA, trim = FALSE, width = 0.28, inherit.aes = FALSE) +
  # Boxplots per signature/tissue position
  geom_boxplot(aes(x = X_Coord, y = Score, group = interaction(Tissue, Sig_Label)),
               fill = "white", color = "grey20", width = 0.12, outlier.shape = NA,
               alpha = 0.8, linewidth = 0.5, inherit.aes = FALSE) +
  # Patient data points with tissue shape encoding & signature color
  geom_point(aes(fill = Sig_Label, shape = Tissue), color = "grey15", size = 2.5, stroke = 0.5, alpha = 0.85) +
  scale_shape_manual(values = tissue_shapes, name = "Tissue Axis") +
  scale_fill_manual(values = sig_colors, guide = "none") +
  scale_color_manual(values = sig_colors, guide = "none") +
  scale_x_continuous(
    breaks = seq_along(all_sig_levels) * sig_center_spacing,
    labels = all_sig_levels,
    expand = expansion(mult = c(0.05, 0.05))
  ) +
  scale_y_continuous(
    breaks = seq(y_min_expanded, y_max_expanded, by = 0.2)
  ) +
  coord_cartesian(ylim = c(y_min_expanded, y_max_expanded)) +
  labs(
    title = "3-Way Paired Signature Activity Profiles Across Tissue Progression",
    subtitle = "GSE50760 Cohort (18 Paired Patients: Normal -> Primary CRC -> Liver Metastasis) | Shared Y-Scale & Cohort Means",
    x = "Gene Signature",
    y = "ssGSEA Score",
    caption = "Each signature displays paired patient trajectories (thin gray lines), cohort means (bold colored lines), and violins for Normal (circle), Primary (triangle), and Metastasis (square)."
  ) +
  pub_theme +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5, color = "grey10"),
    plot.subtitle = element_text(size = 11.5, face = "italic", hjust = 0.5, color = "grey30", margin = margin(b = 14)),
    axis.text.x = element_text(size = 11, face = "bold", angle = 20, hjust = 1, color = "grey15"),
    axis.text.y = element_text(size = 11, face = "bold"),
    axis.title.x = element_text(size = 13, face = "bold", margin = margin(t = 12)),
    axis.title.y = element_text(size = 13, face = "bold"),
    legend.position = "top",
    legend.title = element_text(size = 11.5, face = "bold"),
    legend.text = element_text(size = 10.5, face = "bold"),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.4),
    panel.grid.minor.y = element_blank()
  )

# 4. Save Single-Plot Figures --------------------------------------------------
fig_dir <- "figures"
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
png_path <- file.path(fig_dir, "Fig_paired_violins_3way_single_axis.png")
pdf_path <- file.path(fig_dir, "Fig_paired_violins_3way_single_axis.pdf")

ggsave(png_path, p_single, width = 22, height = 9, dpi = 300, bg = "white")
ggsave(pdf_path, p_single, width = 22, height = 9, bg = "white")

# Also update the primary figure paths
ggsave(file.path(fig_dir, "Fig_paired_violins_3way.png"), p_single, width = 22, height = 9, dpi = 300, bg = "white")
ggsave(file.path(fig_dir, "Fig_paired_violins_3way.pdf"), p_single, width = 22, height = 9, bg = "white")

message("\n--> Figures successfully written to:")
message("    [PNG] ", png_path)
message("    [PDF] ", pdf_path)
message("==============================================================================")
