# Create slide-ready pathway figures from existing pathway spatial outputs.
#
# R console usage:
# setwd("/Users/ly/Documents/Spatial Omics")
# source("scripts/make_pathway_slide_figures.R")

source("R/msi_pipeline.R")

if (!exists("data_dir", inherits = TRUE)) {
  data_dir <- "/Users/ly/Desktop/Jeff Xia/rat_brain_data"
}

test_out_dir <- file.path(data_dir, "spatial_test_outputs")
pathway_dir <- file.path(data_dir, "pathway_inputs")
spatial_dir <- file.path(pathway_dir, "pathway_spatial_maps")
slide_dir <- file.path(pathway_dir, "slide_pathway_figures")
dir.create(slide_dir, recursive = TRUE, showWarnings = FALSE)

required_files <- c(
  file.path(test_out_dir, "preprocessed_matrix.csv"),
  file.path(spatial_dir, "roi_scope_sample_matrix.csv"),
  file.path(spatial_dir, "roi_scope_sample_mapping.csv"),
  file.path(spatial_dir, "path_a_unsaturated_fatty_acids_known_contaminants_removed_linked_features.csv"),
  file.path(pathway_dir, "path_a_feature_chebi_stats.csv"),
  file.path(spatial_dir, "path_a_unsaturated_fatty_acids_known_contaminants_removed_score_map.png"),
  file.path(spatial_dir, "path_a_unsaturated_fatty_acids_known_contaminants_removed_feature_images.png")
)
missing <- required_files[!file.exists(required_files)]
if (length(missing) > 0) {
  stop("Missing required file(s):\n", paste(missing, collapse = "\n"), call. = FALSE)
}

message("1. Loading slide inputs")
pixel_matrix <- read.csv(file.path(test_out_dir, "preprocessed_matrix.csv"), check.names = FALSE, stringsAsFactors = FALSE)
sample_matrix <- read.csv(file.path(spatial_dir, "roi_scope_sample_matrix.csv"), check.names = FALSE, stringsAsFactors = FALSE)
sample_mapping <- read.csv(file.path(spatial_dir, "roi_scope_sample_mapping.csv"), check.names = FALSE, stringsAsFactors = FALSE)
clean_links <- read.csv(
  file.path(spatial_dir, "path_a_unsaturated_fatty_acids_known_contaminants_removed_linked_features.csv"),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
path_a_stats <- read.csv(file.path(pathway_dir, "path_a_feature_chebi_stats.csv"), check.names = FALSE, stringsAsFactors = FALSE)

clean_stats <- merge(
  clean_links,
  path_a_stats[, c("feature", "abbreviation", "signed_log2fc", "p.value", "FDR", "max_group", "min_group"), drop = FALSE],
  by = "feature",
  all.x = TRUE,
  sort = FALSE
)
clean_stats <- clean_stats[order(-abs(clean_stats$signed_log2fc)), , drop = FALSE]
key_feature <- clean_stats$feature[1]
key_label <- paste0(clean_stats$compound_name[1], " (", key_feature, ")")

message("2. Copying core pathway figures")
core_sources <- c(
  Path_A_clean_score_map = file.path(spatial_dir, "path_a_unsaturated_fatty_acids_known_contaminants_removed_score_map.png"),
  Path_A_clean_feature_images = file.path(spatial_dir, "path_a_unsaturated_fatty_acids_known_contaminants_removed_feature_images.png")
)
core_targets <- file.path(slide_dir, paste0(names(core_sources), ".png"))
file.copy(unname(core_sources), core_targets, overwrite = TRUE)

message("3. Plotting key feature region boxplot: ", key_feature)
p_box <- plot_feature_region_boxplot(
  sample_matrix,
  feature = key_feature,
  region_column = "roi_id",
  title = paste0("Region distribution: ", key_label)
) +
  ggplot2::labs(
    subtitle = paste0(
      "signed log2FC = ", signif(clean_stats$signed_log2fc[1], 3),
      "; FDR = ", signif(clean_stats$FDR[1], 3)
    ),
    y = paste0(key_feature, " sample intensity")
  ) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(size = 14),
    plot.subtitle = ggplot2::element_text(size = 10),
    axis.text.x = ggplot2::element_text(angle = 25, hjust = 1)
  )
boxplot_path <- file.path(slide_dir, paste0("Path_A_key_feature_region_boxplot_", key_feature, ".png"))
ggplot2::ggsave(boxplot_path, p_box, width = 7.5, height = 5.5, dpi = 300)

message("4. Plotting ROI/sample region map from pixel_id mapping")
mapped <- do.call(rbind, lapply(seq_len(nrow(sample_mapping)), function(i) {
  ids <- strsplit(as.character(sample_mapping$pixel_ids[i]), ",", fixed = TRUE)[[1]]
  data.frame(
    pixel_id = as.integer(ids),
    sample_id = sample_mapping$sample_id[i],
    roi_id = sample_mapping$roi_id[i],
    stringsAsFactors = FALSE
  )
}))
mapped <- mapped[!duplicated(mapped$pixel_id), , drop = FALSE]
region_data <- merge(
  pixel_matrix[, c("pixel_id", "x", "y"), drop = FALSE],
  mapped[, c("pixel_id", "roi_id"), drop = FALSE],
  by = "pixel_id",
  all.x = FALSE,
  sort = FALSE
)
x_steps <- diff(sort(unique(region_data$x)))
y_steps <- diff(sort(unique(region_data$y)))
tile_width <- stats::median(x_steps[x_steps > 0], na.rm = TRUE)
tile_height <- stats::median(y_steps[y_steps > 0], na.rm = TRUE)
if (!is.finite(tile_width)) tile_width <- 1
if (!is.finite(tile_height)) tile_height <- 1

p_regions <- ggplot2::ggplot(region_data, ggplot2::aes(x = .data[["x"]], y = .data[["y"]], fill = .data[["roi_id"]])) +
  ggplot2::geom_tile(width = tile_width, height = tile_height) +
  ggplot2::coord_fixed() +
  ggplot2::labs(
    title = "ROI labels used for pathway sampling",
    subtitle = "Back-mapped through pixel_id membership",
    x = "x",
    y = "y",
    fill = "ROI"
  ) +
  ggplot2::theme_minimal() +
  ggplot2::theme(
    plot.title = ggplot2::element_text(size = 15),
    plot.subtitle = ggplot2::element_text(size = 10)
  )
regions_path <- file.path(slide_dir, "ROI_sample_regions_from_pixel_mapping.png")
ggplot2::ggsave(regions_path, p_regions, width = 7, height = 7, dpi = 300)

message("5. Plotting optional rank + ion image view for the key feature")
rank_data <- path_a_stats[, c("feature", "FDR", "p.value", "signed_log2fc"), drop = FALSE]
names(rank_data)[names(rank_data) == "feature"] <- "Feature"
rank_path <- file.path(slide_dir, paste0("Path_A_key_feature_signed_log2fc_rank_ion_", key_feature, ".png"))
png(rank_path, width = 2400, height = 1200, res = 240)
plot_metabo_feature_view(
  rank_data,
  pixel_matrix = pixel_matrix,
  feature = key_feature,
  metric = "signed_log2fc",
  top_n = 12,
  transform = "log10"
)
dev.off()

manifest <- data.frame(
  figure = c(
    "Path_A_clean_score_map",
    "Path_A_clean_feature_images",
    "Path_A_key_feature_region_boxplot",
    "ROI_sample_regions_from_pixel_mapping",
    "Path_A_key_feature_rank_ion"
  ),
  path = c(
    core_targets[1],
    core_targets[2],
    boxplot_path,
    regions_path,
    rank_path
  ),
  note = c(
    "Clean Path A ROI-level pathway score map after known contaminant removal.",
    "Ion images for the 4 clean linked features in rno01040.",
    paste0("Key feature chosen by largest abs(signed_log2fc): ", key_label),
    "QC/method figure showing ROI labels mapped through pixel_id membership.",
    paste0("Optional rank plus ion image view for ", key_label)
  ),
  stringsAsFactors = FALSE
)
write.csv(manifest, file.path(slide_dir, "slide_pathway_figures_manifest.csv"), row.names = FALSE)
write.csv(clean_stats, file.path(slide_dir, "Path_A_clean_linked_feature_stats.csv"), row.names = FALSE)

cat("\n=== Slide pathway figures complete ===\n")
cat("Selected key feature:", key_label, "\n")
cat("Output directory:", normalizePath(slide_dir), "\n")
print(manifest)
