# Rat brain end-to-end workflow test
#
# R console usage:
# setwd("/Users/ly/Documents/Spatial Omics")
# Sys.setenv(SPATIALOMICS_RAT_BRAIN_DIR = "/Users/ly/Desktop/Jeff Xia/rat_brain_data")
# source("scripts/test_rat_brain_workflow.R")

library(ggplot2)

if (file.exists("R/msi_pipeline.R")) {
  source("R/msi_pipeline.R")
} else {
  library(SpatialOmicsMSI)
}

data_dir <- Sys.getenv("SPATIALOMICS_RAT_BRAIN_DIR", unset = file.path("data", "rat_brain_data"))
out_dir <- file.path(data_dir, "spatial_test_outputs")
plot_dir <- file.path(out_dir, "plots")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

mzmine_file <- file.path(data_dir, "Filtered_7E4.csv")
imzml_file <- file.path(data_dir, "PNNL05A_V6b_CLMCAFAMM_Lipids_885.imzml")

if (!file.exists(mzmine_file)) stop("Missing MZmine feature file: ", mzmine_file, call. = FALSE)
if (!file.exists(imzml_file)) stop("Missing imzML file: ", imzml_file, call. = FALSE)

message("1. Loading MZmine features and MSI data")
mzmine <- read.csv(mzmine_file, check.names = FALSE)
extracted <- load_msi_target_features(
  imzml_path = imzml_file,
  mzmine_features = mzmine,
  ppm = 10
)

raw_pixel_matrix <- extracted$pixel_matrix
feature_mapping <- extracted$feature_mapping
raw_fcols <- feature_columns(raw_pixel_matrix)
raw_tic <- rowSums(raw_pixel_matrix[raw_fcols], na.rm = TRUE)
raw_tic_plot_data <- data.frame(
  x = raw_pixel_matrix$x,
  y = raw_pixel_matrix$y,
  log10_tic = log10(raw_tic + 1)
)

p_raw_tic <- ggplot(raw_tic_plot_data, aes(x = x, y = y, fill = log10_tic)) +
  geom_raster() +
  coord_fixed() +
  scale_fill_viridis_c(option = "magma") +
  theme_minimal() +
  labs(title = "Raw log10 TIC", x = "x", y = "y", fill = "log10 TIC")

ggsave(file.path(plot_dir, "raw_log10_tic.png"), p_raw_tic, width = 7, height = 9, dpi = 300)

message("2. Cardinal spatial tissue/background segmentation")
background_result <- filter_background_cardinal_spatial_kmeans(
  raw_pixel_matrix,
  k = 2,
  r = 1,
  ncomp = 10,
  weights = "adaptive",
  transform = "log10",
  foreground_cleanup = "largest_component",
  seed = 42
)

tissue_matrix <- background_result$matrix
mask <- data.frame(
  pixel_id = raw_pixel_matrix$pixel_id,
  x = raw_pixel_matrix$x,
  y = raw_pixel_matrix$y,
  cardinal_cluster = background_result$cluster,
  tissue_pixel = background_result$keep
)

write.csv(background_result$background_stats, file.path(out_dir, "cardinal_background_stats.csv"), row.names = FALSE)
write.csv(mask, file.path(out_dir, "cardinal_background_mask.csv"), row.names = FALSE)

p_mask <- ggplot(mask, aes(x = x, y = y, color = tissue_pixel)) +
  geom_point(size = 0.25) +
  coord_equal() +
  scale_color_manual(values = c("FALSE" = "#bdbdbd", "TRUE" = "#d73027")) +
  theme_minimal() +
  labs(title = "Cardinal spatialKMeans tissue mask", color = "Tissue")

ggsave(file.path(plot_dir, "cardinal_tissue_mask.png"), p_mask, width = 7, height = 9, dpi = 300)

p_mask_overlay <- ggplot(raw_tic_plot_data, aes(x = x, y = y)) +
  geom_raster(aes(fill = log10_tic)) +
  geom_point(
    data = mask[mask$tissue_pixel, , drop = FALSE],
    aes(x = x, y = y),
    inherit.aes = FALSE,
    color = "#00ffff",
    size = 0.12,
    alpha = 0.6
  ) +
  coord_fixed() +
  scale_fill_viridis_c(option = "magma") +
  theme_minimal() +
  labs(title = "Raw log10 TIC with Cardinal tissue overlay", x = "x", y = "y", fill = "log10 TIC")

ggsave(file.path(plot_dir, "raw_tic_with_tissue_overlay.png"), p_mask_overlay, width = 7, height = 9, dpi = 300)

message("3. Preprocessing and feature selection")
feature_selection <- preprocess_select_features(
  tissue_matrix,
  feature_mapping = feature_mapping,
  do_background = FALSE,
  do_tic = TRUE,
  do_log = TRUE,
  cv_top_percent = 70,
  mean_min = 0,
  nonzero_min = 0.01,
  manual_columns = character(),
  combine_mode = "union"
)

selected <- feature_selection$selected_features
preprocessed <- feature_selection$normalized_matrix
reduced <- feature_selection$reduced_matrix

if (nrow(selected) == 0 || length(feature_columns(reduced)) == 0) {
  stop(
    "No features were selected after Cardinal background filtering. ",
    "Try lowering nonzero_min or inspecting cardinal_background_mask.csv.",
    call. = FALSE
  )
}

write.csv(feature_mapping, file.path(out_dir, "feature_mapping.csv"), row.names = FALSE)
write.csv(selected, file.path(out_dir, "selected_features.csv"), row.names = FALSE)
write.csv(preprocessed, file.path(out_dir, "preprocessed_matrix.csv"), row.names = FALSE)
write.csv(reduced, file.path(out_dir, "reduced_pixel_matrix.csv"), row.names = FALSE)

filtered_fcols <- feature_columns(tissue_matrix)
filtered_tic <- rowSums(tissue_matrix[filtered_fcols], na.rm = TRUE)
filtered_tic_plot_data <- data.frame(
  x = tissue_matrix$x,
  y = tissue_matrix$y,
  log10_tic = log10(filtered_tic + 1)
)

p_filtered_tic <- ggplot(filtered_tic_plot_data, aes(x = x, y = y, fill = log10_tic)) +
  geom_raster() +
  coord_fixed() +
  scale_fill_viridis_c(option = "magma") +
  theme_minimal() +
  labs(title = "Filtered tissue log10 TIC", x = "x", y = "y", fill = "log10 TIC")

ggsave(file.path(plot_dir, "filtered_tissue_log10_tic.png"), p_filtered_tic, width = 7, height = 9, dpi = 300)

stats <- feature_stats(preprocessed)
p_cv <- ggplot(stats, aes(x = cv)) +
  geom_histogram(bins = 30, fill = "#3b82f6", color = "white") +
  theme_minimal() +
  labs(title = "Feature CV distribution", x = "CV", y = "Number of features")

ggsave(file.path(plot_dir, "feature_cv_distribution.png"), p_cv, width = 8, height = 5, dpi = 300)

test_feature <- selected$column_name[1]
p_ion <- plot_ion_image(
  tissue_matrix,
  feature = test_feature,
  transform = "log10",
  title = paste("Example ion image:", test_feature)
)
ggsave(file.path(plot_dir, "test_ion_image.png"), p_ion, width = 7, height = 9, dpi = 300)

message("4. K-means spatial segmentation")
diagnostics <- cluster_diagnostics(
  reduced,
  max_k = 6,
  pca_components = 10,
  nstart = 10,
  seed = 1,
  max_subsample = 2000
)
write.csv(diagnostics, file.path(out_dir, "cluster_diagnostics.csv"), row.names = FALSE)

p_diagnostics <- ggplot(diagnostics, aes(x = k)) +
  geom_line(aes(y = tot_withinss), color = "#2563eb") +
  geom_point(aes(y = tot_withinss), color = "#2563eb") +
  theme_minimal() +
  labs(title = "Cluster diagnostics: elbow", x = "k", y = "Total within-cluster SS")

ggsave(file.path(plot_dir, "cluster_diagnostics.png"), p_diagnostics, width = 7, height = 5, dpi = 300)

clustered_result <- cluster_pixels(reduced, k = 3, pca_components = 10)
clustered <- clustered_result$matrix
write.csv(clustered, file.path(out_dir, "clustered_matrix.csv"), row.names = FALSE)

p_cluster <- ggplot(clustered, aes(x = x, y = y, fill = factor(cluster))) +
  geom_raster() +
  coord_fixed() +
  theme_minimal() +
  labs(title = "Molecular spatial segmentation", x = "x", y = "y", fill = "Cluster")

ggsave(file.path(plot_dir, "kmeans_segmentation.png"), p_cluster, width = 7, height = 9, dpi = 300)

message("5. Sub-region sampling and MetaboAnalyst export")
if (!exists("target_cluster", inherits = TRUE)) {
  target_cluster <- 1
}
target_roi_id <- paste0("roi_cluster_", target_cluster)
cluster_mapping <- data.frame(
  cluster = target_cluster,
  roi_id = target_roi_id,
  stringsAsFactors = FALSE
)

cluster_roi_matrix <- define_rois(
  clustered,
  mode = "cluster",
  cluster_mapping = cluster_mapping
)
sampled <- sample_subregions(
  cluster_roi_matrix,
  grid_size = 5,
  min_pixels = 30
)
metabo <- make_metaboanalyst_data(sampled$sample_matrix)

write.csv(sampled$sample_matrix, file.path(out_dir, "sample_matrix.csv"), row.names = FALSE)
write.csv(sampled$sample_mapping, file.path(out_dir, "sample_mapping.csv"), row.names = FALSE)
write.csv(metabo, file.path(out_dir, "metaboanalyst_data.csv"), row.names = FALSE)
write.csv(cluster_roi_matrix, file.path(out_dir, "cluster_roi_matrix.csv"), row.names = FALSE)
write.csv(cluster_mapping, file.path(out_dir, "cluster_roi_mapping.csv"), row.names = FALSE)

sample_fcols <- feature_columns(sampled$sample_matrix)
sample_pca <- stats::prcomp(sampled$sample_matrix[sample_fcols], center = TRUE, scale. = FALSE)
sample_scores <- data.frame(
  Sample = sampled$sample_matrix$sample_id,
  PC1 = sample_pca$x[, 1],
  PC2 = if (ncol(sample_pca$x) >= 2) sample_pca$x[, 2] else 0,
  stringsAsFactors = FALSE
)
write.csv(sample_scores, file.path(out_dir, "sample_pca_scores.csv"), row.names = FALSE)

sample_score_backmap <- backmap_sample_scores(sample_scores, sampled$sample_mapping, "PC1")
write.csv(sample_score_backmap, file.path(out_dir, "sample_score_backmap.csv"), row.names = FALSE)

p_sample_score_map <- plot_sample_score_map(
  sample_scores,
  sampled$sample_mapping,
  clustered,
  score_column = "PC1"
) +
  labs(title = paste0("Selected ROI sample PC1 back-map: ", target_roi_id))
ggsave(file.path(plot_dir, "sample_score_map.png"), p_sample_score_map, width = 7, height = 9, dpi = 300)

p_cluster_roi <- ggplot(clustered, aes(x = x, y = y)) +
  geom_raster(fill = "#eeeeee") +
  geom_tile(
    data = cluster_roi_matrix[!is.na(cluster_roi_matrix$roi_id), , drop = FALSE],
    aes(x = x, y = y, fill = roi_id),
    inherit.aes = FALSE
  ) +
  coord_fixed() +
  theme_minimal() +
  labs(title = "Selected cluster-defined ROI", x = "x", y = "y", fill = "ROI")

ggsave(file.path(plot_dir, "cluster_roi_definition.png"), p_cluster_roi, width = 7, height = 9, dpi = 300)

p_subregions <- ggplot(sampled$annotated_pixels, aes(x = x, y = y, fill = roi_id)) +
  geom_tile() +
  facet_wrap(~grid_cell, scales = "free") +
  theme_minimal() +
  theme(strip.text = element_text(size = 7)) +
  labs(title = "Grid-cell sampling QC", x = "x", y = "y", fill = "ROI")

ggsave(file.path(plot_dir, "subregion_sampling.png"), p_subregions, width = 12, height = 9, dpi = 300)

subregion_pixel_counts <- sampled$sample_mapping

p_subregion_counts <- ggplot(subregion_pixel_counts, aes(x = grid_cell, y = n_pixels, fill = roi_id)) +
  geom_col(position = "dodge") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Sub-region pixel counts", x = "Grid cell", y = "Pixels", fill = "ROI")

ggsave(file.path(plot_dir, "subregion_sample_counts.png"), p_subregion_counts, width = 10, height = 5, dpi = 300)

message("6. Coordinate ROI sampling test")
roi_metabo <- NULL

if (!exists("coordinate_roi_table", inherits = TRUE)) {
  x_limits <- stats::quantile(tissue_matrix$x, probs = c(0.35, 0.65), na.rm = TRUE)
  y_limits <- stats::quantile(tissue_matrix$y, probs = c(0.35, 0.65), na.rm = TRUE)
  coordinate_roi_table <- data.frame(
    roi_id = "center_roi",
    x_min = unname(x_limits[1]),
    x_max = unname(x_limits[2]),
    y_min = unname(y_limits[1]),
    y_max = unname(y_limits[2]),
    stringsAsFactors = FALSE
  )
}
if (!"roi_id" %in% names(coordinate_roi_table) && "roi_label" %in% names(coordinate_roi_table)) {
  names(coordinate_roi_table)[names(coordinate_roi_table) == "roi_label"] <- "roi_id"
}

roi_labeled_matrix <- define_rois(
  clustered,
  mode = "coordinate",
  roi_table = coordinate_roi_table
)
roi_sampled <- sample_subregions(
  roi_labeled_matrix,
  grid_size = 2,
  min_pixels = 30
)
roi_metabo <- make_metaboanalyst_data(roi_sampled$sample_matrix, group_column = "roi_id")

write.csv(coordinate_roi_table, file.path(out_dir, "coordinate_roi_table.csv"), row.names = FALSE)
write.csv(roi_sampled$sample_matrix, file.path(out_dir, "coordinate_roi_sample_matrix.csv"), row.names = FALSE)
write.csv(roi_sampled$sample_mapping, file.path(out_dir, "coordinate_roi_sample_mapping.csv"), row.names = FALSE)
write.csv(roi_metabo, file.path(out_dir, "coordinate_roi_metaboanalyst_data.csv"), row.names = FALSE)
write.csv(roi_labeled_matrix, file.path(out_dir, "coordinate_roi_matrix.csv"), row.names = FALSE)

coordinate_cluster_roi_matrix <- roi_labeled_matrix
write.csv(
  coordinate_cluster_roi_matrix,
  file.path(out_dir, "coordinate_cluster_roi_matrix.csv"),
  row.names = FALSE
)

p_coordinate_roi <- ggplot(coordinate_cluster_roi_matrix, aes(x = x, y = y)) +
  geom_tile(fill = "#eeeeee") +
  geom_tile(
    data = coordinate_cluster_roi_matrix[!is.na(coordinate_cluster_roi_matrix$roi_id), , drop = FALSE],
    aes(x = x, y = y, fill = factor(cluster)),
    inherit.aes = FALSE
  ) +
  coord_fixed() +
  theme_minimal() +
  labs(title = "Coordinate center ROI colored by cluster", x = "x", y = "y", fill = "Cluster")

ggsave(file.path(plot_dir, "coordinate_roi_definition.png"), p_coordinate_roi, width = 7, height = 9, dpi = 300)

p_roi <- ggplot(tissue_matrix, aes(x = x, y = y)) +
  geom_raster(fill = "#eeeeee") +
  geom_tile(
    data = roi_sampled$annotated_pixels,
    aes(x = x, y = y, fill = roi_id),
    inherit.aes = FALSE
  ) +
  coord_fixed() +
  theme_minimal() +
  labs(title = "Coordinate center ROI grid sampling QC", x = "x", y = "y", fill = "ROI")

ggsave(file.path(plot_dir, "coordinate_roi_sampling.png"), p_roi, width = 7, height = 9, dpi = 300)

p_roi_counts <- ggplot(roi_sampled$sample_mapping, aes(x = roi_id, y = n_pixels, fill = roi_id)) +
  geom_col() +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none") +
  labs(title = "Coordinate ROI pixel counts", x = "ROI", y = "Pixels")

ggsave(file.path(plot_dir, "coordinate_roi_pixel_counts.png"), p_roi_counts, width = 8, height = 5, dpi = 300)

cat("\n=== Rat brain workflow test complete ===\n")
cat("Pixels before background filtering:", nrow(raw_pixel_matrix), "\n")
cat("Pixels after Cardinal tissue filtering:", nrow(tissue_matrix), "\n")
cat("Selected features:", nrow(selected), "\n")
cat("Cluster sizes:\n")
print(table(clustered$cluster))
cat("Selected cluster ROI:", target_roi_id, "\n")
cat("MetaboAnalyst samples:", nrow(metabo), "\n")
if (!is.null(roi_metabo)) {
  cat("Coordinate ROI samples:", nrow(roi_metabo), "\n")
}
cat("Output directory:", normalizePath(out_dir), "\n")
cat("Plot directory:", normalizePath(plot_dir), "\n")
