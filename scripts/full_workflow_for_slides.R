# Full SpatialOmicsMSI workflow for slide figures
#
# Run this script after installing/loading SpatialOmicsMSI.
# It creates CSV outputs and slide-ready PNG figures under:
# <SPATIALOMICS_RAT_BRAIN_DIR>/spatial_outputs

library(SpatialOmicsMSI)
library(ggplot2)
if (file.exists("R/msi_pipeline.R")) {
  source("R/msi_pipeline.R")
}

data_dir <- Sys.getenv("SPATIALOMICS_RAT_BRAIN_DIR", unset = file.path("data", "rat_brain_data"))
out_dir <- file.path(data_dir, "spatial_outputs")
plot_dir <- file.path(out_dir, "plots_for_slides")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

mzmine_file <- file.path(data_dir, "Filtered_7E4.csv")
imzml_file <- file.path(data_dir, "PNNL05A_V6b_CLMCAFAMM_Lipids_885.imzml")
vip_file <- Sys.getenv("SPATIALOMICS_VIP_FILE", unset = file.path(data_dir, "VIP.csv"))

message("1. Loading MZmine features and MSI data")
mzmine <- read.csv(mzmine_file, check.names = FALSE)

extracted <- load_msi_target_features(
  imzml_path = imzml_file,
  mzmine_features = mzmine,
  ppm = 10
)

pixel_matrix <- extracted$pixel_matrix
feature_mapping <- extracted$feature_mapping
raw_pixel_matrix <- pixel_matrix

message("2. Creating ion image thumbnails")
n_features <- length(feature_columns(pixel_matrix))
n_pages <- ceiling(n_features / 16)

for (page in seq_len(n_pages)) {
  p <- plot_feature_thumbnails(
    pixel_matrix,
    page = page,
    per_page = 16,
    transform = "log10"
  )
  ggsave(
    filename = file.path(plot_dir, sprintf("step2_ion_thumbnails_page_%02d.png", page)),
    plot = p,
    width = 12,
    height = 9,
    dpi = 300
  )
}

message("3. Cardinal spatial tissue/background segmentation")
background_result <- filter_background_cardinal_spatial_kmeans(
  pixel_matrix,
  k = 2,
  r = 1,
  ncomp = 10,
  weights = "adaptive",
  transform = "log10",
  foreground_cleanup = "largest_component",
  seed = 42
)
pixel_matrix <- background_result$matrix

write.csv(
  background_result$background_stats,
  file.path(out_dir, "cardinal_background_stats.csv"),
  row.names = FALSE
)
write.csv(
  data.frame(
    pixel_id = raw_pixel_matrix$pixel_id,
    x = raw_pixel_matrix$x,
    y = raw_pixel_matrix$y,
    cardinal_cluster = background_result$cluster,
    tissue_pixel = background_result$keep
  ),
  file.path(out_dir, "cardinal_background_mask.csv"),
  row.names = FALSE
)

message("4. Preprocessing for feature selection")
feature_selection <- preprocess_select_features(
  pixel_matrix,
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
preprocessed <- feature_selection$normalized_matrix
preprocessed_result <- feature_selection$preprocessing

message("5. Feature statistics and selection")
stats <- feature_stats(preprocessed)

p_cv <- ggplot(stats, aes(x = cv)) +
  geom_histogram(bins = 30, fill = "#3b82f6", color = "white") +
  labs(title = "Feature CV distribution", x = "CV", y = "Number of features") +
  theme_minimal()

ggsave(
  filename = file.path(plot_dir, "step3_feature_cv_distribution.png"),
  plot = p_cv,
  width = 8,
  height = 5,
  dpi = 300
)

selected <- feature_selection$selected_features
reduced <- feature_selection$reduced_matrix
if (nrow(selected) == 0 || length(feature_columns(reduced)) == 0) {
  stop(
    "No features were selected after Cardinal background filtering. ",
    "Try lowering nonzero_min or inspecting cardinal_background_mask.csv.",
    call. = FALSE
  )
}

distribution_rows <- lapply(names(preprocessed_result$distributions), function(step) {
  values <- preprocessed_result$distributions[[step]]
  values <- values[is.finite(values)]
  data.frame(step = step, intensity = values)
})
distribution_data <- do.call(rbind, distribution_rows)

p_pre <- ggplot(distribution_data, aes(x = intensity)) +
  geom_histogram(bins = 60, fill = "#14b8a6", color = "white") +
  facet_wrap(~step, scales = "free") +
  labs(title = "Preprocessing intensity distributions", x = "Value", y = "Count") +
  theme_minimal()

ggsave(
  filename = file.path(plot_dir, "step4_preprocessing_distributions.png"),
  plot = p_pre,
  width = 10,
  height = 7,
  dpi = 300
)

message("6. K-means spatial segmentation")
clustered_result <- cluster_pixels(reduced, k = 3, pca_components = 10)
clustered <- clustered_result$matrix

p_cluster <- ggplot(clustered, aes(x = x, y = y, fill = factor(cluster))) +
  geom_raster() +
  coord_fixed() +
  labs(title = "Molecular spatial segmentation", x = "x", y = "y", fill = "Cluster") +
  theme_minimal()

ggsave(
  filename = file.path(plot_dir, "step5_kmeans_segmentation.png"),
  plot = p_cluster,
  width = 7,
  height = 9,
  dpi = 300
)

diag <- cluster_diagnostics(reduced, max_k = 10, pca_components = 10)

p_elbow <- ggplot(diag, aes(x = k, y = tot_withinss)) +
  geom_line(color = "#3b82f6", linewidth = 1) +
  geom_point(color = "#3b82f6", size = 2) +
  labs(title = "K-means elbow plot", x = "k", y = "Total within-cluster SS") +
  theme_minimal()

ggsave(
  filename = file.path(plot_dir, "step5_elbow_plot.png"),
  plot = p_elbow,
  width = 7,
  height = 5,
  dpi = 300
)

message("7. Sub-region sampling")
cluster_roi_matrix <- define_rois(clustered, mode = "cluster")
sampled <- sample_subregions(
  cluster_roi_matrix,
  grid_size = 5,
  min_pixels = 30
)

p_samples <- ggplot(sampled$annotated_pixels, aes(x = x, y = y, fill = roi_id)) +
  geom_raster() +
  coord_fixed() +
  facet_wrap(~grid_cell) +
  labs(title = "ROI by grid-cell sub-regions", x = "x", y = "y", fill = "ROI") +
  theme_minimal() +
  theme(strip.text = element_text(size = 7))

ggsave(
  filename = file.path(plot_dir, "step6_subregion_sampling.png"),
  plot = p_samples,
  width = 12,
  height = 9,
  dpi = 300
)

message("8. MetaboAnalyst export")
metabo <- make_metaboanalyst_data(sampled$sample_matrix)

written <- write_pipeline_outputs(
  output_dir = out_dir,
  pixel_matrix = raw_pixel_matrix,
  feature_mapping = feature_mapping,
  background_stats = background_result$background_stats,
  selected_features = selected,
  reduced_matrix = reduced,
  preprocessed_matrix = preprocessed,
  clustered_matrix = clustered,
  sample_matrix = sampled$sample_matrix,
  sample_mapping = sampled$sample_mapping,
  metaboanalyst_data = metabo
)

message("CSV files written:")
message(paste(written, collapse = "\n"))

message("9-10. Optional MetaboAnalyst VIP import and back-mapping")
if (file.exists(vip_file)) {
  vip <- read.csv(vip_file, check.names = FALSE)
  vip <- normalize_metaboanalyst_result(vip, source_name = basename(vip_file))
  vip <- vip[order(vip$VIP, decreasing = TRUE), ]
  top_feature <- vip$Feature[1]

  p_vip_rank <- plot_metabo_feature_rank(vip, top_n = 20, source_name = basename(vip_file))
  ggsave(
    filename = file.path(plot_dir, "step9_vip_rank.png"),
    plot = p_vip_rank,
    width = 8,
    height = 6,
    dpi = 300
  )

  p_ion <- plot_ion_image(
    pixel_matrix,
    feature = top_feature,
    transform = "log10",
    title = paste0(top_feature, " | VIP=", round(vip$VIP[1], 2))
  )
  ggsave(
    filename = file.path(plot_dir, "step9_top_vip_ion_image.png"),
    plot = p_ion,
    width = 7,
    height = 9,
    dpi = 300
  )

  p_box <- plot_feature_region_boxplot(
    sample_matrix = sampled$sample_matrix,
    feature = top_feature,
    title = paste0(top_feature, " | VIP=", round(vip$VIP[1], 2))
  )
  ggsave(
    filename = file.path(plot_dir, "step9_top_vip_region_boxplot.png"),
    plot = p_box,
    width = 8,
    height = 6,
    dpi = 300
  )

  png(
    filename = file.path(plot_dir, "step9_vip_feature_view.png"),
    width = 3200,
    height = 1800,
    res = 300
  )
  plot_metabo_feature_view(
    result_data = vip,
    pixel_matrix = pixel_matrix,
    rank = 1,
    top_n = 20,
    source_name = basename(vip_file),
    transform = "log10"
  )
  dev.off()
} else {
  message("VIP file not found; skipped Step 8-9 VIP figures: ", vip_file)
}

message("Done.")
message("Slide figures are in: ", plot_dir)
message("MetaboAnalyst upload file is: ", file.path(out_dir, "metaboanalyst_data.csv"))
