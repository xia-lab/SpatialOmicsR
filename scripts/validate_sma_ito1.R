source("R/msi_pipeline.R")

data_root <- Sys.getenv("SPATIALOMICS_ITO1_DIR", unset = file.path("data", "ITO1"))
output_dir <- "outputs/sma_ito1_validation"

paths <- list.files(
  data_root,
  pattern = "csv$",
  recursive = TRUE,
  full.names = TRUE
)

if (length(paths) != 3) {
  stop("Expected 3 ITO1 CSV files under: ", data_root, call. = FALSE)
}

paths <- paths[order(basename(paths))]

series <- load_peakpicked_msi_series(
  paths,
  section_ids = c("ITO1_A1", "ITO1_B1", "ITO1_C1"),
  sample_ids = c("mPD1", "mPD3", "mPD4")
)

feature_selection <- preprocess_select_features(
  series$pixel_matrix,
  feature_mapping = series$feature_mapping,
  serial = TRUE,
  min_section_fraction = 1,
  cv_top_percent = 1,
  nonzero_min = 0.05,
  do_background = FALSE,
  do_tic = TRUE,
  do_log = TRUE
)

selected <- feature_selection$selected_features
reduced <- feature_selection$reduced_matrix

# Temporary region labels for code validation only.
# Replace this with anatomical, registered, or manually annotated regions for real analysis.
labels <- data.frame(
  pixel_id = reduced$pixel_id,
  matched_region_label = ifelse(reduced$x <= stats::median(reduced$x), "Left", "Right")
)

labeled <- apply_matched_region_labels(reduced, labels)
clustered <- cluster_pixels(
  labeled,
  k = 3,
  pca_components = 10
)
matched <- sample_matched_regions(labeled, min_pixels = 100)
metabo <- make_metaboanalyst_data(
  matched$sample_matrix,
  group_column = "matched_region_label"
)

write_pipeline_outputs(
  output_dir,
  pixel_matrix = series$pixel_matrix,
  feature_mapping = series$feature_mapping,
  section_mapping = series$section_mapping,
  background_stats = feature_selection$background_stats,
  selected_features = selected,
  reduced_matrix = labeled,
  preprocessed_matrix = labeled,
  clustered_matrix = clustered$matrix,
  sample_matrix = matched$sample_matrix,
  sample_mapping = matched$sample_mapping,
  metaboanalyst_data = metabo
)

plot_paths <- write_validation_plots(
  output_dir,
  pixel_matrix = series$pixel_matrix,
  sample_matrix = matched$sample_matrix,
  sample_mapping = matched$sample_mapping,
  metaboanalyst_data = metabo,
  feature = selected$column_name[1]
)

cat("Validation complete\n")
cat("Input CSV files:", length(paths), "\n")
cat("Pixels:", nrow(series$pixel_matrix), "\n")
cat("Shared features:", length(feature_columns(series$pixel_matrix)), "\n")
cat("Selected features:", nrow(selected), "\n")
cat("MetaboAnalyst samples:", nrow(metabo), "\n")
cat("Output directory:", normalizePath(output_dir), "\n")
cat("Plots:\n")
cat(paste(normalizePath(plot_paths), collapse = "\n"), "\n")
