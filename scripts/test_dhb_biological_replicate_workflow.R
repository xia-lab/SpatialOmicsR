# DHB matrix condition: cross-biological-replicate robust feature selection test
#
# IMPORTANT: mPD1 / mPD3 / mPD4 are three DIFFERENT MICE (biological replicates),
# not consecutive slices of the same tissue block. The existing serial-section
# machinery is reused here because it operates on a generic section_column.
# In this script, section_id means mouse/replicate ID, not tissue depth.
#
# Run from the repository root after setting SPATIALOMICS_SMA_DIR.

library(ggplot2)
source("scripts/_bootstrap.R")
load_spatialomics_code()

data_dir <- spatialomics_data_dir(
  "SPATIALOMICS_SMA_DIR", "data_raw/sma_data", "SMA/DHB data"
)
out_dir <- file.path(data_dir, "dhb_replicate_test_outputs")
plot_dir <- file.path(out_dir, "plots")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

replicate_ids <- c("mPD1", "mPD3", "mPD4")

find_dhb_file <- function(root, replicate_id) {
  direct_candidates <- c(
    file.path(root, paste0("ITO2_", replicate_id, "_DHB_220826.csv")),
    file.path(root, paste0("ITO2.", replicate_id, ".DHB.220826.csv"))
  )
  direct_candidates <- direct_candidates[file.exists(direct_candidates)]
  if (length(direct_candidates) > 0) return(direct_candidates[1])

  csv_files <- list.files(root, pattern = "\\.csv$", recursive = TRUE, full.names = TRUE)
  pattern <- paste0("ITO2.*", replicate_id, ".*DHB.*220826\\.csv$")
  matches <- csv_files[grepl(pattern, basename(csv_files), ignore.case = TRUE)]
  if (length(matches) == 0) return(NA_character_)
  matches[1]
}

dhb_files <- vapply(replicate_ids, function(id) find_dhb_file(data_dir, id), character(1))
missing_files <- dhb_files[!file.exists(dhb_files)]
if (length(missing_files) > 0) {
  stop("Missing DHB CSV file(s): ", paste(missing_files, collapse = ", "), call. = FALSE)
}

message("1. Loading peak-picked DHB pixel matrices for 3 biological replicates")
extracted <- load_peakpicked_msi_series(
  paths = dhb_files,
  section_ids = replicate_ids,
  section_order = seq_along(dhb_files)
)

raw_pixel_matrix <- extracted$pixel_matrix
feature_mapping <- extracted$feature_mapping
write.csv(extracted$section_mapping, file.path(out_dir, "replicate_mapping.csv"), row.names = FALSE)

cat("Shared mz features across all 3 mice (DHB):", length(feature_columns(raw_pixel_matrix)), "\n")
cat("Pixels per mouse:\n")
print(table(raw_pixel_matrix$section_id))

raw_fcols <- feature_columns(raw_pixel_matrix)
raw_tic <- rowSums(raw_pixel_matrix[raw_fcols], na.rm = TRUE)
raw_tic_plot_data <- data.frame(
  mouse = raw_pixel_matrix$section_id,
  x = raw_pixel_matrix$x,
  y = raw_pixel_matrix$y,
  log10_tic = log10(raw_tic + 1)
)

p_raw_tic <- ggplot(raw_tic_plot_data, aes(x = x, y = y, fill = log10_tic)) +
  geom_raster() +
  facet_wrap(~mouse, scales = "free") +
  scale_fill_viridis_c(option = "magma") +
  theme_minimal() +
  labs(title = "Raw log10 TIC by mouse (DHB)", x = "x", y = "y", fill = "log10 TIC")

ggsave(file.path(plot_dir, "raw_log10_tic_by_mouse.png"), p_raw_tic, width = 13, height = 5, dpi = 300)

message("2. Cardinal spatial tissue/background segmentation per mouse")
mice <- unique(raw_pixel_matrix$section_id)
replicate_background_results <- lapply(mice, function(mouse_id) {
  mouse_matrix <- raw_pixel_matrix[raw_pixel_matrix$section_id == mouse_id, , drop = FALSE]
  result <- filter_background_cardinal_spatial_kmeans(
    mouse_matrix,
    k = 2,
    r = 1,
    ncomp = 10,
    weights = "adaptive",
    transform = "log10",
    foreground_cleanup = "largest_component",
    seed = 42
  )
  result$background_stats$section_id <- mouse_id
  result
})
names(replicate_background_results) <- mice

tissue_matrix <- do.call(rbind, lapply(replicate_background_results, function(r) r$matrix))
rownames(tissue_matrix) <- NULL
background_stats_all <- do.call(rbind, lapply(replicate_background_results, function(r) r$background_stats))
write.csv(background_stats_all, file.path(out_dir, "cardinal_background_stats_by_mouse.csv"), row.names = FALSE)

mask_all <- do.call(rbind, lapply(mice, function(mouse_id) {
  mouse_matrix <- raw_pixel_matrix[raw_pixel_matrix$section_id == mouse_id, , drop = FALSE]
  r <- replicate_background_results[[mouse_id]]
  data.frame(
    mouse = mouse_id,
    pixel_id = mouse_matrix$pixel_id,
    x = mouse_matrix$x,
    y = mouse_matrix$y,
    tissue_pixel = r$keep
  )
}))
write.csv(mask_all, file.path(out_dir, "cardinal_background_mask_by_mouse.csv"), row.names = FALSE)

p_mask <- ggplot(mask_all, aes(x = x, y = y, color = tissue_pixel)) +
  geom_point(size = 0.3) +
  facet_wrap(~mouse, scales = "free") +
  scale_color_manual(values = c("FALSE" = "#bdbdbd", "TRUE" = "#d73027")) +
  theme_minimal() +
  labs(title = "Cardinal spatialKMeans tissue mask by mouse (DHB)", color = "Tissue")

ggsave(file.path(plot_dir, "cardinal_tissue_mask_by_mouse.png"), p_mask, width = 13, height = 5, dpi = 300)

message("3. Cross-replicate robust feature selection")
feature_selection <- preprocess_select_features(
  tissue_matrix,
  feature_mapping = feature_mapping,
  serial = TRUE,
  section_column = "section_id",
  min_section_fraction = 1,
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
    "No cross-replicate-robust features survived selection. ",
    "Try lowering min_section_fraction, nonzero_min, or cv_top_percent.",
    call. = FALSE
  )
}

write.csv(feature_mapping, file.path(out_dir, "feature_mapping.csv"), row.names = FALSE)
write.csv(selected, file.path(out_dir, "selected_cross_replicate_features.csv"), row.names = FALSE)
write.csv(preprocessed, file.path(out_dir, "preprocessed_matrix.csv"), row.names = FALSE)
write.csv(reduced, file.path(out_dir, "reduced_pixel_matrix.csv"), row.names = FALSE)

cat("Selected-replicate-count per feature:\n")
print(table(selected$selected_section_count))

stats <- feature_stats(preprocessed)
p_cv <- ggplot(stats, aes(x = cv)) +
  geom_histogram(bins = 30, fill = "#3b82f6", color = "white") +
  theme_minimal() +
  labs(title = "Feature CV distribution (all 3 mice combined, DHB)", x = "CV", y = "Number of features")

ggsave(file.path(plot_dir, "feature_cv_distribution.png"), p_cv, width = 8, height = 5, dpi = 300)

message("4. Joint k-means clustering across the shared feature space")
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
  labs(title = "Cluster diagnostics: elbow (DHB)", x = "k", y = "Total within-cluster SS")

ggsave(file.path(plot_dir, "cluster_diagnostics.png"), p_diagnostics, width = 7, height = 5, dpi = 300)

clustered_result <- cluster_pixels(reduced, k = 3, pca_components = 10)
clustered <- clustered_result$matrix
write.csv(clustered, file.path(out_dir, "clustered_matrix.csv"), row.names = FALSE)

p_cluster <- ggplot(clustered, aes(x = x, y = y, fill = factor(cluster))) +
  geom_raster() +
  facet_wrap(~section_id, scales = "free") +
  theme_minimal() +
  labs(title = "Molecular spatial segmentation by mouse (DHB)", x = "x", y = "y", fill = "Cluster")

ggsave(file.path(plot_dir, "kmeans_segmentation_by_mouse.png"), p_cluster, width = 13, height = 5, dpi = 300)

message("5. Region definition from joint clustering and cross-replicate sampling")
clustered$matched_region_label <- paste0("region_", clustered$cluster)

matched <- sample_matched_regions(
  clustered,
  region_column = "matched_region_label",
  section_column = "section_id",
  min_pixels = 30
)

metabo <- make_metaboanalyst_data(matched$sample_matrix, group_column = "matched_region_label")

write.csv(matched$sample_matrix, file.path(out_dir, "matched_sample_matrix.csv"), row.names = FALSE)
write.csv(matched$sample_mapping, file.path(out_dir, "matched_sample_mapping.csv"), row.names = FALSE)
write.csv(metabo, file.path(out_dir, "dhb_metaboanalyst_data.csv"), row.names = FALSE)

p_matched_counts <- ggplot(
  matched$sample_mapping,
  aes(x = matched_region_label, y = n_pixels, fill = section_id)
) +
  geom_col(position = "dodge") +
  theme_minimal() +
  labs(title = "Region pixel counts by mouse (DHB)", x = "Region", y = "Pixels", fill = "Mouse")

ggsave(file.path(plot_dir, "region_counts_by_mouse.png"), p_matched_counts, width = 8, height = 5, dpi = 300)

cat("\n=== DHB cross-replicate workflow test complete ===\n")
cat("Mice (biological replicates):", paste(mice, collapse = ", "), "\n")
cat("Pixels before background filtering (all mice):", nrow(raw_pixel_matrix), "\n")
cat("Pixels after Cardinal tissue filtering (all mice):", nrow(tissue_matrix), "\n")
cat("Cross-replicate-robust features selected:", nrow(selected), "\n")
cat("Cluster sizes (combined):\n")
print(table(clustered$cluster))
cat("Region x mouse sample table:\n")
print(table(matched$sample_matrix$matched_region_label, matched$sample_matrix$section_id))
cat("Output directory:", normalizePath(out_dir), "\n")
cat("Plot directory:", normalizePath(plot_dir), "\n")
