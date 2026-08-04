# Link MetaboAnalystR pathway hits back to MSI spatial images.
#
# R console usage:
# setwd("/Users/ly/Documents/Spatial Omics")
# source("scripts/link_pathway_spatial_maps.R")

source("R/msi_pipeline.R")

if (!exists("data_dir", inherits = TRUE)) {
  data_dir <- "/Users/ly/Desktop/Jeff Xia/rat_brain_data"
}

test_out_dir <- file.path(data_dir, "spatial_test_outputs")
pathway_dir <- file.path(data_dir, "pathway_inputs")
out_dir <- file.path(pathway_dir, "pathway_spatial_maps")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

pixel_matrix_csv <- file.path(test_out_dir, "preprocessed_matrix.csv")
clustered_matrix_csv <- file.path(test_out_dir, "clustered_matrix.csv")
sample_matrix_csv <- file.path(pathway_dir, "all_feature_sample_matrix.csv")
sample_mapping_csv <- file.path(test_out_dir, "metaboanalyst_all_clusters_sample_mapping.csv")
path_a_mset_rds <- file.path(pathway_dir, "path_a_metaboanalyst_fa_name_results", "path_a_fa_name_pathora_mSet.rds")
path_a_feature_csv <- file.path(pathway_dir, "path_a_feature_chebi_stats.csv")
path_b_feature_csv <- file.path(pathway_dir, "path_b_feature_stats.csv")
path_b_mset_rds <- file.path(
  pathway_dir,
  "path_b_metaboanalyst_mummichog_results_update_instrument",
  "path_b_mummichog_mSet.rds"
)

for (path in c(pixel_matrix_csv, clustered_matrix_csv, sample_matrix_csv, sample_mapping_csv, path_a_mset_rds, path_a_feature_csv, path_b_feature_csv, path_b_mset_rds)) {
  if (!file.exists(path)) stop("Missing required file: ", path, call. = FALSE)
}

message("1. Loading spatial matrix and MetaboAnalystR results")
pixel_matrix <- read.csv(pixel_matrix_csv, check.names = FALSE, stringsAsFactors = FALSE)
clustered_matrix <- read.csv(clustered_matrix_csv, check.names = FALSE, stringsAsFactors = FALSE)
path_a_mset <- readRDS(path_a_mset_rds)
path_b_mset <- readRDS(path_b_mset_rds)

message("1b. Re-sampling subregions within each ROI")
cluster_labels <- clustered_matrix[, c("pixel_id", "cluster"), drop = FALSE]
all_feature_clustered <- merge(
  pixel_matrix,
  cluster_labels,
  by = "pixel_id",
  all.x = FALSE,
  sort = FALSE
)
roi_matrix <- define_rois(all_feature_clustered, mode = "cluster")
roi_scope_sampled <- sample_subregions(
  roi_matrix,
  grid_size = 5,
  min_pixels = 30,
  grid_scope = "roi"
)
sample_matrix <- roi_scope_sampled$sample_matrix
sample_mapping <- roi_scope_sampled$sample_mapping
write.csv(sample_matrix, file.path(out_dir, "roi_scope_sample_matrix.csv"), row.names = FALSE)
write.csv(sample_mapping, file.path(out_dir, "roi_scope_sample_mapping.csv"), row.names = FALSE)

message("2. Linking Path A hits back to features")
path_a_feature_meta <- read.csv(path_a_feature_csv, check.names = FALSE, stringsAsFactors = FALSE)
path_a_kegg_map <- data.frame(
  swiss_name_master = c(
    "hexadecanoate",
    "octadecanoate",
    "octadecenoate",
    "octadecadienoate",
    "eicosatetraenoate",
    "docosahexaenoate"
  ),
  compound_id = c("C00249", "C01530", "C00712", "C01595", "C00219", "C06429"),
  compound_name = c(
    "Hexadecanoic acid",
    "Octadecanoic acid",
    "(9Z)-Octadecenoic acid",
    "Linoleate",
    "Arachidonate",
    "(4Z,7Z,10Z,13Z,16Z,19Z)-Docosahexaenoic acid"
  ),
  stringsAsFactors = FALSE
)
path_a_feature_meta <- merge(
  path_a_feature_meta,
  path_a_kegg_map,
  by = "swiss_name_master",
  all.x = TRUE,
  sort = FALSE
)

path_a_name <- "Biosynthesis of unsaturated fatty acids"
path_a_links <- link_metaboanalyst_pathway_features(
  path_a_mset,
  pathway = path_a_name,
  pixel_matrix = pixel_matrix,
  feature_metadata = path_a_feature_meta,
  feature_col = "feature",
  compound_id_col = "compound_id",
  compound_name_col = "compound_name"
)
write.csv(path_a_links, file.path(out_dir, "path_a_unsaturated_fatty_acids_linked_features.csv"), row.names = FALSE)

known_contaminants <- c("mz_255.23298", "mz_283.26430", "mz_554.26209")
path_a_links_clean <- path_a_links[!path_a_links$feature %in% known_contaminants, , drop = FALSE]
write.csv(
  path_a_links_clean,
  file.path(out_dir, "path_a_unsaturated_fatty_acids_known_contaminants_removed_linked_features.csv"),
  row.names = FALSE
)

message("3. Plotting Path A feature ion images and pathway score map")
p_path_a_features <- plot_pathway_feature_images(
  pixel_matrix,
  features = path_a_links$feature,
  pathway_name = paste0("Path A: ", path_a_name),
  transform = "log10",
  ncol = 3
)
ggplot2::ggsave(
  file.path(out_dir, "path_a_unsaturated_fatty_acids_feature_images.png"),
  p_path_a_features,
  width = 9,
  height = 7,
  dpi = 300
)

p_path_a_score <- plot_pathway_score_map(
  roi_matrix,
  features = path_a_links$feature,
  pathway_name = paste0("Path A score: ", path_a_name),
  sample_matrix = sample_matrix,
  sample_mapping = sample_mapping,
  feature_weights = unique(path_a_feature_meta[, c("feature", "signed_log2fc"), drop = FALSE]),
  weight_col = "signed_log2fc",
  feature_col = "feature",
  sample_transform = "identity",
  score_method = "weighted_mean",
  map_level = "roi"
)
ggplot2::ggsave(
  file.path(out_dir, "path_a_unsaturated_fatty_acids_score_map.png"),
  p_path_a_score,
  width = 7,
  height = 7,
  dpi = 300
)

message("3b. Plotting Path A maps after removing known contaminant features")
p_path_a_clean_features <- plot_pathway_feature_images(
  pixel_matrix,
  features = path_a_links_clean$feature,
  pathway_name = paste0("Path A clean: ", path_a_name),
  transform = "log10",
  ncol = 2
)
ggplot2::ggsave(
  file.path(out_dir, "path_a_unsaturated_fatty_acids_known_contaminants_removed_feature_images.png"),
  p_path_a_clean_features,
  width = 8,
  height = 7,
  dpi = 300
)

p_path_a_clean_score <- plot_pathway_score_map(
  roi_matrix,
  features = path_a_links_clean$feature,
  pathway_name = paste0("Path A clean score: ", path_a_name),
  sample_matrix = sample_matrix,
  sample_mapping = sample_mapping,
  feature_weights = unique(path_a_feature_meta[
    path_a_feature_meta$feature %in% path_a_links_clean$feature,
    c("feature", "signed_log2fc"),
    drop = FALSE
  ]),
  weight_col = "signed_log2fc",
  feature_col = "feature",
  sample_transform = "identity",
  score_method = "weighted_mean",
  map_level = "roi"
)
ggplot2::ggsave(
  file.path(out_dir, "path_a_unsaturated_fatty_acids_known_contaminants_removed_score_map.png"),
  p_path_a_clean_score,
  width = 7,
  height = 7,
  dpi = 300
)

message("4. Linking Path B top mummichog pathway back to features")
path_b_feature_meta <- read.csv(path_b_feature_csv, check.names = FALSE, stringsAsFactors = FALSE)
path_b_name <- path_b_mset$path.nms[[1]]
path_b_links <- link_metaboanalyst_pathway_features(
  path_b_mset,
  pathway = path_b_name,
  pixel_matrix = pixel_matrix
)
write.csv(path_b_links, file.path(out_dir, "path_b_top_mummichog_linked_features.csv"), row.names = FALSE)

message("5. Plotting Path B feature ion images and pathway score map")
p_path_b_features <- plot_pathway_feature_images(
  pixel_matrix,
  features = path_b_links$feature,
  pathway_name = paste0("Path B: ", path_b_name),
  transform = "log10",
  ncol = 3
)
ggplot2::ggsave(
  file.path(out_dir, "path_b_top_mummichog_feature_images.png"),
  p_path_b_features,
  width = 9,
  height = 7,
  dpi = 300
)

p_path_b_score <- plot_pathway_score_map(
  roi_matrix,
  features = path_b_links$feature,
  pathway_name = paste0("Path B score: ", path_b_name),
  sample_matrix = sample_matrix,
  sample_mapping = sample_mapping,
  feature_weights = path_b_feature_meta[path_b_feature_meta$feature %in% path_b_links$feature, c("feature", "signed_log2fc"), drop = FALSE],
  weight_col = "signed_log2fc",
  feature_col = "feature",
  sample_transform = "identity",
  score_method = "weighted_mean",
  map_level = "roi"
)
ggplot2::ggsave(
  file.path(out_dir, "path_b_top_mummichog_score_map.png"),
  p_path_b_score,
  width = 7,
  height = 7,
  dpi = 300
)

cat("\n=== Pathway spatial linking complete ===\n")
cat("Path A linked features:", nrow(path_a_links), "\n")
cat("Path A linked features after known contaminant removal:", nrow(path_a_links_clean), "\n")
cat("Path B pathway:", path_b_name, "\n")
cat("Path B linked features:", nrow(path_b_links), "\n")
cat("Output directory:", normalizePath(out_dir), "\n")
