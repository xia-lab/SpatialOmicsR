# Validate spatial-niche functionality on the three real-data collections kept
# under data_raw. Each analysis branch is saved as one self-contained RDS file.

source("scripts/_bootstrap.R")
load_spatialomics_code()

output_dir <- spatialomics_path("outputs", "niche_real_data_validation")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
set.seed(20260810)

make_exploratory_domains <- function(pixel_matrix, k_domains = 4L, max_features = 200L) {
  fcols <- feature_columns(pixel_matrix)
  values <- as.matrix(pixel_matrix[fcols])
  storage.mode(values) <- "double"
  values[!is.finite(values) | values < 0] <- 0
  totals <- rowSums(values)
  normalized <- values / pmax(totals, .Machine$double.eps)
  normalized <- log1p(normalized * 1e6)
  feature_variance <- apply(normalized, 2L, stats::var)
  keep <- order(feature_variance, decreasing = TRUE)[
    seq_len(min(max_features, sum(is.finite(feature_variance) & feature_variance > 0)))
  ]
  if (length(keep) < 2L) stop("Too few variable MSI features for exploratory domains.", call. = FALSE)
  reduced <- stats::prcomp(normalized[, keep, drop = FALSE], center = TRUE, scale. = TRUE, rank. = 10)$x
  fit <- stats::kmeans(reduced, centers = k_domains, nstart = 25, iter.max = 100)
  out <- pixel_matrix[c("pixel_id", "sample_id", "section_id", "x", "y")]
  out$domain <- paste0("domain_", fit$cluster)
  attr(out, "domain_method") <- "joint PCA/k-means on TIC-normalized log1p intensities"
  out
}

infer_grid_radius <- function(data) {
  positive_steps <- unlist(lapply(split(seq_len(nrow(data)), data$section_id), function(index) {
    c(diff(sort(unique(data$x[index]))), diff(sort(unique(data$y[index]))))
  }), use.names = FALSE)
  step <- stats::median(positive_steps[is.finite(positive_steps) & positive_steps > 0])
  if (!is.finite(step)) stop("Could not infer a positive coordinate increment.", call. = FALSE)
  2.01 * step
}

load_full_tissue <- function() {
  root <- spatialomics_path("data_raw", "full_tissue_mouse_brain", "OMIX016317")
  loaded <- load_variable_mz_msi_features(
    file.path(root, "OMIX016317-02.imzML"),
    ibd_path = file.path(root, "OMIX016317-01.ibd"),
    sample_id = "OMIX016317", section_id = "OMIX016317",
    ion_mode = "negative", ion_mode_source = "dataset_metadata",
    processing = "profile_diff", alignment_ppm = 10, peak_pick_snr = 6,
    min_detection_fraction = 0.10
  )
  tissue <- build_msi_tissue_mask(
    loaded$coordinates, loaded$pixel_qc$raw_tic, loaded$pixel_qc$raw_peak_count,
    method = "kmeans_log_tic_peak_count", seed = 20260810
  )$mask$tissue
  make_exploratory_domains(loaded$pixel_feature_matrix[tissue, , drop = FALSE])
}

load_mouse_brain <- function() {
  root <- spatialomics_path("data_raw", "mouse_brain_he_msi", "metaspace_brain01")
  loaded <- load_variable_mz_msi_features(
    file.path(root, "Brain01_Bregma-1-46_centroid.imzML"),
    sample_id = "Brain01", section_id = "Brain01",
    ion_mode = "positive", ion_mode_source = "METASPACE dataset metadata",
    processing = "processed_peak_lists", alignment_ppm = 10,
    min_detection_fraction = 0.10
  )
  tissue <- build_msi_tissue_mask(
    loaded$coordinates, loaded$pixel_qc$raw_tic, loaded$pixel_qc$raw_peak_count,
    method = "kmeans_log_tic_peak_count", seed = 20260810
  )$mask$tissue
  make_exploratory_domains(loaded$pixel_feature_matrix[tissue, , drop = FALSE])
}

load_msiflow <- function() {
  path <- spatialomics_path("data_raw", "msiflow", "ly6g_heterogeneity_signatures", "umap_data.csv")
  domains <- utils::read.csv(path, stringsAsFactors = FALSE)
  domains <- domains[as.character(domains$label) != "-1", , drop = FALSE]
  data.frame(
    pixel_id = paste(domains$sample, domains$x, domains$y, sep = "::"),
    sample_id = domains$sample,
    section_id = domains$sample,
    x = domains$x,
    y = domains$y,
    domain = paste0("domain_", domains$label),
    stringsAsFactors = FALSE
  )
}

run_branches <- function(dataset_name, domains) {
  dataset_dir <- file.path(output_dir, dataset_name)
  dir.create(dataset_dir, recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(domains, file.path(dataset_dir, "joint_domain_input.csv"), row.names = FALSE)
  radius <- infer_grid_radius(domains)
  branches <- list(
    radius_hellinger = list(neighborhood = "radius", transform = "hellinger", data = domains),
    radius_none = list(neighborhood = "radius", transform = "none", data = domains),
    radius_clr = list(neighborhood = "radius", transform = "clr", data = domains)
  )
  # The compatibility kNN implementation is quadratic. Exercise it on a fixed
  # real-data subset while radius branches use all retained tissue pixels.
  knn_index <- unlist(lapply(split(seq_len(nrow(domains)), domains$section_id), function(index) {
    if (length(index) <= 2000L) index else sort(sample(index, 2000L))
  }), use.names = FALSE)
  branches$knn_hellinger <- list(
    neighborhood = "knn", transform = "hellinger", data = domains[knn_index, , drop = FALSE]
  )
  summaries <- lapply(names(branches), function(branch_name) {
    specification <- branches[[branch_name]]
    result <- detect_spatial_niches(
      specification$data,
      domain_column = "domain", x_col = "x", y_col = "y",
      subject_column = "sample_id", section_column = "section_id",
      neighborhood = specification$neighborhood,
      radius = if (specification$neighborhood == "radius") radius else NULL,
      k = 10, include_self = TRUE, min_neighbors = 5,
      transform = specification$transform, clr_pseudocount = 1e-6,
      k_niches = 4, domain_alignment = "joint", seed = 20260810
    )
    result$validation <- list(
      dataset = dataset_name, branch = branch_name,
      input_pixels = nrow(specification$data), full_dataset_pixels = nrow(domains),
      radius_coordinate_units = if (specification$neighborhood == "radius") radius else NA_real_
    )
    saveRDS(result, file.path(dataset_dir, paste0(branch_name, ".rds")), compress = "xz")
    data.frame(
      dataset = dataset_name, branch = branch_name,
      neighborhood = specification$neighborhood, transform = specification$transform,
      n_input = nrow(specification$data), n_eligible = sum(result$matrix$eligible),
      n_excluded = sum(!result$matrix$eligible),
      n_niches = length(unique(stats::na.omit(result$matrix$niche_id))),
      mean_ambiguity_ratio = mean(result$matrix$ambiguity_ratio, na.rm = TRUE),
      result_file = file.path(dataset_name, paste0(branch_name, ".rds")),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, summaries)
}

loaders <- list(
  full_tissue_mouse_brain = load_full_tissue,
  mouse_brain = load_mouse_brain,
  msiflow = load_msiflow
)
all_summaries <- lapply(names(loaders), function(dataset_name) {
  message("Loading and testing ", dataset_name)
  domains <- loaders[[dataset_name]]()
  run_branches(dataset_name, domains)
})
manifest <- do.call(rbind, all_summaries)
utils::write.csv(manifest, file.path(output_dir, "result_manifest.csv"), row.names = FALSE)
saveRDS(manifest, file.path(output_dir, "result_manifest.rds"), compress = "xz")
print(manifest)
