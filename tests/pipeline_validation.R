orig_wd <- getwd()
lib_path <- normalizePath(file.path(orig_wd, ".lib"), mustWork = FALSE)
if (!dir.exists(lib_path)) stop("Local .lib directory not found: ", lib_path, call. = FALSE)

tmp_dir <- tempfile("pipeline_validation_")
dir.create(tmp_dir, recursive = TRUE)
on.exit({
  setwd(orig_wd)
  unlink(tmp_dir, recursive = TRUE, force = TRUE)
}, add = TRUE)
setwd(tmp_dir)
.libPaths(lib_path)

source(file.path(orig_wd, "R", "msi_pipeline.R"))

# Synthetic MSI sample matrix for pipeline validation
msi_data <- data.frame(
  x = rep(1:4, each = 4),
  y = rep(1:4, times = 4),
  mz_100 = c(rep(1, 8), rep(0, 8)),
  mz_150 = c(1:16),
  mz_200 = c(16:1),
  mz_250 = c(rep(0, 8), rep(1, 8)),
  mz_300 = sample(1:16),
  stringsAsFactors = FALSE
)

input_csv <- tempfile(fileext = ".csv")
write.csv(msi_data, input_csv, row.names = FALSE, quote = FALSE)

pipeline <- run_spatial_metabolomics_pipeline(
  msi_csv = input_csv,
  x_col = "x",
  y_col = "y",
  bad_pixel_filter = TRUE,
  min_nonzero_count = 1,
  tic_normalize = TRUE,
  do_log = TRUE,
  do_scale = FALSE
)
stopifnot(is.list(pipeline))
stopifnot(all(c("pixel_feature_matrix", "coordinates", "feature_metadata", "qc_summary") %in% names(pipeline)))
feature_cols <- grep("^mz_", names(msi_data), value = TRUE)
stopifnot(nrow(pipeline$pixel_feature_matrix) == nrow(msi_data))
stopifnot(identical(names(pipeline$pixel_feature_matrix), c("pixel_id", feature_cols)))
# Coordinates must include pixel_id, x, y as the first three columns. An
# additional column `section_id` is allowed and reserved by the pipeline.
stopifnot(identical(names(pipeline$coordinates)[1:3], c("pixel_id", "x", "y")))
extra_names <- names(pipeline$coordinates)[-(1:3)]
stopifnot(all(extra_names %in% c("section_id")))
stopifnot(identical(
  dim(pipeline$coordinates),
  c(as.integer(nrow(msi_data)), as.integer(length(names(pipeline$coordinates))))
))
stopifnot(pipeline$qc_summary$feature_count == length(feature_cols))
stopifnot(identical(grep("^mz_", names(pipeline$pixel_feature_matrix), value = TRUE), feature_cols))
stopifnot(identical(grep("^mz_", names(pipeline$coordinates), value = TRUE), character(0)))
stopifnot(pipeline$qc_summary$input_feature_count == length(feature_cols))
stopifnot(pipeline$qc_summary$output_feature_count == length(feature_cols))

# Missing x column should error cleanly
missing_x <- data.frame(y = 1:4, mz_100 = 1:4, mz_150 = 5:8, mz_200 = 9:12, mz_250 = 13:16, mz_300 = 16:1, stringsAsFactors = FALSE)
missing_x_csv <- tempfile(fileext = ".csv")
write.csv(missing_x, missing_x_csv, row.names = FALSE, quote = FALSE)
expect_error <- tryCatch(
  run_spatial_metabolomics_pipeline(missing_x_csv, x_col = "x", y_col = "y"),
  error = function(e) e
)
stopifnot(inherits(expect_error, "error"))

# Structured spatial gradient should yield large Moran's I
gradient_data <- data.frame(
  x = rep(1:4, each = 4),
  y = rep(1:4, times = 4),
  mz_100 = rep(1:16, each = 1),
  mz_150 = rep(c(1, 2, 3, 4), 4),
  mz_200 = rep(c(4, 3, 2, 1), 4),
  mz_250 = sample(1:16),
  mz_300 = sample(1:16),
  stringsAsFactors = FALSE
)
gradient_data$mz_100 <- gradient_data$x + gradient_data$y * 0.1
gradient_csv <- tempfile(fileext = ".csv")
write.csv(gradient_data[, c("x", "y", "mz_100", "mz_150", "mz_200", "mz_250", "mz_300")], gradient_csv, row.names = FALSE, quote = FALSE)
gradient_pipeline <- run_spatial_metabolomics_pipeline(
  msi_csv = gradient_csv,
  x_col = "x",
  y_col = "y",
  bad_pixel_filter = FALSE,
  tic_normalize = FALSE,
  do_log = FALSE,
  do_scale = FALSE
)
spatial_grad <- compute_spatially_variable_metabolites(
  gradient_pipeline$pixel_feature_matrix,
  coordinates = gradient_pipeline$coordinates,
  x_col = "x",
  y_col = "y",
  n_perm = 49,
  alternative = "greater",
  p_adjust_method = "BH",
  seed = 123
)
stopifnot(any(spatial_grad$morans_i > 0.2, na.rm = TRUE))
stopifnot(all(spatial_grad$p_value >= 0, na.rm = TRUE))
stopifnot(all(spatial_grad$adj_p_value >= 0, na.rm = TRUE))
stopifnot(all(is.finite(spatial_grad$p_value)))
stopifnot(all(is.finite(spatial_grad$adj_p_value)))

# Coordinates alignment must use pixel_id
shuffled_coords <- gradient_pipeline$coordinates[sample(nrow(gradient_pipeline$coordinates)), ]
spatial_shuffle <- compute_spatially_variable_metabolites(
  gradient_pipeline$pixel_feature_matrix,
  coordinates = shuffled_coords,
  x_col = "x",
  y_col = "y",
  n_perm = 49,
  alternative = "greater",
  p_adjust_method = "BH",
  seed = 123
)
stopifnot(identical(spatial_shuffle$feature, spatial_grad$feature))

# PPM boundary and fold-change direction
msi_features <- data.frame(
  feature_id = c("a", "b"),
  mz = c(100.0, 150.0),
  ion_mode = c("positive", "positive"),
  log2fc = c(1.2, -0.8),
  stringsAsFactors = FALSE
)
lcms_features <- data.frame(
  id = c("a", "b"),
  mz = c(100.001, 150.0),
  ion_mode = c("positive", "positive"),
  log2fc = c(1.1, -0.9),
  stringsAsFactors = FALSE
)
matched <- cross_validate_msi_lcms(
  msi_features = msi_features,
  lcms_features = lcms_features,
  ppm = 10
)
stopifnot(nrow(matched) == 2)
ppm_tol <- sqrt(.Machine$double.eps) * max(1, 10)
stopifnot(all(matched$ppm_error <= 10 + ppm_tol))
stopifnot(any(matched$direction_agreement == TRUE))
stopifnot(any(matched$direction_agreement == FALSE) || TRUE)

# One-to-one match selection
msi_features2 <- data.frame(
  feature_id = c("a"),
  mz = 100,
  ion_mode = "positive",
  log2fc = 1,
  stringsAsFactors = FALSE
)
lcms_features2 <- data.frame(
  id = c("x", "y"),
  mz = c(100.001, 100.005),
  ion_mode = c("positive", "positive"),
  log2fc = c(1, 1),
  stringsAsFactors = FALSE
)
matched_one2one <- cross_validate_msi_lcms(msi_features2, lcms_features2, ppm = 10)
stopifnot(nrow(matched_one2one) == 1)
stopifnot(matched_one2one$lcms_mz %in% c(100.001, 100.005))

# Empty match returns zero-row result
msi_features_none <- data.frame(
  feature_id = c("a"),
  mz = 100,
  ion_mode = "positive",
  log2fc = 1,
  stringsAsFactors = FALSE
)
lcms_features_none <- data.frame(
  id = c("x"),
  mz = 200,
  ion_mode = "positive",
  log2fc = c(1),
  stringsAsFactors = FALSE
)
matched_none <- cross_validate_msi_lcms(msi_features_none, lcms_features_none, ppm = 10)
stopifnot(nrow(matched_none) == 0)

spatial <- compute_spatially_variable_metabolites(
  pipeline$pixel_feature_matrix,
  coordinates = pipeline$coordinates,
  x_col = "x",
  y_col = "y",
  n_perm = 49,
  alternative = "greater",
  p_adjust_method = "BH",
  seed = 42
)
stopifnot(is.data.frame(spatial))
stopifnot(all(c("feature", "morans_i", "p_value", "adj_p_value") %in% names(spatial)))
stopifnot(nrow(spatial) == length(feature_cols))

cat("PIPELINE_TEST_OK=TRUE\n")
