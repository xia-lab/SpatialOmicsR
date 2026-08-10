testthat::test_that("regression: test-variable_mz_adapter", {
checking_installed_package <- nzchar(Sys.getenv("_R_CHECK_PACKAGE_NAME_"))
if (dir.exists(".lib")) .libPaths(c(normalizePath(".lib"), .libPaths()))
if (checking_installed_package) library(SpatialOmicsMSI) else {
  source(file.path(normalizePath(getwd()), "R", "msi_pipeline.R"))
  source(file.path(normalizePath(getwd()), "R", "real_data_adapters.R"))
}

expect_error <- function(x, pattern) {
  e <- tryCatch(x, error = identity)
  stopifnot(inherits(e, "error"), grepl(pattern, conditionMessage(e), fixed = TRUE))
}

d <- tempfile("imz_pair_"); dir.create(d)
imz <- file.path(d, "a.imzML"); ibd <- file.path(d, "a.ibd")
writeLines("<mzML><fileDescription/></mzML>", imz)
writeBin(as.raw(seq_len(16)), ibd)
pair <- validate_imzml_ibd_pair(imz, ibd)
stopifnot(identical(pair$imzml_path, normalizePath(imz)))
expect_error(validate_imzml_ibd_pair(imz, file.path(d, "missing.ibd")), "companion ibd_path")

# Slow local real-data alignment is explicitly opt-in and never runs in routine CI.
run_real_data_tests <- identical(
  tolower(Sys.getenv("SPATIALOMICS_RUN_REAL_DATA_TESTS", unset = "false")), "true"
)
if (!checking_installed_package && run_real_data_tests && requireNamespace("Cardinal", quietly = TRUE)) {
  # Real-data contract checks are below; unit behavior is exercised through public inputs.
  root <- file.path(getwd(), "data_raw", "mouse_brain_he_msi", "metaspace_brain01")
  path <- file.path(root, "Brain01_Bregma-1-46_centroid.imzML")
  if (file.exists(path)) {
    x <- load_variable_mz_msi_features(path,
      sample_id = "Brain01", section_id = "Brain01", ion_mode = "positive",
      ion_mode_source = "METASPACE_dataset_metadata",
      processing = "processed_peak_lists", alignment_ppm = 10,
      min_detection_fraction = 0)
    stopifnot(nrow(x$pixel_feature_matrix) == 17809L,
              nrow(x$coordinates) == 17809L,
              !anyDuplicated(x$feature_metadata$column_name),
              identical(x$parameters$peak_pick_method, "not_applied"),
              all(x$coordinates$x == x$pixel_feature_matrix$x))
  }
}
cat("VARIABLE_MZ_ADAPTER_TEST_OK=TRUE\n")

  testthat::succeed()
})
