testthat::test_that("regression: test-spatial_lmm", {
if (nzchar(Sys.getenv("_R_CHECK_PACKAGE_NAME_"))) {
  library(SpatialOmicsMSI)
} else {
  source(file.path(getwd(), "R", "msi_pipeline.R"))
  source(file.path(getwd(), "R", "spatial_lmm.R"))
}

spatial_samples <- expand.grid(
  subject_id = paste0("subject_", 1:6),
  section_id = c("section_1", "section_2"),
  x_center = 1:4,
  y_center = 1:3,
  stringsAsFactors = FALSE
)
spatial_samples$roi_id <- ifelse(spatial_samples$x_center <= 2, "region_a", "region_b")
spatial_samples$sample_id <- paste(
  spatial_samples$subject_id, spatial_samples$section_id,
  spatial_samples$x_center, spatial_samples$y_center, sep = "::"
)
subject_effect <- stats::setNames(seq(-0.25, 0.25, length.out = 6), paste0("subject_", 1:6))
set.seed(42)
spatial_samples$mz_100 <-
  unname(subject_effect[spatial_samples$subject_id]) +
  ifelse(spatial_samples$roi_id == "region_b", 2, 0) +
  0.05 * spatial_samples$x_center + stats::rnorm(nrow(spatial_samples), sd = 0.15)
# Give the second feature an explicit, well-identified exponential spatial
# covariance within every subject-section field. A nearly pure-noise fixture
# makes the range parameter weakly identifiable and lets nlme land on different
# optimizer boundaries across BLAS/compiler builds.
spatial_samples$mz_200 <- NA_real_
field_indices <- split(
  seq_len(nrow(spatial_samples)),
  interaction(spatial_samples$subject_id, spatial_samples$section_id, drop = TRUE)
)
for (index in field_indices) {
  physical_coordinates <- cbind(
    50 * spatial_samples$x_center[index],
    25 * spatial_samples$y_center[index]
  )
  distance_matrix <- as.matrix(stats::dist(physical_coordinates))
  covariance <- 0.2^2 * (
    0.9 * exp(-distance_matrix / 75) +
      0.1 * diag(length(index))
  )
  spatial_error <- drop(t(chol(covariance)) %*% stats::rnorm(length(index)))
  spatial_samples$mz_200[index] <-
    unname(subject_effect[spatial_samples$subject_id[index]]) + spatial_error
}

fit <- differential_region_analysis_spatial_lmm(
  spatial_samples,
  group_column = "roi_id",
  subject_column = "subject_id",
  section_column = "section_id",
  x_col = "x_center",
  y_col = "y_center",
  x_resolution = 50,
  y_resolution = 25,
  distance_unit = "um",
  correlation_structure = "exponential",
  reference_group = "region_a",
  min_subjects_per_group = 5,
  min_unique_coordinates_per_field = 4
)
stopifnot(is.list(fit))
stopifnot(all(c("features", "contrasts", "field_qc", "settings") %in% names(fit)))
stopifnot(all(fit$features$converged))
stopifnot(all(fit$features$range_parameter > 0))
stopifnot(all(fit$features$nugget >= 0 & fit$features$nugget <= 1))
stopifnot(all(is.finite(fit$features$AIC_spatial)))
stopifnot(!"likelihood_ratio_p" %in% names(fit$features))
effect <- fit$contrasts$estimate[
  fit$contrasts$feature == "mz_100" &
    fit$contrasts$group_a == "region_a" & fit$contrasts$group_b == "region_b"
]
stopifnot(length(effect) == 1L, effect > 1)
stopifnot(all(fit$field_qc$n_duplicate_coordinates == 0L))
stopifnot(grepl("Gaussian spatial LMM", fit$interpretation, fixed = TRUE))

duplicate_samples <- rbind(spatial_samples, spatial_samples[1, ])
duplicate_samples$sample_id[nrow(duplicate_samples)] <- "duplicate"
duplicate_failed <- tryCatch({
  differential_region_analysis_spatial_lmm(
    duplicate_samples,
    subject_column = "subject_id", section_column = "section_id",
    x_col = "x_center", y_col = "y_center",
    x_resolution = 1, y_resolution = 1, distance_unit = "grid",
    min_subjects_per_group = 5
  )
  FALSE
}, error = function(e) TRUE)
stopifnot(duplicate_failed)

alias_warning <- NULL
alias_fit <- withCallingHandlers(
  differential_region_analysis_glmm(
    spatial_samples,
    subject_column = "subject_id", section_column = "section_id",
    x_col = "x_center", y_col = "y_center",
    x_resolution = 50, y_resolution = 25, distance_unit = "um",
    min_subjects_per_group = 5
  ),
  warning = function(w) {
    alias_warning <<- conditionMessage(w)
    invokeRestart("muffleWarning")
  }
)
stopifnot(grepl("compatibility alias", alias_warning, fixed = TRUE))
stopifnot(is.list(alias_fit))

  testthat::succeed()
})
