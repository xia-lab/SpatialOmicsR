testthat::test_that("regression: test-hmrf_clustering", {
if (nzchar(Sys.getenv("_R_CHECK_PACKAGE_NAME_"))) {
  library(SpatialOmicsMSI)
} else {
  source(file.path(getwd(), "R", "msi_pipeline.R"))
}

pixel <- data.frame(
  pixel_id = 1:16,
  x = rep(1:4, each = 4),
  y = rep(1:4, times = 4),
  mz_100 = c(rep(0, 8), rep(10, 8)),
  mz_200 = c(rep(10, 8), rep(0, 8)),
  check.names = FALSE
)
rook <- build_spatial_neighbors(pixel$x, pixel$y, method = "rook")

fit_a <- cluster_pixels_hmrf(
  pixel, k = 2, beta = 1, neighbors = rook,
  update_order = "fixed", seed = 9
)
fit_b <- cluster_pixels_hmrf(
  pixel, k = 2, beta = 1, neighbors = rook,
  update_order = "fixed", seed = 9
)
stopifnot(identical(fit_a$labels, fit_b$labels))
stopifnot(identical(fit_a$full_bayesian_hmrf, FALSE))
stopifnot(all(fit_a$fit$size > 0))
stopifnot(all(diff(fit_a$iteration_log$total_energy) <= 1e-8))

# With beta zero and identical preprocessing, the converged k-means
# initialization is already an ICM fixed point.
nonspatial <- cluster_pixels_hmrf(
  pixel, k = 2, beta = 0, neighbors = rook, scale_features = FALSE,
  data_term_scale = "raw", update_order = "fixed", seed = 9
)
ordinary <- cluster_pixels(pixel, k = 2, seed = 9)
stopifnot(identical(nonspatial$labels, ordinary$matrix$cluster))

diagnostics <- cluster_diagnostics_hmrf(
  pixel, k = 2, beta_grid = c(0, 0.5, 2), neighbor_method = "rook",
  update_order = "fixed", seed = 9
)
stopifnot(nrow(diagnostics) == 3L)
stopifnot(all(diagnostics$adjacent_pair_agreement >= 0 & diagnostics$adjacent_pair_agreement <= 1))
stopifnot(nzchar(attr(diagnostics, "interpretation")))

# A deliberately discordant interior pixel is absorbed at sufficiently high
# beta in this simple two-domain example, reducing boundary-edge fraction.
noisy <- expand.grid(x = 1:7, y = 1:7)
noisy$pixel_id <- seq_len(nrow(noisy))
noisy$mz_100 <- ifelse(noisy$x <= 3, 0, 10)
noisy$mz_200 <- 10 - noisy$mz_100
noisy[noisy$x == 2 & noisy$y == 4, c("mz_100", "mz_200")] <- c(10, 0)
smoothness <- cluster_diagnostics_hmrf(
  noisy[, c("pixel_id", "x", "y", "mz_100", "mz_200")],
  k = 2, beta_grid = c(0, 4), neighbor_method = "rook",
  update_order = "fixed", seed = 9
)
stopifnot(smoothness$boundary_edge_fraction[2] < smoothness$boundary_edge_fraction[1])

  testthat::succeed()
})
