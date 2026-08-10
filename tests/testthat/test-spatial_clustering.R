testthat::test_that("regression: test-spatial_clustering", {
if (nzchar(Sys.getenv("_R_CHECK_PACKAGE_NAME_"))) {
  library(SpatialOmicsMSI)
} else {
  source(file.path(getwd(), "R", "msi_pipeline.R"))
}

pixel <- data.frame(
  pixel_id = 1:9,
  x = rep(1:3, each = 3),
  y = rep(1:3, times = 3),
  mz_100 = c(1, 1, 1, 5, 5, 5, 9, 9, 9),
  mz_200 = c(9, 9, 9, 5, 5, 5, 1, 1, 1),
  check.names = FALSE
)

rook <- build_spatial_neighbors(pixel$x, pixel$y, method = "rook")
h0 <- compute_neighborhood_average(pixel[c("mz_100", "mz_200")], rook)
stopifnot(identical(dim(h0), c(9L, 2L)))
stopifnot(isTRUE(all.equal(unname(h0[1, ]), c(3, 7))))
stopifnot(!anyNA(h0))

fit_a <- cluster_pixels_spatial(pixel, k = 3, lambda = 0.5, neighbors = rook, seed = 17)
fit_b <- cluster_pixels_spatial(pixel, k = 3, lambda = 0.5, neighbors = rook, seed = 17)
stopifnot(identical(fit_a$matrix$cluster, fit_b$matrix$cluster))
stopifnot(identical(fit_a$full_banksy, FALSE))
stopifnot(ncol(fit_a$augmented_features) == 4L)

isolated_pixel <- rbind(pixel[1:2, ], transform(pixel[3, ], pixel_id = 10, x = 20, y = 20))
isolated_graph <- build_spatial_neighbors(isolated_pixel$x, isolated_pixel$y, method = "rook")
isolated_h0 <- compute_neighborhood_average(
  isolated_pixel[c("mz_100", "mz_200")], isolated_graph, isolate_action = "self"
)
stopifnot(isTRUE(all.equal(unname(isolated_h0[3, ]), unname(as.numeric(isolated_pixel[3, c("mz_100", "mz_200")])))))

diagnostics <- cluster_diagnostics_spatial(
  pixel, k = 3, lambda_grid = c(0, 0.5, 1), neighbor_method = "rook", seed = 17
)
stopifnot(nrow(diagnostics) == 3L)
stopifnot(all(diagnostics$adjacent_pair_agreement >= 0 & diagnostics$adjacent_pair_agreement <= 1))
stopifnot(nzchar(attr(diagnostics, "interpretation")))

  testthat::succeed()
})
