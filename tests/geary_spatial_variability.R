if (nzchar(Sys.getenv("_R_CHECK_PACKAGE_NAME_"))) {
  library(SpatialOmicsMSI)
} else {
  source(file.path(getwd(), "R", "msi_pipeline.R"))
}

chain <- build_spatial_neighbors(1:3, rep(1, 3), method = "rook")
hand <- compute_gearys_c_grid(c(1, 2, 3), neighbors = chain)
stopifnot(isTRUE(all.equal(hand$C, 0.5)))
stopifnot(identical(hand$expected_C, 1))

alternating <- compute_gearys_c_grid(c(1, 3, 1), neighbors = chain)
stopifnot(isTRUE(all.equal(alternating$C, 1.5)))

# An isolated extreme value cannot alter the connected-subgraph statistic.
with_isolate <- build_spatial_neighbors(c(1:3, 20), c(1, 1, 1, 20), method = "rook")
isolated <- compute_gearys_c_grid(c(1, 2, 3, 1e9), neighbors = with_isolate)
stopifnot(isTRUE(all.equal(isolated$C, hand$C)))
stopifnot(isolated$n_isolated == 1L)

seed_a <- compute_gearys_c_grid(c(1, 2, 4, 8), x = 1:4, y = rep(1, 4), n_perm = 49, seed = 22)
seed_b <- compute_gearys_c_grid(c(1, 2, 4, 8), x = 1:4, y = rep(1, 4), n_perm = 49, seed = 22)
stopifnot(identical(seed_a$p_value, seed_b$p_value))
stopifnot(seed_a$p_value >= 1 / 50)

pixel <- data.frame(
  pixel_id = 1:9,
  x = rep(1:3, each = 3), y = rep(1:3, times = 3),
  mz_100 = rep(1:3, each = 3),
  mz_200 = c(1, 8, 2, 9, 3, 7, 4, 6, 5),
  check.names = FALSE
)
geary <- compute_spatially_variable_metabolites_geary(
  pixel, n_perm = 49, alternative = "less", seed = 8
)
moran <- compute_spatially_variable_metabolites(
  pixel, n_perm = 49, alternative = "greater", seed = 8
)
stopifnot(all(c("feature", "gearys_c", "p_value", "adj_p_value") %in% names(geary)))
comparison <- compare_moran_geary(moran, geary)
stopifnot(nrow(comparison) == 2L)
stopifnot(nzchar(attr(comparison, "interpretation")))
