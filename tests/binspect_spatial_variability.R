if (nzchar(Sys.getenv("_R_CHECK_PACKAGE_NAME_"))) {
  library(SpatialOmicsMSI)
} else {
  source(file.path(getwd(), "R", "msi_pipeline.R"))
}

grid <- expand.grid(x = 1:6, y = 1:6)
grid$pixel_id <- seq_len(nrow(grid))
grid$mz_block <- as.numeric(grid$x <= 3)
grid$mz_mixed <- as.numeric((grid$x + grid$y) %% 2 == 0)
grid$mz_constant <- 1
grid <- grid[, c("pixel_id", "x", "y", "mz_block", "mz_mixed", "mz_constant")]
neighbors <- build_spatial_neighbors(grid$x, grid$y, method = "rook")

block_a <- compute_binspect_feature(
  grid$mz_block, neighbors, bin_method = "rank", percentage_rank = 50,
  inference = "permutation", n_perm = 99, seed = 14
)
block_b <- compute_binspect_feature(
  grid$mz_block, neighbors, bin_method = "rank", percentage_rank = 50,
  inference = "permutation", n_perm = 99, seed = 14
)
stopifnot(identical(block_a$p_value, block_b$p_value))
stopifnot(block_a$p_value >= 0.01)
stopifnot(block_a$odds_ratio > 1)

fisher <- compute_binspect_feature(
  grid$mz_block, neighbors, bin_method = "rank", percentage_rank = 50,
  inference = "fisher_approx", seed = 14
)
stopifnot(identical(fisher$inference, "fisher_approx"))
stopifnot(isTRUE(all.equal(fisher$contingency_table, t(fisher$contingency_table))))

batch <- compute_spatially_variable_metabolites_binspect(
  grid, bin_method = "rank", percentage_rank = 50,
  inference = "permutation", n_perm = 49, seed = 14
)
stopifnot(batch$status[match("mz_constant", batch$feature)] == "not_tested")
stopifnot(nzchar(batch$error_message[match("mz_constant", batch$feature)]))

moran <- compute_spatially_variable_metabolites(
  grid[, setdiff(names(grid), "mz_constant")], n_perm = 49, seed = 14
)
geary <- compute_spatially_variable_metabolites_geary(
  grid[, setdiff(names(grid), "mz_constant")], n_perm = 49, seed = 14
)
comparison <- compare_svg_methods(
  moran, geary, batch[batch$feature != "mz_constant", ]
)
stopifnot(nrow(comparison) == 2L)
stopifnot(all(comparison$n_methods_tested == 3L))
stopifnot(nzchar(attr(comparison, "interpretation")))
