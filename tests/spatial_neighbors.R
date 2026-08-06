checking_installed_package <- nzchar(Sys.getenv("_R_CHECK_PACKAGE_NAME_"))
if (checking_installed_package) {
  library(SpatialOmicsMSI)
} else {
  project_root <- normalizePath(getwd(), mustWork = TRUE)
  source(file.path(project_root, "R", "msi_pipeline.R"))
}

expect_error <- function(expression, pattern = NULL) {
  error <- tryCatch(expression, error = function(condition) condition)
  stopifnot(inherits(error, "error"))
  if (!is.null(pattern)) stopifnot(grepl(pattern, conditionMessage(error), fixed = TRUE))
  invisible(error)
}

# Complete 3 x 3 grid: 12 rook edges and 20 queen edges.
grid <- expand.grid(x = 0:2, y = 0:2)
rook <- build_spatial_neighbors(
  grid$x, grid$y,
  method = "rook", distance_threshold = 1,
  weights = "binary", symmetric = TRUE
)
queen <- build_spatial_neighbors(
  grid$x, grid$y,
  method = "queen", distance_threshold = sqrt(2),
  weights = "binary", symmetric = TRUE
)
stopifnot(nrow(rook$undirected_edges) == 12L)
stopifnot(nrow(queen$undirected_edges) == 20L)
stopifnot(rook$summary$n_isolated == 0L, queen$summary$n_isolated == 0L)
stopifnot(!rook$summary$includes_diagonal_neighbors)
stopifnot(queen$summary$includes_diagonal_neighbors)

# Symmetric edge list contains every reverse edge and never a self-loop.
stopifnot(!any(queen$edges$from == queen$edges$to))
edge_keys <- paste(queen$edges$from, queen$edges$to, sep = "->")
reverse_keys <- paste(queen$edges$to, queen$edges$from, sep = "->")
stopifnot(all(reverse_keys %in% edge_keys))
stopifnot(nrow(queen$edges) == 2L * nrow(queen$undirected_edges))

# Purely diagonal pixels connect under queen but not rook.
diagonal_x <- c(0, 1)
diagonal_y <- c(0, 1)
diagonal_rook <- build_spatial_neighbors(diagonal_x, diagonal_y, method = "rook")
diagonal_queen <- build_spatial_neighbors(diagonal_x, diagonal_y, method = "queen")
stopifnot(nrow(diagonal_rook$undirected_edges) == 0L)
stopifnot(nrow(diagonal_queen$undirected_edges) == 1L)

# Explicit distance threshold is supported and remains binary/symmetric.
distance_graph <- build_spatial_neighbors(
  c(0, 1, 3), c(0, 1, 0),
  method = "distance", distance_threshold = sqrt(2),
  weights = "binary", symmetric = TRUE
)
stopifnot(nrow(distance_graph$undirected_edges) == 1L)
stopifnot(distance_graph$summary$n_isolated == 1L)

# An isolated extreme value cannot change Moran's I on the connected subgraph.
connected_coords <- expand.grid(x = 0:1, y = 0:1)
connected_values <- c(1, 2, 3, 4)
connected_graph <- build_spatial_neighbors(connected_coords$x, connected_coords$y, method = "rook")
connected_stat <- compute_morans_i_grid(
  connected_values, neighbors = connected_graph,
  n_perm = 99, alternative = "two.sided", seed = 77
)
with_isolate_graph <- build_spatial_neighbors(
  c(connected_coords$x, 100), c(connected_coords$y, 100), method = "rook"
)
with_isolate_stat <- compute_morans_i_grid(
  c(connected_values, 1e12), neighbors = with_isolate_graph,
  n_perm = 99, alternative = "two.sided", seed = 77
)
stopifnot(identical(connected_stat$I, with_isolate_stat$I))
stopifnot(identical(connected_stat$p_value, with_isolate_stat$p_value))
stopifnot(with_isolate_stat$n_total == 5L)
stopifnot(with_isolate_stat$n_effective == 4L)
stopifnot(with_isolate_stat$n_isolated == 1L)

# No-edge graphs must fail rather than return a fabricated numeric statistic.
no_edge_graph <- build_spatial_neighbors(c(0, 10, 20), c(0, 10, 20), method = "rook")
expect_error(
  compute_morans_i_grid(c(1, 2, 3), neighbors = no_edge_graph),
  "effective graph has no edges"
)

# Fixed seed gives exactly reproducible permutation results.
gradient_values <- grid$x + 2 * grid$y
seed_a <- compute_morans_i_grid(
  gradient_values, neighbors = queen,
  n_perm = 199, alternative = "two.sided", seed = 1234
)
seed_b <- compute_morans_i_grid(
  gradient_values, neighbors = queen,
  n_perm = 199, alternative = "two.sided", seed = 1234
)
stopifnot(identical(seed_a, seed_b))

diagnostics <- spatial_neighbor_diagnostics(queen)
stopifnot(all(c("summary", "degree_distribution", "component_sizes") %in% unique(diagnostics$section)))
stopifnot(all(c("n_nodes", "n_effective", "n_isolated", "n_edges_undirected") %in% diagnostics$metric))

cat("SPATIAL_NEIGHBOR_TEST_OK=TRUE\n")
