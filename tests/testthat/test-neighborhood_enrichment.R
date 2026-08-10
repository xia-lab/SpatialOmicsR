testthat::test_that("regression: test-neighborhood_enrichment", {
if (nzchar(Sys.getenv("_R_CHECK_PACKAGE_NAME_"))) {
  library(SpatialOmicsMSI)
} else {
  source(file.path(getwd(), "R", "msi_pipeline.R"))
  source(file.path(getwd(), "R", "neighborhood_enrichment.R"))
}

# A 2 x 3 rook grid with symmetric directed adjacency entries.
x <- rep(1:3, each = 2)
y <- rep(1:2, times = 3)
labels <- c("A", "A", "A", "B", "B", "B")
graph <- build_spatial_neighbors(x, y, method = "rook")
result_a <- compute_neighborhood_enrichment(labels, graph, n_perm = 99, seed = 12)
result_b <- compute_neighborhood_enrichment(labels, graph, n_perm = 99, seed = 12)
stopifnot(identical(result_a$observed_count, result_b$observed_count))
stopifnot(identical(result_a$z_score, result_b$z_score))
stopifnot(isTRUE(all.equal(result_a$observed_count, t(result_a$observed_count))))

# Each undirected same-type edge contributes both adjacency orientations.
edges <- graph$undirected_edges
aa_edges <- sum(labels[edges$from] == "A" & labels[edges$to] == "A")
stopifnot(result_a$observed_count["A", "A"] == 2 * aa_edges)
stopifnot(all(result_a$pair_table$adj_p_two_sided >= 0 & result_a$pair_table$adj_p_two_sided <= 1))

# Shuffling is region-stratified and cross-region graph edges are rejected.
region <- rep(c("r1", "r2"), each = 3)
cross_region_failed <- tryCatch({
  compute_neighborhood_enrichment(labels, graph, n_perm = 9, region_id = region)
  FALSE
}, error = function(e) TRUE)
stopifnot(cross_region_failed)

synthetic <- result_a
synthetic$cell_types <- c("A", "B", "C")
synthetic$z_score <- matrix(c(
  3, 2, -2,
  2, 3, -2,
  -2, -2, 4
), 3, 3, byrow = TRUE, dimnames = list(synthetic$cell_types, synthetic$cell_types))
grouped <- group_cell_types_from_enrichment(synthetic, k_groups = 2)
stopifnot(grouped$groups$interaction_profile_group[1] == grouped$groups$interaction_profile_group[2])
stopifnot(grouped$groups$interaction_profile_group[1] != grouped$groups$interaction_profile_group[3])
stopifnot(grepl("not spatial niches", grouped$interpretation, fixed = TRUE))

  testthat::succeed()
})
