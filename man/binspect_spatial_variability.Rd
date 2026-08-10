\name{binspect_spatial_variability}
\alias{compute_binspect_feature}
\alias{compute_spatially_variable_metabolites_binspect}
\alias{compare_svg_methods}
\title{binSpect-Inspired Spatial Feature Enrichment}
\description{
Binarizes each MSI feature and tests whether high-state pixels are enriched
among spatial neighbors. The default uses a node-label permutation test;
an explicitly labeled Fisher approximation is available for compatibility
with the contingency-table style of Giotto binSpect.
}
\usage{
compute_binspect_feature(values, neighbors,
  bin_method = c("kmeans", "rank"), percentage_rank = 30,
  inference = c("permutation", "fisher_approx"), n_perm = 199,
  kmeans_nstart = 10, kmeans_iter_max = 100,
  hub_min_neighbors = 1, seed = NULL)

compute_spatially_variable_metabolites_binspect(pixel_matrix,
  coordinates = NULL, x_col = "x", y_col = "y",
  fcols = feature_columns(pixel_matrix),
  bin_method = c("kmeans", "rank"), percentage_rank = 30,
  inference = c("permutation", "fisher_approx"), n_perm = 199,
  p_adjust_method = "BH", seed = NULL, neighbors = NULL,
  neighbor_method = c("rook", "queen", "distance"),
  distance_threshold = NULL, weights = "binary", symmetric = TRUE)

compare_svg_methods(moran_result, geary_result, binspect_result,
  alpha = 0.05)
}
\arguments{
  \item{values}{Numeric feature value for every graph node.}
  \item{neighbors}{Graph returned by \code{build_spatial_neighbors()}.}
  \item{bin_method}{Two-means or upper-rank binarization.}
  \item{percentage_rank}{Percentage targeted as high under rank binarization; ties may increase the realized fraction.}
  \item{inference}{Node-label permutation test or Fisher contingency-table approximation.}
  \item{n_perm}{Number of permutations.}
  \item{kmeans_nstart}{Number of starts for two-means binarization.}
  \item{kmeans_iter_max}{Maximum two-means iterations.}
  \item{hub_min_neighbors}{Minimum high-state neighbors for a high-state pixel to count as a hub.}
  \item{seed}{Optional random seed.}
  \item{pixel_matrix}{Pixel table with \code{pixel_id} and \code{mz_} features.}
  \item{coordinates}{Optional separate coordinate table keyed by \code{pixel_id}.}
  \item{x_col, y_col}{Coordinate column names.}
  \item{fcols}{Features to test.}
  \item{p_adjust_method}{Method passed to \code{stats::p.adjust()}.}
  \item{neighbor_method}{Rook, Queen, or distance-threshold graph.}
  \item{distance_threshold}{Positive distance threshold.}
  \item{weights}{Weight type passed to graph construction.}
  \item{symmetric}{Whether a newly constructed graph is symmetric.}
  \item{moran_result, geary_result, binspect_result}{One-row-per-feature results.}
  \item{alpha}{Adjusted-p-value threshold used for comparison labels.}
}
\details{
The permutation test fixes the number of high pixels and randomly reallocates
their labels over effective graph nodes. Its upper-tail statistic is the total
weight of high--high undirected edges. The p-value uses the plus-one correction.

The Fisher approximation constructs a symmetric 2 by 2 edge table. Because
edges share nodes and each undirected edge contributes both orientations, its
observations are not independent; its p-value must not be described as an
exact graph-level test. Constant or otherwise untestable features are retained
in batch output with a status and error message.
}
\value{
The core function returns a named list. Batch and comparison functions return
data frames with adjusted p-values and audit columns.
}
\references{
Dries R, et al. Giotto, a toolbox for integrative analysis and visualization
of spatial expression data. Genome Biology 22, 78 (2021).
}
