\name{hmrf_clustering}
\alias{cluster_pixels_hmrf}
\alias{cluster_diagnostics_hmrf}
\title{Potts-Regularized Spatial Clustering with Sequential ICM}
\description{
Fits a Gaussian/Potts-inspired spatially regularized k-means objective using
sequential iterated conditional modes (ICM) label updates and centroid updates.
This is a local objective optimizer, not full Bayesian HMRF inference.
}
\usage{
cluster_pixels_hmrf(pixel_matrix, k = 3, beta = 1,
  neighbor_method = c("queen", "rook", "distance"),
  distance_threshold = NULL, neighbors = NULL, scale_features = TRUE,
  data_term_scale = c("per_feature", "raw"), pca_components = NULL,
  max_iter = 20, update_order = c("random", "fixed"),
  init_nstart = 25, init_iter_max = 100, init_algorithm = "Lloyd",
  seed = 1, energy_tolerance = 1e-10)

cluster_diagnostics_hmrf(pixel_matrix, k = 3,
  beta_grid = c(0, 0.5, 1, 2, 4),
  neighbor_method = c("queen", "rook", "distance"),
  distance_threshold = NULL, scale_features = TRUE,
  data_term_scale = c("per_feature", "raw"), pca_components = NULL,
  max_iter = 20, update_order = c("random", "fixed"), seed = 1)
}
\arguments{
  \item{pixel_matrix}{Pixel table containing \code{x}, \code{y}, and \code{mz_} features.}
  \item{k}{Number of clusters.}
  \item{beta}{Non-negative Potts disagreement penalty.}
  \item{neighbor_method}{Rook, Queen, or distance-threshold neighborhood.}
  \item{distance_threshold}{Positive threshold for distance neighborhoods.}
  \item{neighbors}{Optional symmetric graph from \code{build_spatial_neighbors()}.}
  \item{scale_features}{Feature-wise z-score before PCA and clustering.}
  \item{data_term_scale}{Divide squared distance by the retained feature dimension or leave it raw.}
  \item{pca_components}{Optional retained principal components.}
  \item{max_iter}{Maximum ICM sweeps.}
  \item{update_order}{Random seeded or fixed node order within each sweep.}
  \item{init_nstart}{K-means initialization starts.}
  \item{init_iter_max}{Maximum initialization iterations.}
  \item{init_algorithm}{Initialization k-means algorithm.}
  \item{seed}{Optional random seed.}
  \item{energy_tolerance}{Numerical tolerance for ties and energy auditing.}
  \item{beta_grid}{Non-negative beta values for sensitivity analysis.}
}
\details{
The minimized objective is the within-cluster squared distance, optionally
divided by feature dimension, plus beta times the weighted number of
disagreeing undirected neighbor pairs. Labels are updated one node at a time;
cluster centers are recomputed after each sweep. Moves that would empty a
cluster are rejected.

Beta has no universal scale. It depends on preprocessing, retained dimension,
and graph degree. Adjacent-pair agreement is a smoothness diagnostic rather
than an accuracy measure and must not be used alone to select beta.
}
\value{
\code{cluster_pixels_hmrf()} returns clustered pixels, final centers and fit,
the graph, an objective trace, convergence information, and
\code{full_bayesian_hmrf = FALSE}. \code{cluster_diagnostics_hmrf()} returns
one diagnostic row per beta.
}
\references{
Besag J. On the statistical analysis of dirty pictures. Journal of the Royal
Statistical Society Series B 48, 259--302 (1986).
}
