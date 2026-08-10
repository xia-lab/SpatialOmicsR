\name{geary_spatial_variability}
\alias{compute_gearys_c_grid}
\alias{compute_spatially_variable_metabolites_geary}
\alias{compare_moran_geary}
\title{Global Geary's C Spatial-Variability Analysis}
\description{
Computes global Geary's C on an existing spatial-neighbor graph, performs
Monte Carlo random-label permutation tests, applies multiple-testing
correction across MSI features, and compares significance calls with Moran's I.
}
\usage{
compute_gearys_c_grid(values, x = NULL, y = NULL, n_perm = 0,
  alternative = c("less", "greater", "two.sided"), seed = NULL,
  neighbors = NULL, neighbor_method = c("rook", "queen", "distance"),
  distance_threshold = NULL, weights = "binary", symmetric = TRUE)

compute_spatially_variable_metabolites_geary(pixel_matrix,
  coordinates = NULL, x_col = "x", y_col = "y",
  fcols = feature_columns(pixel_matrix), n_perm = 199,
  alternative = c("less", "greater", "two.sided"),
  p_adjust_method = "BH", seed = NULL, neighbors = NULL,
  neighbor_method = c("rook", "queen", "distance"),
  distance_threshold = NULL, weights = "binary", symmetric = TRUE)

compare_moran_geary(moran_result, geary_result, alpha = 0.05)
}
\arguments{
  \item{values}{Numeric feature value for every graph node.}
  \item{x, y}{Spatial coordinates used when \code{neighbors} is absent.}
  \item{n_perm}{Number of random-label permutations.}
  \item{alternative}{Positive autocorrelation (less), negative autocorrelation (greater), or two-sided departure from one.}
  \item{seed}{Optional random seed.}
  \item{neighbors}{Graph returned by \code{build_spatial_neighbors()}.}
  \item{neighbor_method}{Rook, Queen, or distance-threshold graph.}
  \item{distance_threshold}{Positive distance threshold.}
  \item{weights}{Weight type passed to \code{build_spatial_neighbors()}.}
  \item{symmetric}{Whether a newly constructed graph is symmetric.}
  \item{pixel_matrix}{Pixel table with \code{pixel_id} and \code{mz_} features.}
  \item{coordinates}{Optional separate coordinate table keyed by \code{pixel_id}.}
  \item{x_col, y_col}{Coordinate column names.}
  \item{fcols}{Features to test.}
  \item{p_adjust_method}{Method passed to \code{stats::p.adjust()}.}
  \item{moran_result}{One-row-per-feature Moran result.}
  \item{geary_result}{One-row-per-feature Geary result.}
  \item{alpha}{Adjusted-p-value threshold used only for comparison labels.}
}
\details{
For symmetric weights represented once per undirected edge, Geary's C is
\deqn{C = \frac{n-1}{2\sum w_{ij}}
\frac{\sum w_{ij}(x_i-x_j)^2}{\sum_i(x_i-\bar{x})^2}.}
Its random-label expectation is one. Values below one indicate positive spatial
autocorrelation and values above one indicate negative spatial autocorrelation.
The default lower-tail test therefore targets positive spatial structure.

Global Geary's C emphasizes squared neighbor differences but is not a local
hotspot statistic. Moran/Geary agreement is a sensitivity analysis because the
two tests use the same data and graph; it is not independent validation.
}
\value{
The core function returns a named list. The batch and comparison functions
return data frames with adjusted p-values and graph diagnostics.
}
\references{
Geary RC. The contiguity ratio and statistical mapping. The Incorporated
Statistician 5, 115--145 (1954). \doi{10.2307/2986645}
}
