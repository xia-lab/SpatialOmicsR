\name{spatial_clustering}
\alias{compute_neighborhood_average}
\alias{cluster_pixels_spatial}
\alias{cluster_diagnostics_spatial}
\title{BANKSY-Inspired Mean-Neighborhood Spatial Clustering}
\description{
Constructs a mean-neighborhood (H0) feature block, combines it with each
pixel's own features using BANKSY-style lambda weights, and applies optional
PCA followed by k-means. This is an H0-only approximation and not the complete
BANKSY algorithm: it does not compute the azimuthal Gabor feature, use the
default BANKSY k-nearest-neighbor Gaussian kernel, or perform Leiden clustering.
}
\usage{
compute_neighborhood_average(feature_matrix, neighbors,
  isolate_action = c("self", "error"))

cluster_pixels_spatial(pixel_matrix, k = 3, lambda = 0.5,
  neighbor_method = c("queen", "rook", "distance"),
  distance_threshold = NULL, neighbors = NULL,
  isolate_action = c("self", "error"), pca_components = NULL,
  nstart = 25, iter.max = 100, algorithm = "Lloyd", seed = 1)

cluster_diagnostics_spatial(pixel_matrix, k = 3,
  lambda_grid = seq(0, 1, by = 0.25),
  neighbor_method = c("queen", "rook", "distance"),
  distance_threshold = NULL, isolate_action = c("self", "error"),
  pca_components = NULL, nstart = 25, iter.max = 100,
  algorithm = "Lloyd", seed = 1)
}
\arguments{
  \item{feature_matrix}{Numeric matrix with spatial nodes in rows and features in columns.}
  \item{neighbors}{A symmetric graph returned by \code{build_spatial_neighbors()}.}
  \item{isolate_action}{Use an isolated node's own features or raise an error.}
  \item{pixel_matrix}{Pixel table containing \code{x}, \code{y}, and \code{mz_} feature columns.}
  \item{k}{Number of k-means clusters.}
  \item{lambda}{Spatial mixing weight in the interval zero to one.}
  \item{neighbor_method}{Rook, Queen, or distance-threshold neighborhood.}
  \item{distance_threshold}{Positive spatial distance threshold.}
  \item{pca_components}{Optional number of principal components.}
  \item{nstart}{Number of k-means random starts.}
  \item{iter.max}{Maximum k-means iterations.}
  \item{algorithm}{K-means algorithm passed to \code{stats::kmeans()}.}
  \item{seed}{Optional random seed.}
  \item{lambda_grid}{Lambda values evaluated by the diagnostic helper.}
}
\details{
The own and mean-neighborhood blocks are independently feature-wise z-scaled
and multiplied by \eqn{\sqrt{1-\lambda}} and \eqn{\sqrt{\lambda}}, respectively.
Consequently, lambda zero is nonspatial clustering of z-scaled own features;
it is not necessarily identical to \code{cluster_pixels()} when that function
receives unscaled features.

Adjacent-pair agreement is a descriptive smoothness statistic. It generally
favors smoother partitions and must not be used as the sole criterion for
choosing lambda or as evidence of clustering accuracy.
}
\value{
\code{compute_neighborhood_average()} returns a numeric matrix.
\code{cluster_pixels_spatial()} returns the clustered matrix, k-means fit,
neighbor graph, feature blocks, settings, and an explicit \code{full_banksy = FALSE}
flag. \code{cluster_diagnostics_spatial()} returns one diagnostic row per lambda.
}
\references{
Singhal V, et al. BANKSY unifies cell typing and tissue domain segmentation
for scalable spatial omics data analysis. Nature Genetics 56, 431--441 (2024).
\doi{10.1038/s41588-024-01664-3}
}
