\name{niche_analysis}
\alias{compute_neighborhood_composition}
\alias{define_niches}
\title{Cell-Type Composition Neighborhoods and Niche Clustering}
\description{
Builds fixed-size, region-aware spatial windows, represents each position by
the cell-type composition of its window, and clusters those compositions with
k-means to obtain cellular-neighborhood (niche) labels.
}
\usage{
compute_neighborhood_composition(cell_type_labels, x, y, k = 10,
  cell_types = NULL, pixel_id = seq_along(cell_type_labels),
  region_id = rep("region_1", length(cell_type_labels)),
  include_self = TRUE)

define_niches(composition, k_niches, nstart = 25,
  iter.max = 100, algorithm = "Lloyd", seed = 1)
}
\arguments{
  \item{cell_type_labels}{Hard cell-type label for each spatial position.}
  \item{x, y}{Finite spatial coordinates.}
  \item{k}{Fixed number of positions in each neighborhood window.}
  \item{cell_types}{Optional complete, ordered set of cell-type labels.}
  \item{pixel_id}{Unique position identifiers.}
  \item{region_id}{Tissue region, section, or sample identifier. Neighbors are
  found only within the same region.}
  \item{include_self}{Whether the focal position is one member of its fixed-size
  window. The default follows the reproduced Schurch et al. workflow.}
  \item{composition}{Object returned by
  \code{compute_neighborhood_composition()}, or a data frame whose composition
  columns start with \code{composition__}.}
  \item{k_niches}{Number of k-means niche clusters. This is dataset-specific.}
  \item{nstart, iter.max, algorithm}{Arguments passed to \code{stats::kmeans()}.}
  \item{seed}{Optional random seed.}
}
\details{
Equal-size windows make count vectors and proportion vectors differ only by
the constant \code{k}; therefore they give the same Euclidean k-means
partition. Both are retained: proportions in \code{matrix} and integer counts
in \code{counts}. Distance ties are resolved deterministically by input row
index. The implementation is dependency-free and uses an exact, quadratic-time
nearest-neighbor search, so a dedicated nearest-neighbor backend may be needed
for very large data sets.

Do not mix different tissue sections in one nearest-neighbor search. Supply a
meaningful \code{region_id}. Niche labels are arbitrary cluster identifiers;
interpret niches from their composition centers rather than from the numbers.
Hard labels transferred from another modality should be checked for registration
and deconvolution uncertainty before biological interpretation.
}
\value{
\code{compute_neighborhood_composition()} returns a
\code{neighborhood_composition} list containing the composition table, count
matrix, row-index neighbor matrix, cell-type/column mapping, and settings.
\code{define_niches()} returns the labeled composition table, k-means fit and
centers, settings, and an interpretation note.
}
\references{
Schurch CM, et al. Coordinated Cellular Neighborhoods Orchestrate Antitumoral
Immunity at the Colorectal Cancer Invasive Front. Cell 182, 1341--1359 (2020).
\doi{10.1016/j.cell.2020.07.005}
}
