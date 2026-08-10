\name{neighborhood_enrichment}
\alias{compute_neighborhood_enrichment}
\alias{group_cell_types_from_enrichment}
\title{Permutation-Based Cell-Type Neighborhood Enrichment}
\description{
Counts directed adjacency entries between cell types on a symmetric spatial
graph, compares counts with within-region label permutations, and optionally
groups cell types by the Euclidean distance between their complete enrichment
profiles for heatmap organization.
}
\usage{
compute_neighborhood_enrichment(cell_type_labels, neighbors,
  n_perm = 1000, seed = NULL,
  region_id = rep("region_1", length(cell_type_labels)),
  p_adjust_method = "BH", return_permutations = FALSE)

group_cell_types_from_enrichment(enrichment_result, k_groups,
  method = "ward.D2", zero_sd_action = c("error", "zero"),
  exclude_diagonal = FALSE)
}
\arguments{
  \item{cell_type_labels}{Hard cell-type or morphology label per graph node.}
  \item{neighbors}{Symmetric graph returned by \code{build_spatial_neighbors()}.}
  \item{n_perm}{Number of label permutations; at least two.}
  \item{seed}{Optional random seed.}
  \item{region_id}{Section/library identifier. Labels are shuffled only within region.}
  \item{p_adjust_method}{Adjustment method passed to \code{stats::p.adjust()}.}
  \item{return_permutations}{Whether to retain the potentially large permutation array.}
  \item{enrichment_result}{Result from \code{compute_neighborhood_enrichment()}.}
  \item{k_groups}{Number of cell-type interaction-profile groups.}
  \item{method}{Linkage method passed to \code{stats::hclust()}.}
  \item{zero_sd_action}{Stop on undefined z-scores or explicitly replace them with zero for clustering.}
  \item{exclude_diagonal}{Replace self-type z-scores with zero before profile clustering.}
}
\details{
The count convention follows a symmetric connectivity matrix: every undirected
edge contributes two directed adjacency entries. Consequently, a same-type edge
adds two to the diagonal, while a mixed-type edge adds one to each symmetric
off-diagonal cell. Permutations preserve label abundance separately within each
region. Graphs containing cross-region edges are rejected.

The z-score uses the population standard deviation of permutation counts. A
zero-variance null distribution produces an undefined z-score and is retained
as \code{NA}, rather than silently being reported as no enrichment. Empirical
one- and two-sided p-values use a plus-one correction; BH correction is applied
once to unique unordered type pairs and mirrored to the symmetric entry.

Interaction-profile grouping computes Euclidean distances between z-score rows
before hierarchical clustering. It mirrors heatmap reordering logic; it does
not turn a similarity matrix into a Ward distance and does not define or assign
spatial niches.
}
\value{
\code{compute_neighborhood_enrichment()} returns observed and null matrices,
z-scores, empirical p-values, a long pair table, settings, and optional
permutations. \code{group_cell_types_from_enrichment()} returns interaction
profile groups, hierarchical fit, distances, profiles, and an interpretation.
}
\references{
Palla G, et al. Squidpy: a scalable framework for spatial omics analysis.
Nature Methods 19, 171--178 (2022).
\doi{10.1038/s41592-021-01358-2}
}
