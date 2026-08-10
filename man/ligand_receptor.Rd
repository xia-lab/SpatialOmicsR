\name{ligand_receptor_analysis}
\alias{cell_proximity_enrichment}
\alias{expr_cell_cell_communication}
\alias{spat_cell_cell_communication}
\title{Cell-type proximity and ligand-receptor sensitivity analyses}
\description{
Computes stratified cell-type proximity enrichment and expression-only or
spatial ligand-receptor communication-potential sensitivity analyses.
}
\usage{
cell_proximity_enrichment(
  cell_type_labels, neighbors, permutation_strata = NULL,
  n_perm = 1000L, p_adjust_method = "BH",
  cross_stratum_action = c("error", "drop"), seed = 1234
)

expr_cell_cell_communication(
  expression_matrix, cell_type_labels, lr_database,
  ligand_col = "ligand", receptor_col = "receptor",
  complex_delimiter = "|", permutation_strata = NULL,
  expression_threshold = 0, min_expression_fraction = 0.1,
  min_cells_per_type = 5L, n_perm = 1000L,
  p_adjust_method = "BH", seed = 1234
)

spat_cell_cell_communication(
  expression_matrix, cell_type_labels, neighbors, lr_database,
  ligand_col = "ligand", receptor_col = "receptor",
  complex_delimiter = "|", permutation_strata = NULL,
  expression_threshold = 0, min_expression_fraction = 0.1,
  min_interacting_cells = 5L, n_perm = 1000L,
  p_adjust_method = "BH",
  cross_stratum_action = c("error", "drop"), seed = 1234
)
}
\arguments{
  \item{cell_type_labels}{One cell-type label per node.}
  \item{neighbors}{A graph returned by \code{build_spatial_neighbors()}.}
  \item{permutation_strata}{Optional vector, matrix, or data frame defining
  exchangeability blocks such as subject by section.}
  \item{n_perm}{Positive number of conditional randomizations.}
  \item{p_adjust_method}{Method passed to \code{p.adjust}.}
  \item{cross_stratum_action}{Whether graph edges joining distinct strata are
  errors or are explicitly removed.}
  \item{seed}{Optional local seed. The caller's random state is restored.}
  \item{expression_matrix}{Numeric nodes-by-genes normalized expression matrix
  with unique gene column names.}
  \item{lr_database}{Data frame whose rows are known ligand-receptor pairs.}
  \item{ligand_col,receptor_col}{Columns containing ligand and receptor names.}
  \item{complex_delimiter}{Fixed separator between required subunits. All
  subunits must be present; per-node complex expression is their minimum.}
  \item{expression_threshold}{Value above which a gene or complex is expressed.}
  \item{min_expression_fraction}{Minimum expressing fraction required on both
  sides before the combined score is nonzero.}
  \item{min_cells_per_type}{Minimum cells in each sender and receiver type.}
  \item{min_interacting_cells}{Minimum unique interacting cells on each side.}
}
\value{
\code{cell_proximity_enrichment()} returns count and ratio matrices, a long
table with upper, lower, and two-sided Monte Carlo p-values, and settings.
The communication functions return one row per database pair and directed
cell-type combination, including expression summaries, expressing fractions,
\code{LR_expr}, \code{rand_expr}, \code{av_diff}, cell counts, Monte Carlo
p-values, and globally adjusted p-values. Spatial output also includes edge
counts, log2 fold change, and a Giotto-style PI summary.
}
\details{
The proximity null shuffles labels within each stratum on a fixed graph. The
expression-only null similarly tests cell-type-specific expression of known
paired database rows.

The spatial null keeps expression, cell types, strata, and participating-cell
counts fixed. Within each stratum and cell type, it samples the same number of
cells as found on observed cross-type edges. It asks whether spatially
interacting subsets have greater ligand-receptor expression than matched cells
of the same types.

These are Giotto-inspired conditional-randomization sensitivity analyses, not
an exact reimplementation of Giotto. Significance represents communication
potential supported by expression and proximity; it does not establish protein
binding, receptor activation, causal signaling, or population-level biological
replication. Mixed-cell spots must be described as spot-level communication.
}
