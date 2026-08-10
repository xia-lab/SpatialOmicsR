if (nzchar(Sys.getenv("_R_CHECK_PACKAGE_NAME_"))) {
  library(SpatialOmicsMSI)
} else {
  source(file.path(getwd(), "R", "msi_pipeline.R"))
  source(file.path(getwd(), "R", "ligand_receptor.R"))
}

# Two independent tissue fields with the same cell-type composition.
field <- rep(c("subject_1::section_1", "subject_2::section_1"), each = 12L)
x <- c(rep(1:6, each = 2L), rep(101:106, each = 2L))
y <- rep(rep(1:2, times = 6L), 2L)
cell_type <- rep(c("A", "B"), 12L)
neighbors <- build_spatial_neighbors(x, y, method = "queen")

expression <- cbind(
  L = ifelse(cell_type == "A", 5, 0),
  L2 = ifelse(cell_type == "A", 4, 0),
  R = ifelse(cell_type == "B", 6, 0),
  background = rep(c(0, 1, 2), length.out = length(cell_type))
)
lr_database <- data.frame(
  ligand = c("L", "L|L2", "missing"),
  receptor = c("R", "R", "R"),
  stringsAsFactors = FALSE
)

proximity <- cell_proximity_enrichment(
  cell_type, neighbors, permutation_strata = field,
  n_perm = 49, seed = 7
)
stopifnot(is.list(proximity))
stopifnot(nrow(proximity$enrichment_table) == 3L)
stopifnot(all(c(
  "p_enrichment", "p_depletion", "p_two_sided",
  "adj_p_enrichment", "adj_p_depletion", "adj_p_two_sided"
) %in% names(proximity$enrichment_table)))
stopifnot(proximity$settings$stratified_permutation)
stopifnot(proximity$settings$n_strata == 2L)
stopifnot(proximity$settings$p_value_resolution == 1 / 50)

# A subject/section data frame is combined column-wise into exchangeability blocks.
strata_frame <- data.frame(
  subject = rep(c("subject_1", "subject_2"), each = 12L),
  section = "section_1",
  stringsAsFactors = FALSE
)
proximity_frame <- cell_proximity_enrichment(
  cell_type, neighbors, permutation_strata = strata_frame,
  n_perm = 9, seed = 7
)
stopifnot(proximity_frame$settings$n_strata == 2L)

expression_only <- expr_cell_cell_communication(
  expression, cell_type, lr_database,
  permutation_strata = field, min_cells_per_type = 5,
  n_perm = 49, seed = 7
)
stopifnot(nrow(expression_only) == 8L) # 2 complete LR rows x 4 directed type pairs
stopifnot(!any(expression_only$ligand == "missing"))
stopifnot(all(c("LR_expr", "rand_expr", "av_diff", "lig_nr", "rec_nr") %in% names(expression_only)))
signal <- expression_only[
  expression_only$ligand == "L" & expression_only$receptor == "R" &
    expression_only$lig_cell_type == "A" & expression_only$rec_cell_type == "B",
  , drop = FALSE
]
stopifnot(nrow(signal) == 1L, signal$LR_expr == 11)
stopifnot(signal$lig_expression_fraction == 1, signal$rec_expression_fraction == 1)
stopifnot(signal$p_enrichment <= 0.05)

spatial <- spat_cell_cell_communication(
  expression, cell_type, neighbors, lr_database,
  permutation_strata = field, min_interacting_cells = 5,
  n_perm = 49, seed = 7
)
stopifnot(nrow(spatial) == 8L)
stopifnot(all(c("log2fc", "PI", "n_edges", "status") %in% names(spatial)))
spatial_signal <- spatial[
  spatial$ligand == "L" & spatial$receptor == "R" &
    spatial$lig_cell_type == "A" & spatial$rec_cell_type == "B",
  , drop = FALSE
]
stopifnot(nrow(spatial_signal) == 1L)
stopifnot(spatial_signal$status == "tested")
stopifnot(spatial_signal$lig_nr >= 5L, spatial_signal$rec_nr >= 5L)
stopifnot(attr(spatial, "settings")$stratified_permutation)
stopifnot(grepl("same number of cells", attr(spatial, "settings")$null_model, fixed = TRUE))

# A database is a set of paired rows, not the ligand x receptor Cartesian product.
pair_only <- expr_cell_cell_communication(
  expression, cell_type,
  data.frame(ligand = "L", receptor = "R"),
  permutation_strata = field, min_cells_per_type = 5,
  n_perm = 9, seed = 1
)
stopifnot(nrow(pair_only) == 4L)

# Reproducible local seeding must not mutate the caller's RNG stream.
set.seed(99)
before <- .Random.seed
invisible(cell_proximity_enrichment(cell_type, neighbors, field, n_perm = 9, seed = 4))
stopifnot(identical(before, .Random.seed))

# Cross-field edges are rejected by default and may only be dropped explicitly.
cross_neighbors <- build_spatial_neighbors(1:4, rep(1, 4), method = "rook")
cross_error <- tryCatch(
  cell_proximity_enrichment(
    c("A", "B", "A", "B"), cross_neighbors,
    permutation_strata = c("field_1", "field_1", "field_2", "field_2"),
    n_perm = 9
  ),
  error = function(e) conditionMessage(e)
)
stopifnot(grepl("crossing permutation strata", cross_error, fixed = TRUE))
cross_dropped <- cell_proximity_enrichment(
  c("A", "B", "A", "B"), cross_neighbors,
  permutation_strata = c("field_1", "field_1", "field_2", "field_2"),
  n_perm = 9, cross_stratum_action = "drop"
)
stopifnot(cross_dropped$settings$n_cross_stratum_edges_dropped == 1L)

cat("LIGAND_RECEPTOR_TEST_OK=TRUE\n")
