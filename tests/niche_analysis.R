if (nzchar(Sys.getenv("_R_CHECK_PACKAGE_NAME_"))) {
  library(SpatialOmicsMSI)
} else {
  source(file.path(getwd(), "R", "msi_pipeline.R"))
}

niche_input <- data.frame(
  pixel_id = paste0("p", 1:8),
  region_id = rep(c("section_a", "section_b"), each = 4),
  x = rep(0:3, 2),
  y = c(rep(0, 4), rep(100, 4)),
  cell_type = rep(c("T", "T", "Tumor", "Tumor"), 2),
  stringsAsFactors = FALSE
)

composition <- compute_neighborhood_composition(
  niche_input$cell_type,
  niche_input$x,
  niche_input$y,
  k = 2,
  pixel_id = niche_input$pixel_id,
  region_id = niche_input$region_id
)
stopifnot(inherits(composition, "neighborhood_composition"))
stopifnot(all(rowSums(composition$counts) == 2L))
stopifnot(all(rowSums(composition$matrix[grep("^composition__", names(composition$matrix))]) == 1))
stopifnot(all(composition$neighbor_indices[1:4, ] <= 4L))
stopifnot(all(composition$neighbor_indices[5:8, ] >= 5L))
stopifnot(identical(composition$neighbor_indices[, 1], 1:8))

# Equal-distance ties are resolved by input row number and duplicate coordinates
# remain valid because the focal position is tracked by row identity.
tied <- compute_neighborhood_composition(
  c("A", "B", "C"), c(0, 0, 0), c(0, 0, 0), k = 2
)
stopifnot(identical(unname(tied$neighbor_indices[3, ]), c(3L, 1L)))

without_self <- compute_neighborhood_composition(
  niche_input$cell_type,
  niche_input$x,
  niche_input$y,
  k = 2,
  region_id = niche_input$region_id,
  include_self = FALSE
)
stopifnot(all(vapply(seq_len(8), function(i) !i %in% without_self$neighbor_indices[i, ], logical(1))))

niches_a <- define_niches(composition, k_niches = 2, seed = 19)
niches_b <- define_niches(composition, k_niches = 2, seed = 19)
stopifnot(identical(niches_a$matrix$niche_id, niches_b$matrix$niche_id))
stopifnot(nrow(niches_a$centers) == 2L)
stopifnot(grepl("arbitrary", niches_a$interpretation, fixed = TRUE))

too_small_failed <- tryCatch({
  compute_neighborhood_composition(c("T", "Tumor"), 1:2, c(1, 1), k = 2, include_self = FALSE)
  FALSE
}, error = function(e) TRUE)
stopifnot(too_small_failed)
