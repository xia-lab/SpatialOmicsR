testthat::test_that("regression: test-cell_type_labels", {
if (nzchar(Sys.getenv("_R_CHECK_PACKAGE_NAME_"))) {
  library(SpatialOmicsMSI)
} else {
  source(file.path(getwd(), "R", "msi_pipeline.R"))
  source(file.path(getwd(), "R", "spatial_registration.R"))
  source(file.path(getwd(), "R", "cell_type_labels.R"))
}

counts <- matrix(c(
  9, 1, 8, 2, 7, 3, 6, 4, 5, 5,
  1, 9, 2, 8, 3, 7, 4, 6, 5, 5,
  5, 5, 5, 5, 5, 5, 5, 5, 5, 5
), nrow = 3, byrow = TRUE)
colnames(counts) <- paste0("g", 1:10)
rownames(counts) <- c("s1", "s2", "s3")
reference <- cbind(
  type_a = c(9, 1, 8, 2, 7, 3, 6, 4, 5, 5),
  type_b = c(1, 9, 2, 8, 3, 7, 4, 6, 5, 5)
)
rownames(reference) <- colnames(counts)
spots <- data.frame(spot_id = c("s2", "s1", "s3"), x = 1:3, y = 0, section_id = "sec1")
cosine <- suppressWarnings(deconvolve_spatial_transcriptomics(
  counts, spots, reference, method = "cosine_label_transfer",
  region_column = "section_id"
))
stopifnot(inherits(cosine, "cell_type_labels"))
stopifnot(identical(cosine$labels$object_id, spots$spot_id))
stopifnot(all(abs(rowSums(cosine$proportions) - 1) < 1e-12))
stopifnot(identical(cosine$labels$label_mode, rep("similarity_score", 3)))

soft_proportions <- rbind(c(1, 0), c(0.5, 0.5), c(0, 1))
colnames(soft_proportions) <- c("a", "b")
soft <- compute_neighborhood_composition_soft(
  soft_proportions, x = 1:3, y = c(0, 0, 0),
  k = 2, region_id = rep("sec1", 3)
)
stopifnot(isTRUE(all.equal(unname(as.numeric(soft$matrix[1, c("composition__a", "composition__b")])), c(0.75, 0.25))))

labels <- prepare_cell_type_labels(data.frame(
  object_id = c("a", "b"), x = c(0, 2), y = c(0, 0), region_id = "sec1",
  cell_type = c("A", "B"), confidence = c(0.8, 0.9), source_method = "external"
))
pixels <- data.frame(pixel_id = 1:3, x = 0:2, y = 0, section_id = "sec1")
transferred <- transfer_cell_type_to_msi_pixels(
  pixels, labels, method = "nearest_neighbor", max_distance = 0.6
)
stopifnot(identical(transferred$matrix$cell_type, c("A", NA_character_, "B")))
stopifnot(transferred$qc$n_matched == 2L)

failed_without_distance <- tryCatch({
  transfer_cell_type_to_msi_pixels(pixels, labels, method = "nearest_neighbor")
  FALSE
}, error = function(e) TRUE)
stopifnot(failed_without_distance)

if (requireNamespace("EBImage", quietly = TRUE)) {
  image_path <- system.file("images", "nuclei.tif", package = "EBImage")
  morphology <- segment_cells_from_histology(
    image_path, min_nucleus_size = 5, max_nucleus_size = 5000,
    n_morphology_classes = 2, seed = 3
  )
  stopifnot(inherits(morphology, "cell_type_labels"))
  stopifnot(all(grepl("^morphology_class_", morphology$labels$cell_type)))
  stopifnot(all(morphology$labels$confidence >= 0.5 & morphology$labels$confidence <= 1))
}

  testthat::succeed()
})
