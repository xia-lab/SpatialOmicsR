checking_installed_package <- nzchar(Sys.getenv("_R_CHECK_PACKAGE_NAME_"))
if (checking_installed_package) {
  library(SpatialOmicsMSI)
} else {
  source(file.path(getwd(), "R", "msi_pipeline.R"))
  source(file.path(getwd(), "R", "roi_selection.R"))
  source(file.path(getwd(), "R", "roi_annotation_import.R"))
}

grid <- expand.grid(x = 0:10, y = 0:10)
grid$qc_pass <- TRUE
grid$histology_cluster <- ifelse(grid$x <= 5, "target", "other")

# One ROI part with an exclusion hole.
outer <- data.frame(
  roi_id = "with_hole", polygon_part_id = 1L, ring_id = 1L,
  ring_role = "outer", vertex_order = 1:5,
  x = c(1, 9, 9, 1, 1), y = c(1, 1, 9, 9, 1)
)
hole <- data.frame(
  roi_id = "with_hole", polygon_part_id = 1L, ring_id = 2L,
  ring_role = "hole", vertex_order = 1:5,
  x = c(4, 6, 6, 4, 4), y = c(4, 4, 6, 6, 4)
)
with_hole <- select_rois(
  grid, selection_mode = "manual", manual_method = "polygon",
  polygon_vertices = rbind(outer, hole)
)
stopifnot(!any(with_hole$selected_pixels$x >= 4 & with_hole$selected_pixels$x <= 6 &
                 with_hole$selected_pixels$y >= 4 & with_hole$selected_pixels$y <= 6))
stopifnot(any(with_hole$selected_pixels$x == 2 & with_hole$selected_pixels$y == 2))

# Two disjoint polygon parts retain one ROI identity without a connecting edge.
part_a <- transform(outer, roi_id = "multipart", ring_id = 1L,
                    polygon_part_id = 1L, x = c(0, 2, 2, 0, 0), y = c(0, 0, 2, 2, 0))
part_b <- transform(outer, roi_id = "multipart", ring_id = 1L,
                    polygon_part_id = 2L, x = c(8, 10, 10, 8, 8), y = c(8, 8, 10, 10, 8))
multipart <- select_rois(
  grid, "manual", "polygon", polygon_vertices = rbind(part_a, part_b)
)
stopifnot(identical(unique(multipart$selected_pixels$roi_id), "multipart"))
stopifnot(any(multipart$selected_pixels$x <= 2), any(multipart$selected_pixels$x >= 8))
stopifnot(!any(multipart$selected_pixels$x == 5 & multipart$selected_pixels$y == 5))

# Missing labels count against target enrichment rather than inflating purity.
candidate <- structure(list(
  candidates = data.frame(candidate_id = "c1", x = 1, y = 1),
  membership = list(1:4)
), class = "roi_candidates")
score_data <- data.frame(
  x = 1:4, y = 1, qc_pass = TRUE,
  label = c("target", NA, NA, NA)
)
target_score <- score_roi_indices(
  score_data, candidate, 1L, cluster_column = "label",
  objective = "target_enrichment", target_cluster = "target"
)
stopifnot(target_score[["target_purity_among_labeled"]] == 1)
stopifnot(target_score[["target_fraction_among_valid"]] == 0.25)
stopifnot(target_score[["label_coverage"]] == 0.25)

# Domain components respect missing grid cells and rook connectivity.
domain_grid <- rbind(
  data.frame(x = 0:2, y = 0, label = "a", qc_pass = TRUE),
  data.frame(x = 5:7, y = 0, label = "a", qc_pass = TRUE),
  data.frame(x = 0:2, y = 2, label = "b", qc_pass = TRUE)
)
domain_candidates <- generate_roi_candidates_from_domains(
  domain_grid, cluster_column = "label", min_domain_size = 3,
  connectivity = "rook", topology_x_step = 1, topology_y_step = 1
)
stopifnot(nrow(domain_candidates$candidates) == 3L)
stopifnot(identical(sort(lengths(domain_candidates$membership)), c(3L, 3L, 3L)))

# Only pixels whose actual label pair is accepted are called corroborated.
agreement <- corroborate_cluster_labels(
  label_a = rep("a1", 10),
  label_b = c(rep("b1", 6), rep("b2", 4)),
  min_cooccurrence_fraction = 0.5,
  mapping = "one_way"
)
stopifnot(sum(agreement$corroborated) == 6L)
stopifnot(all(agreement$combined_label[7:10] == "uncorroborated"))

# QuPath Polygon and MultiPolygon coordinates are normalized to the common ROI schema.
if (requireNamespace("jsonlite", quietly = TRUE)) {
  geojson_path <- tempfile(fileext = ".geojson")
  writeLines(paste0(
    '{"type":"FeatureCollection","features":[',
    '{"type":"Feature","id":"q1","properties":{"name":"tumor","classification":{"name":"Tumor"}},',
    '"geometry":{"type":"Polygon","coordinates":[[[0,0],[3,0],[3,3],[0,3],[0,0]],',
    '[[1,1],[2,1],[2,2],[1,2],[1,1]]]}}]}'
  ), geojson_path)
  imported <- import_qupath_geojson(geojson_path, section_id = "s1")
  stopifnot(all(c("polygon_part_id", "ring_id", "ring_role", "coordinate_unit") %in% names(imported)))
  stopifnot(identical(unique(imported$ring_role), c("outer", "hole")))
  stopifnot(all(imported$coordinate_unit == "full_resolution_pixel"))
}

# LMD validation reports violations without altering input geometry.
lmd_a <- transform(part_a, roi_id = "a")
lmd_b <- transform(part_a, roi_id = "b", x = x + 2.5)
lmd_qc <- validate_lmd_shapes(
  rbind(lmd_a, lmd_b), max_area = 10, min_spacing = 1,
  min_collectable_area = 1, coordinate_unit = "um"
)
stopifnot(nrow(lmd_qc) == 2L)
stopifnot(all(lmd_qc$area == 4))
stopifnot(all(lmd_qc$spacing_below_minimum))
stopifnot(!any(lmd_qc$qc_pass))
stopifnot(grepl("Read-only", attr(lmd_qc, "scope"), fixed = TRUE))

cat("ROI_UPGRADE_TEST_OK=TRUE\n")
