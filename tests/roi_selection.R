if (nzchar(Sys.getenv("_R_CHECK_PACKAGE_NAME_"))) {
  library(SpatialOmicsMSI)
} else {
  source(file.path(getwd(), "R", "msi_pipeline.R"))
  source(file.path(getwd(), "R", "spatial_registration.R"))
  source(file.path(getwd(), "R", "roi_selection.R"))
}

grid <- expand.grid(x = 0:19, y = 0:19)
grid$pixel_id <- seq_len(nrow(grid))
grid$qc_pass <- !(grid$x < 2 & grid$y < 2)
grid$histology_cluster <- ifelse(grid$x < 10, 1, 2)
grid$mz_100 <- seq_len(nrow(grid))

automatic <- select_rois(
  grid,
  selection_mode = "automatic",
  roi_size = 6,
  n_candidates = 60,
  max_rois = 3,
  improvement_threshold = 0.01,
  beam_width = 8,
  candidate_pool = 30,
  seed = 3
)
stopifnot(
  nrow(automatic$selected_pixels) > 0,
  all(c("roi_score", "balance_score", "coverage_score", "size_score") %in%
        names(automatic$optimization$score))
)

geometry <- data.frame(
  roi_id = "box",
  shape = "rectangle",
  x_min = 1,
  x_max = 4,
  y_min = 2,
  y_max = 5
)
manual <- select_rois(grid, "manual", "geometry", roi_table = geometry)
stopifnot(nrow(manual$selected_pixels) == 16)

vertices <- data.frame(
  roi_id = "polygon",
  vertex_order = 1:4,
  x = c(1, 5, 5, 1),
  y = c(1, 1, 5, 5)
)
polygon <- select_rois(grid, "manual", "polygon", polygon_vertices = vertices)
stopifnot(nrow(polygon$selected_pixels) > 0)

combined <- select_rois(
  grid,
  "manual",
  "combined",
  roi_table = geometry,
  polygon_vertices = transform(vertices, roi_id = "polygon_2"),
  overlap = "first"
)
stopifnot(all(c("box", "polygon_2") %in% unique(combined$selected_pixels$roi_id)))

serial <- rbind(
  transform(grid, section_id = "A"),
  transform(grid, pixel_id = pixel_id + 1000, section_id = "B")
)
serial_geometry <- data.frame(
  section_id = c("A", "B"),
  roi_id = c("a", "b"),
  shape = "rectangle",
  x_min = c(0, 10),
  x_max = c(3, 13),
  y_min = 0,
  y_max = 3
)
serial_manual <- select_rois(
  serial,
  "manual",
  "geometry",
  roi_table = serial_geometry,
  section_column = "section_id"
)
stopifnot(
  all(serial_manual$selected_pixels$section_id[serial_manual$selected_pixels$roi_id == "a"] == "A"),
  all(serial_manual$selected_pixels$section_id[serial_manual$selected_pixels$roi_id == "b"] == "B")
)

labeled <- apply_roi_labels(grid, manual$selected_pixels)
samples <- sample_subregions(labeled, grid_size = 2, min_pixels = 1, grid_scope = "roi")
stopifnot(nrow(samples$sample_matrix) > 0)

control_points <- data.frame(
  histology_x = c(0, 1, 0, 2),
  histology_y = c(0, 0, 1, 2),
  msi_x = c(10, 12, 10, 14),
  msi_y = c(20, 20, 23, 26)
)
registration <- fit_histology_msi_registration(control_points)
registered_vertices <- transform_histology_coordinates(
  data.frame(roi_id = "he_polygon", vertex_order = 1:3, x = c(0, 1, 1), y = c(0, 0, 1)),
  registration
)
stopifnot(
  isTRUE(all.equal(registered_vertices$x, c(10, 12, 12), tolerance = 1e-10)),
  isTRUE(all.equal(registered_vertices$y, c(20, 20, 23), tolerance = 1e-10)),
  registration_diagnostics(registration)$rmse < 1e-10
)

# Control-point error cases
too_few_points <- data.frame(
  histology_x = c(0, 1),
  histology_y = c(0, 0),
  msi_x = c(10, 11),
  msi_y = c(20, 20)
)
err_too_few <- tryCatch(
  fit_histology_msi_registration(too_few_points),
  error = function(e) e
)
stopifnot(inherits(err_too_few, "error"))

collinear_points <- data.frame(
  histology_x = c(0, 1, 2),
  histology_y = c(0, 0, 0),
  msi_x = c(10, 11, 12),
  msi_y = c(20, 20, 20)
)
err_collinear <- tryCatch(
  fit_histology_msi_registration(collinear_points),
  error = function(e) e
)
stopifnot(inherits(err_collinear, "error"))

# Section-aware histology registration with polygon ROI coordinates
section_grid <- rbind(
  transform(grid, section_id = "A", x = x + 10, y = y + 20, pixel_id = pixel_id),
  transform(grid, section_id = "B", x = x + 20, y = y + 30, pixel_id = pixel_id + 1000)
)
control_points_section <- data.frame(
  section_id = rep(c("A", "B"), each = 3),
  histology_x = c(0, 1, 0, 0, 1, 0),
  histology_y = c(0, 0, 1, 0, 0, 1),
  msi_x = c(10, 11, 10, 20, 21, 20),
  msi_y = c(20, 20, 21, 30, 30, 31)
)
registration_section <- fit_histology_msi_registration(control_points_section, section_column = "section_id")
polygon_vertices <- data.frame(
  section_id = c("A", "A", "A", "A", "B", "B", "B", "B"),
  roi_id = c(rep("polyA", 4), rep("polyB", 4)),
  vertex_order = rep(1:4, 2),
  x = c(0, 2, 2, 0, 0, 2, 2, 0),
  y = c(0, 0, 2, 2, 0, 0, 2, 2)
)
transformed_polygons <- transform_histology_coordinates(
  polygon_vertices,
  registration_section,
  x_column = "x",
  y_column = "y",
  section_column = "section_id"
)
selected_section <- select_rois(
  section_grid,
  "manual",
  "polygon",
  polygon_vertices = transformed_polygons,
  section_column = "section_id"
)
stopifnot(
  nrow(selected_section$selected_pixels) > 0,
  all(selected_section$selected_pixels$section_id %in% c("A", "B")),
  all(c("polyA", "polyB") %in% unique(selected_section$selected_pixels$roi_id))
)
selected_samples <- sample_subregions(
  selected_section$annotated_pixels,
  grid_size = 2,
  min_pixels = 1,
  roi_column = "roi_id",
  grid_scope = "roi"
)
stopifnot(
  nrow(selected_samples$sample_matrix) > 0,
  all(c("polyA", "polyB") %in% unique(selected_samples$sample_mapping$roi_id))
)
