if (nzchar(Sys.getenv("_R_CHECK_PACKAGE_NAME_"))) {
  library(SpatialOmicsMSI)
} else {
  source(file.path(getwd(), "R", "msi_pipeline.R"))
  source(file.path(getwd(), "R", "distance_gradient.R"))
}

grid <- expand.grid(x = 1:9, y = 1:9)
one_section <- transform(
  grid,
  pixel_id = seq_len(nrow(grid)),
  section_id = "section_1",
  subject_id = "subject_1",
  roi_id = ifelse(x >= 3 & x <= 7 & y >= 3 & y <= 7, "target", "other")
)
one_section$mz_100 <- abs(one_section$x - 5) + abs(one_section$y - 5)

distance <- compute_reference_distance(
  one_section, domain_column = "roi_id", target_domain = "target",
  x_resolution = 50, y_resolution = 25,
  topology_x_step = 1, topology_y_step = 1, distance_unit = "um",
  reference = "boundary", boundary_type = "domain_interface",
  section_column = "section_id", neighbor_method = "rook"
)
stopifnot(inherits(distance, "reference_distance"))
stopifnot(distance$matrix$distance[one_section$x == 5 & one_section$y == 5] == -50)
stopifnot(all(distance$matrix$distance[one_section$roi_id == "target"] <= 0))
stopifnot(all(distance$matrix$distance[one_section$roi_id == "other"] >= 0))
stopifnot(all(is.finite(distance$matrix$reference_x)))

rings <- bin_distance_to_rings(
  distance, method = "fixed_breaks", breaks = c(0, 25, 50, 100, 200),
  separate_sides = TRUE
)
stopifnot(inherits(rings, "distance_rings"))
stopifnot(!any(rings$labels == "ring_NA", na.rm = TRUE))
stopifnot(all(rings$labels[distance$matrix$distance == 0] == "boundary"))
stopifnot(any(grepl("^inner_", rings$labels, useBytes = TRUE), na.rm = TRUE))
stopifnot(any(grepl("^outer_", rings$labels, useBytes = TRUE), na.rm = TRUE))

# Duplicate coordinate grids are valid across sections because topology is built per section.
multi <- do.call(rbind, lapply(1:5, function(subject) {
  out <- one_section
  out$pixel_id <- paste0("s", subject, "_", seq_len(nrow(out)))
  out$section_id <- paste0("section_", subject)
  out$subject_id <- paste0("subject_", subject)
  out$mz_100 <- out$mz_100 + subject / 10
  out
}))
multi_distance <- compute_reference_distance(
  multi, domain_column = "roi_id", target_domain = "target",
  x_resolution = 1, y_resolution = 1,
  topology_x_step = 1, topology_y_step = 1, distance_unit = "grid_step",
  reference = "boundary", section_column = "section_id", neighbor_method = "rook"
)
multi_rings <- bin_distance_to_rings(
  multi_distance, method = "fixed_breaks", breaks = c(0, 1, 2, 3, 5)
)
profiles <- aggregate_distance_profiles(
  multi, multi_distance, multi_rings,
  subject_column = "subject_id", section_column = "section_id",
  features = "mz_100", min_pixels_per_subject_ring = 1,
  min_subjects_per_ring = 5
)
stopifnot(all(profiles$ring_qc$included))
stopifnot(all(profiles$included_profiles$n_pixels >= 1))

if (requireNamespace("mgcv", quietly = TRUE)) {
  exploratory <- fit_distance_gam(
    multi$mz_100, multi_distance$matrix$distance,
    analysis_mode = "exploratory_pixel", k = 5
  )
  stopifnot(is.na(exploratory$p_value), exploratory$exploratory_only)

  population <- fit_distance_gam(
    profiles$included_profiles$mz_100,
    profiles$included_profiles$distance_midpoint,
    analysis_mode = "population_subject",
    subject_id = profiles$included_profiles$subject_id,
    section_id = profiles$included_profiles$section_id,
    data_are_aggregated = TRUE, k = 4, min_subjects = 5
  )
  stopifnot(!population$exploratory_only, population$n_subjects == 5L)
}

# A disconnected target fails by default and can be restricted explicitly.
disconnected <- one_section
disconnected$roi_id <- "other"
disconnected$roi_id[disconnected$x == 2 & disconnected$y == 2] <- "target"
disconnected$roi_id[disconnected$x == 8 & disconnected$y == 8] <- "target"
component_failed <- tryCatch({
  compute_reference_distance(
    disconnected, "roi_id", "target",
    x_resolution = 1, y_resolution = 1,
    topology_x_step = 1, topology_y_step = 1, distance_unit = "pixel",
    section_column = "section_id", component_action = "error"
  )
  FALSE
}, error = function(e) TRUE)
stopifnot(component_failed)

largest <- compute_reference_distance(
  disconnected, "roi_id", "target",
  x_resolution = 1, y_resolution = 1,
  topology_x_step = 1, topology_y_step = 1, distance_unit = "pixel",
  section_column = "section_id", component_action = "largest",
  boundary_type = "domain_interface"
)
stopifnot(sum(largest$matrix$inside_target, na.rm = TRUE) == 1L)

if (requireNamespace("mgcv", quietly = TRUE)) {
  integrated <- analyze_distance_gradient(
    multi, domain_column = "roi_id", target_domain = "target",
    ring_method = "fixed_breaks",
    x_resolution = 1, y_resolution = 1,
    topology_x_step = 1, topology_y_step = 1,
    distance_unit = "grid_step", domain_source_features = character(0),
    test_features = "mz_100", section_column = "section_id",
    subject_column = "subject_id", ring_breaks = c(0, 1, 2, 3, 5),
    analysis_mode = "population_subject", gam_k = 4,
    min_subjects = 5, min_pixels_per_subject_ring = 1,
    min_subjects_per_ring = 5
  )
  stopifnot(nrow(integrated$continuous_result) == 1L)
  stopifnot(identical(integrated$continuous_result$status, "fitted"))
  stopifnot(!is.null(integrated$discrete_result))
}
