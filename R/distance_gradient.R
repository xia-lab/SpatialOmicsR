# Distance-to-domain gradient analysis -------------------------------------

find_connected_components <- function(node_indices, edges) {
  node_indices <- as.integer(node_indices)
  if (!length(node_indices)) return(integer())
  position <- stats::setNames(seq_along(node_indices), as.character(node_indices))
  adjacency <- vector("list", length(node_indices))
  if (nrow(edges)) {
    for (i in seq_len(nrow(edges))) {
      a <- unname(position[as.character(edges$from[i])])
      b <- unname(position[as.character(edges$to[i])])
      if (!is.na(a) && !is.na(b)) {
        adjacency[[a]] <- c(adjacency[[a]], b)
        adjacency[[b]] <- c(adjacency[[b]], a)
      }
    }
  }
  component <- integer(length(node_indices))
  component_number <- 0L
  for (start in seq_along(node_indices)) {
    if (component[start] != 0L) next
    component_number <- component_number + 1L
    queue <- integer(length(node_indices))
    queue[1] <- start
    component[start] <- component_number
    head <- 1L
    tail <- 1L
    while (head <= tail) {
      current <- queue[head]
      head <- head + 1L
      for (neighbor in adjacency[[current]]) {
        if (component[neighbor] == 0L) {
          tail <- tail + 1L
          queue[tail] <- neighbor
          component[neighbor] <- component_number
        }
      }
    }
  }
  component
}

compute_reference_distance <- function(pixel_matrix,
                                       domain_column,
                                       target_domain,
                                       x_resolution,
                                       y_resolution,
                                       topology_x_step,
                                       topology_y_step,
                                       distance_unit,
                                       reference = c("boundary", "centroid"),
                                       boundary_type = c(
                                         "domain_interface",
                                         "observed_domain_perimeter",
                                         "domain_interface_and_perimeter"
                                       ),
                                       x_col = "x",
                                       y_col = "y",
                                       section_column = NULL,
                                       neighbor_method = c("queen", "rook"),
                                       component_action = c("error", "nearest", "largest")) {
  reference <- match.arg(reference)
  boundary_type <- match.arg(boundary_type)
  neighbor_method <- match.arg(neighbor_method)
  component_action <- match.arg(component_action)
  required_columns(pixel_matrix, c(x_col, y_col, domain_column, section_column), "Pixel matrix")
  if (!is.numeric(pixel_matrix[[x_col]]) || !is.numeric(pixel_matrix[[y_col]]) ||
      any(!is.finite(pixel_matrix[[x_col]])) || any(!is.finite(pixel_matrix[[y_col]]))) {
    stop("Spatial coordinates must be finite numeric values.", call. = FALSE)
  }
  positive_scalar <- function(value, name) {
    if (length(value) != 1L || !is.finite(value) || value <= 0) {
      stop(name, " must be supplied as one positive finite value.", call. = FALSE)
    }
    as.numeric(value)
  }
  x_resolution <- positive_scalar(x_resolution, "x_resolution")
  y_resolution <- positive_scalar(y_resolution, "y_resolution")
  topology_x_step <- positive_scalar(topology_x_step, "topology_x_step")
  topology_y_step <- positive_scalar(topology_y_step, "topology_y_step")
  if (length(distance_unit) != 1L || is.na(distance_unit) || !nzchar(as.character(distance_unit))) {
    stop("distance_unit must be supplied explicitly, for example 'um' or 'pixel'.", call. = FALSE)
  }
  if (length(target_domain) != 1L || is.na(target_domain)) {
    stop("target_domain must be one non-missing value.", call. = FALSE)
  }
  section <- if (is.null(section_column)) {
    rep(".__single_section__", nrow(pixel_matrix))
  } else as.character(pixel_matrix[[section_column]])
  if (anyNA(section) || any(!nzchar(section))) {
    stop("Section identifiers must be non-missing and non-empty.", call. = FALSE)
  }
  domain <- as.character(pixel_matrix[[domain_column]])
  target <- as.character(target_domain)
  physical_x <- pixel_matrix[[x_col]] * x_resolution
  physical_y <- pixel_matrix[[y_col]] * y_resolution
  distance <- rep(NA_real_, nrow(pixel_matrix))
  inside <- rep(NA, nrow(pixel_matrix))
  component_id <- rep(NA_integer_, nrow(pixel_matrix))
  nearest_component_id <- rep(NA_integer_, nrow(pixel_matrix))
  reference_x <- rep(NA_real_, nrow(pixel_matrix))
  reference_y <- rep(NA_real_, nrow(pixel_matrix))
  section_qc <- list()
  component_qc <- list()

  for (section_number in seq_along(unique(section))) {
    section_name <- unique(section)[section_number]
    global_index <- which(section == section_name)
    local_x <- pixel_matrix[[x_col]][global_index] / topology_x_step
    local_y <- pixel_matrix[[y_col]][global_index] / topology_y_step
    if (any(abs(local_x - round(local_x)) > sqrt(.Machine$double.eps)) ||
        any(abs(local_y - round(local_y)) > sqrt(.Machine$double.eps))) {
      stop(
        "Coordinates in section '", section_name,
        "' are not aligned to topology_x_step/topology_y_step.",
        call. = FALSE
      )
    }
    graph <- build_spatial_neighbors(round(local_x), round(local_y), method = neighbor_method)
    edges <- graph$undirected_edges
    active_domain <- !is.na(domain[global_index]) & domain[global_index] == target
    original_target_count <- sum(active_domain)
    if (!original_target_count) {
      section_qc[[section_number]] <- data.frame(
        section_id = section_name, n_pixels = length(global_index),
        n_target_pixels_original = 0L, n_target_pixels_used = 0L,
        n_components_original = 0L, target_present = FALSE,
        boundary_pixels = 0L, status = "target_absent", stringsAsFactors = FALSE
      )
      next
    }
    original_domain_index <- which(active_domain)
    domain_edges <- edges[
      edges$from %in% original_domain_index & edges$to %in% original_domain_index,
      , drop = FALSE
    ]
    original_components <- find_connected_components(original_domain_index, domain_edges)
    n_components <- max(original_components)
    original_component_sizes <- table(original_components)
    if (n_components > 1L && component_action == "error") {
      stop(
        "target_domain '", target, "' in section '", section_name, "' has ",
        n_components, " disconnected components. Choose 'nearest' or 'largest' explicitly.",
        call. = FALSE
      )
    }
    if (n_components > 1L && component_action == "largest") {
      largest <- as.integer(names(original_component_sizes)[which.max(original_component_sizes)])
      keep_local <- original_domain_index[original_components == largest]
      active_domain <- seq_along(active_domain) %in% keep_local
    }
    domain_index <- which(active_domain)
    active_edges <- edges[edges$from %in% domain_index & edges$to %in% domain_index, , drop = FALSE]
    active_components <- find_connected_components(domain_index, active_edges)
    global_domain_index <- global_index[domain_index]
    component_id[global_domain_index] <- active_components
    component_qc[[section_number]] <- data.frame(
      section_id = section_name,
      component_id = seq_len(max(original_components)),
      n_pixels = as.integer(original_component_sizes[as.character(seq_len(max(original_components)))]),
      retained = if (component_action == "largest" && n_components > 1L) {
        seq_len(max(original_components)) == as.integer(names(original_component_sizes)[which.max(original_component_sizes)])
      } else TRUE,
      stringsAsFactors = FALSE
    )
    local_px <- physical_x[global_index]
    local_py <- physical_y[global_index]

    if (reference == "centroid") {
      active_component_labels <- sort(unique(active_components))
      centroids <- do.call(rbind, lapply(active_component_labels, function(component) {
        index <- domain_index[active_components == component]
        data.frame(component_id = component, x = mean(local_px[index]), y = mean(local_py[index]))
      }))
      distance_to_centroid <- vapply(seq_len(nrow(centroids)), function(i) {
        sqrt((local_px - centroids$x[i])^2 + (local_py - centroids$y[i])^2)
      }, numeric(length(global_index)))
      if (is.null(dim(distance_to_centroid))) distance_to_centroid <- matrix(distance_to_centroid, ncol = 1L)
      selected <- max.col(-distance_to_centroid, ties.method = "first")
      distance[global_index] <- distance_to_centroid[cbind(seq_along(global_index), selected)]
      nearest_component_id[global_index] <- centroids$component_id[selected]
      reference_x[global_index] <- centroids$x[selected]
      reference_y[global_index] <- centroids$y[selected]
      inside[global_index] <- active_domain
      section_qc[[section_number]] <- data.frame(
        section_id = section_name, n_pixels = length(global_index),
        n_target_pixels_original = original_target_count,
        n_target_pixels_used = sum(active_domain),
        n_components_original = n_components, target_present = TRUE,
        boundary_pixels = NA_integer_, status = "computed", stringsAsFactors = FALSE
      )
      next
    }

    expected_degree <- if (neighbor_method == "rook") 4L else 8L
    boundary_flag <- logical(length(global_index))
    for (i in domain_index) {
      adjacent <- unique(c(edges$to[edges$from == i], edges$from[edges$to == i]))
      interface <- length(adjacent) > 0L && any(!active_domain[adjacent])
      observed_perimeter <- graph$degree[i] < expected_degree
      boundary_flag[i] <- switch(
        boundary_type,
        domain_interface = interface,
        observed_domain_perimeter = observed_perimeter,
        domain_interface_and_perimeter = interface || observed_perimeter
      )
    }
    boundary_index <- which(boundary_flag)
    if (!length(boundary_index)) {
      warning("No requested boundary pixels found in section '", section_name, "'.", call. = FALSE)
      inside[global_index] <- active_domain
      section_qc[[section_number]] <- data.frame(
        section_id = section_name, n_pixels = length(global_index),
        n_target_pixels_original = original_target_count,
        n_target_pixels_used = sum(active_domain),
        n_components_original = n_components, target_present = TRUE,
        boundary_pixels = 0L, status = "boundary_absent", stringsAsFactors = FALSE
      )
      next
    }
    boundary_component <- active_components[match(boundary_index, domain_index)]
    nearest_boundary <- vapply(boundary_index, function(i) {
      sqrt((local_px - local_px[i])^2 + (local_py - local_py[i])^2)
    }, numeric(length(global_index)))
    if (is.null(dim(nearest_boundary))) nearest_boundary <- matrix(nearest_boundary, ncol = 1L)
    selected <- max.col(-nearest_boundary, ties.method = "first")
    local_distance <- nearest_boundary[cbind(seq_along(global_index), selected)]
    local_distance[active_domain] <- -local_distance[active_domain]
    distance[global_index] <- local_distance
    inside[global_index] <- active_domain
    nearest_component_id[global_index] <- boundary_component[selected]
    reference_x[global_index] <- local_px[boundary_index[selected]]
    reference_y[global_index] <- local_py[boundary_index[selected]]
    section_qc[[section_number]] <- data.frame(
      section_id = section_name, n_pixels = length(global_index),
      n_target_pixels_original = original_target_count,
      n_target_pixels_used = sum(active_domain),
      n_components_original = n_components, target_present = TRUE,
      boundary_pixels = length(boundary_index), status = "computed", stringsAsFactors = FALSE
    )
  }
  if (!any(vapply(section_qc, function(x) isTRUE(x$target_present[1]), logical(1)))) {
    stop("target_domain is absent from every section.", call. = FALSE)
  }
  matrix_output <- data.frame(
    pixel_index = seq_len(nrow(pixel_matrix)), distance = distance,
    inside_target = inside, component_id = component_id,
    nearest_component_id = nearest_component_id,
    reference_x = reference_x, reference_y = reference_y,
    reference_type = reference, section_id = section,
    distance_unit = as.character(distance_unit), stringsAsFactors = FALSE
  )
  structure(list(
    matrix = matrix_output,
    section_qc = do.call(rbind, section_qc),
    component_qc = if (length(component_qc)) do.call(rbind, component_qc) else data.frame(),
    settings = list(
      domain_column = domain_column, target_domain = target,
      reference = reference, boundary_type = boundary_type,
      neighbor_method = neighbor_method, component_action = component_action,
      x_resolution = x_resolution, y_resolution = y_resolution,
      topology_x_step = topology_x_step, topology_y_step = topology_y_step,
      distance_unit = as.character(distance_unit)
    ),
    interpretation = if (boundary_type == "observed_domain_perimeter") {
      "Perimeter means the edge of the observed grid/mask, including holes or acquisition gaps; it is not automatically a biological tissue outline."
    } else NULL
  ), class = "reference_distance")
}

bin_distance_to_rings <- function(distance_result,
                                  method,
                                  n_rings = 5,
                                  ring_width = NULL,
                                  breaks = NULL,
                                  separate_sides = TRUE) {
  if (!inherits(distance_result, "reference_distance")) {
    stop("distance_result must come from compute_reference_distance().", call. = FALSE)
  }
  method <- match.arg(method, c("quantile", "fixed_width", "fixed_breaks"))
  if (length(n_rings) != 1L || !is.finite(n_rings) || n_rings < 1 || n_rings != as.integer(n_rings)) {
    stop("n_rings must be one positive integer.", call. = FALSE)
  }
  if (method == "fixed_width" &&
      (is.null(ring_width) || length(ring_width) != 1L || !is.finite(ring_width) || ring_width <= 0)) {
    stop("ring_width must be one positive finite value for fixed_width.", call. = FALSE)
  }
  if (method == "fixed_breaks") {
    if (is.null(breaks) || length(breaks) < 2L || any(!is.finite(breaks)) ||
        any(breaks < 0) || breaks[1] != 0 || is.unsorted(breaks, strictly = TRUE)) {
      stop("breaks must be strictly increasing finite non-negative values beginning at zero.", call. = FALSE)
    }
  }
  distance <- distance_result$matrix$distance
  reference_type <- distance_result$settings$reference
  labels <- rep(NA_character_, length(distance))
  break_rows <- list()
  make_side <- function(index, prefix) {
    values <- abs(distance[index])
    finite <- is.finite(values)
    if (!any(finite)) return(invisible(NULL))
    maximum <- max(values[finite])
    side_breaks <- switch(method,
      quantile = unique(as.numeric(stats::quantile(
        values[finite], probs = seq(0, 1, length.out = as.integer(n_rings) + 1L),
        names = FALSE, na.rm = TRUE
      ))),
      fixed_width = unique(c(seq(0, maximum, by = ring_width), maximum)),
      fixed_breaks = breaks
    )
    if (method != "fixed_breaks" && tail(side_breaks, 1) < maximum) {
      side_breaks <- c(side_breaks, maximum)
    }
    side_breaks <- sort(unique(side_breaks))
    if (length(side_breaks) < 2L) return(invisible(NULL))
    ring <- cut(values, breaks = side_breaks, include.lowest = TRUE, labels = FALSE)
    valid_ring <- !is.na(ring)
    labels[index[valid_ring]] <<- sprintf("%s_%02d", prefix, ring[valid_ring])
    break_rows[[length(break_rows) + 1L]] <<- data.frame(
      side = prefix, ring_index = seq_len(length(side_breaks) - 1L),
      lower = head(side_breaks, -1L), upper = tail(side_breaks, -1L),
      stringsAsFactors = FALSE
    )
    invisible(NULL)
  }
  zero <- is.finite(distance) & distance == 0
  if (reference_type == "boundary" && isTRUE(separate_sides)) {
    make_side(which(is.finite(distance) & distance < 0), "inner")
    make_side(which(is.finite(distance) & distance > 0), "outer")
    labels[zero] <- "boundary"
  } else {
    make_side(which(is.finite(distance)), "ring")
  }
  ring_qc <- data.frame(
    distance_ring = sort(unique(labels[!is.na(labels)])),
    stringsAsFactors = FALSE
  )
  ring_qc$n_pixels <- as.integer(table(factor(labels, levels = ring_qc$distance_ring)))
  unassigned_finite <- sum(is.finite(distance) & is.na(labels))
  structure(list(
    labels = labels,
    breaks = if (length(break_rows)) do.call(rbind, break_rows) else data.frame(),
    qc = ring_qc,
    settings = list(
      method = method, n_rings = as.integer(n_rings), ring_width = ring_width,
      supplied_breaks = breaks, separate_sides = isTRUE(separate_sides),
      distance_unit = distance_result$settings$distance_unit,
      n_finite_unassigned = unassigned_finite
    )
  ), class = "distance_rings")
}

aggregate_distance_profiles <- function(pixel_matrix,
                                        distance_result,
                                        rings,
                                        subject_column,
                                        section_column = NULL,
                                        features = feature_columns(pixel_matrix),
                                        summary_method = c("mean", "median"),
                                        min_pixels_per_subject_ring = 10,
                                        min_subjects_per_ring = 2) {
  if (!inherits(distance_result, "reference_distance") || !inherits(rings, "distance_rings")) {
    stop("distance_result and rings must come from the distance-gradient helpers.", call. = FALSE)
  }
  summary_method <- match.arg(summary_method)
  required_columns(pixel_matrix, c(subject_column, section_column, features), "Pixel matrix")
  if (length(features) == 0L || any(!vapply(pixel_matrix[features], is.numeric, logical(1)))) {
    stop("features must name at least one numeric pixel-matrix column.", call. = FALSE)
  }
  if (nrow(pixel_matrix) != nrow(distance_result$matrix) || length(rings$labels) != nrow(pixel_matrix)) {
    stop("Distance, ring, and pixel rows must align one-to-one.", call. = FALSE)
  }
  if (length(min_pixels_per_subject_ring) != 1L || !is.finite(min_pixels_per_subject_ring) ||
      min_pixels_per_subject_ring < 1 || min_pixels_per_subject_ring != as.integer(min_pixels_per_subject_ring) ||
      length(min_subjects_per_ring) != 1L || !is.finite(min_subjects_per_ring) ||
      min_subjects_per_ring < 2 || min_subjects_per_ring != as.integer(min_subjects_per_ring)) {
    stop("Minimum pixel and subject counts must be positive integers; min_subjects_per_ring must be at least two.", call. = FALSE)
  }
  subject <- as.character(pixel_matrix[[subject_column]])
  section <- if (is.null(section_column)) distance_result$matrix$section_id else as.character(pixel_matrix[[section_column]])
  if (anyNA(subject) || any(!nzchar(subject)) || anyNA(section) || any(!nzchar(section))) {
    stop("Subject and section identifiers must be non-missing and non-empty.", call. = FALSE)
  }
  work <- data.frame(
    subject_id = subject, section_id = section,
    distance_ring = rings$labels,
    distance = distance_result$matrix$distance,
    pixel_matrix[, features, drop = FALSE],
    check.names = FALSE, stringsAsFactors = FALSE
  )
  work <- work[!is.na(work$distance_ring) & is.finite(work$distance), , drop = FALSE]
  if (!nrow(work)) stop("No finite ring-assigned pixels are available for aggregation.", call. = FALSE)
  group_key <- interaction(work$subject_id, work$section_id, work$distance_ring, drop = TRUE, lex.order = TRUE)
  split_rows <- split(seq_len(nrow(work)), group_key)
  summarize <- if (summary_method == "mean") {
    function(x) mean(x, na.rm = TRUE)
  } else function(x) stats::median(x, na.rm = TRUE)
  profiles <- do.call(rbind, lapply(split_rows, function(index) {
    values <- vapply(features, function(feature) {
      finite <- is.finite(work[[feature]][index])
      if (!any(finite)) NA_real_ else summarize(work[[feature]][index][finite])
    }, numeric(1))
    data.frame(
      sample_id = paste(work$subject_id[index[1]], work$section_id[index[1]], work$distance_ring[index[1]], sep = "::"),
      subject_id = work$subject_id[index[1]], section_id = work$section_id[index[1]],
      distance_ring = work$distance_ring[index[1]],
      distance_midpoint = mean(work$distance[index], na.rm = TRUE),
      n_pixels = length(index), t(values), check.names = FALSE,
      stringsAsFactors = FALSE
    )
  }))
  rownames(profiles) <- NULL
  # Eligibility is based on total pixels per subject and ring across sections,
  # while profiles remain subject x section x ring for model auditing.
  subject_ring <- stats::aggregate(
    profiles$n_pixels,
    by = list(subject_id = profiles$subject_id, distance_ring = profiles$distance_ring),
    FUN = sum
  )
  names(subject_ring)[3] <- "n_pixels"
  ring_names <- sort(unique(profiles$distance_ring))
  ring_qc <- do.call(rbind, lapply(ring_names, function(ring) {
    rows <- subject_ring$distance_ring == ring
    passing <- subject_ring$n_pixels[rows] >= min_pixels_per_subject_ring
    data.frame(
      distance_ring = ring,
      n_subjects_observed = sum(rows),
      n_subjects_passing_pixels = sum(passing),
      total_pixels = sum(subject_ring$n_pixels[rows]),
      included = sum(passing) >= min_subjects_per_ring,
      exclusion_reason = if (sum(passing) >= min_subjects_per_ring) NA_character_ else "too_few_subjects_with_minimum_pixels",
      stringsAsFactors = FALSE
    )
  }))
  eligible_rings <- ring_qc$distance_ring[ring_qc$included]
  eligible_subject_ring <- subject_ring[
    subject_ring$distance_ring %in% eligible_rings &
      subject_ring$n_pixels >= min_pixels_per_subject_ring,
    c("subject_id", "distance_ring"), drop = FALSE
  ]
  eligibility_key <- paste(eligible_subject_ring$subject_id, eligible_subject_ring$distance_ring, sep = "\r")
  profiles$included <- paste(profiles$subject_id, profiles$distance_ring, sep = "\r") %in% eligibility_key
  list(
    profiles = profiles,
    included_profiles = profiles[profiles$included, , drop = FALSE],
    subject_ring_qc = subject_ring,
    ring_qc = ring_qc,
    settings = list(
      subject_column = subject_column, section_column = section_column,
      summary_method = summary_method,
      min_pixels_per_subject_ring = as.integer(min_pixels_per_subject_ring),
      min_subjects_per_ring = as.integer(min_subjects_per_ring)
    )
  )
}

fit_distance_gam <- function(intensity,
                             distance,
                             analysis_mode = c("exploratory_pixel", "population_subject"),
                             subject_id = NULL,
                             section_id = NULL,
                             data_are_aggregated = FALSE,
                             k = 10,
                             subject_k = NULL,
                             min_subjects = 5,
                             family = stats::gaussian()) {
  analysis_mode <- match.arg(analysis_mode)
  if (!requireNamespace("mgcv", quietly = TRUE)) {
    stop("Package 'mgcv' is required for distance-based GAM fitting.", call. = FALSE)
  }
  if (length(k) != 1L || !is.finite(k) || k < 3 || k != as.integer(k)) {
    stop("k must be one integer of at least three.", call. = FALSE)
  }
  work <- data.frame(intensity = as.numeric(intensity), distance = as.numeric(distance))
  if (analysis_mode == "population_subject") {
    if (!isTRUE(data_are_aggregated)) {
      stop("population_subject mode only accepts subject-level aggregated distance profiles.", call. = FALSE)
    }
    if (is.null(subject_id) || length(subject_id) != nrow(work)) {
      stop("subject_id must contain one identifier per aggregated profile.", call. = FALSE)
    }
    work$subject <- factor(subject_id)
    if (!is.null(section_id)) {
      if (length(section_id) != nrow(work)) stop("section_id must align with aggregated profiles.", call. = FALSE)
      work$subject_section <- factor(interaction(subject_id, section_id, drop = TRUE, sep = "::"))
    }
  }
  finite <- is.finite(work$intensity) & is.finite(work$distance)
  if (analysis_mode == "population_subject") finite <- finite & !is.na(work$subject)
  work <- droplevels(work[finite, , drop = FALSE])
  unique_distance <- length(unique(work$distance))
  k_used <- min(as.integer(k), unique_distance - 1L)
  if (nrow(work) < 2L * max(3L, k_used) || k_used < 3L) {
    stop("Insufficient finite observations or distinct distance values for the requested GAM.", call. = FALSE)
  }
  if (analysis_mode == "exploratory_pixel") {
    fit <- mgcv::gam(intensity ~ s(distance, k = k_used), data = work, method = "REML", family = family)
    summary_fit <- summary(fit)
    return(list(
      fit = fit, mode = analysis_mode,
      edf = summary_fit$s.table[1, "edf"], p_value = NA_real_,
      exploratory_only = TRUE, r_squared = summary_fit$r.sq,
      deviance_explained = summary_fit$dev.expl, n = nrow(work),
      n_subjects = NA_integer_, k_used = k_used,
      interpretation = "Pixel-level curve is descriptive; no inferential p-value is returned."
    ))
  }
  if (length(min_subjects) != 1L || !is.finite(min_subjects) || min_subjects < 3 ||
      min_subjects != as.integer(min_subjects)) {
    stop("min_subjects must be an integer of at least three.", call. = FALSE)
  }
  n_subjects <- nlevels(work$subject)
  if (n_subjects < min_subjects) {
    stop("Population GAM requires at least ", min_subjects, " subjects; found ", n_subjects, ".", call. = FALSE)
  }
  formula_terms <- c(sprintf("s(distance, k = %d)", k_used), "s(subject, bs = 're')")
  if (!is.null(section_id) && nlevels(work$subject_section) > n_subjects) {
    formula_terms <- c(formula_terms, "s(subject_section, bs = 're')")
  }
  if (!is.null(subject_k)) {
    if (length(subject_k) != 1L || !is.finite(subject_k) || subject_k < 3 || subject_k != as.integer(subject_k)) {
      stop("subject_k must be NULL or one integer of at least three.", call. = FALSE)
    }
    formula_terms <- c(formula_terms, sprintf("s(distance, subject, bs = 'fs', k = %d)", as.integer(subject_k)))
  }
  formula <- stats::as.formula(paste("intensity ~", paste(formula_terms, collapse = " + ")))
  fit <- mgcv::gam(formula, data = work, method = "REML", family = family)
  summary_fit <- summary(fit)
  list(
    fit = fit, mode = analysis_mode,
    edf = summary_fit$s.table[1, "edf"],
    p_value = summary_fit$s.table[1, "p-value"],
    exploratory_only = FALSE, r_squared = summary_fit$r.sq,
    deviance_explained = summary_fit$dev.expl, n = nrow(work),
    n_subjects = n_subjects, k_used = k_used,
    interpretation = paste(
      "Inference uses aggregated subject-by-section distance profiles.",
      "The smooth-term p-value remains an mgcv approximation and requires model diagnostics."
    )
  )
}

compute_distance_variable_metabolites <- function(pixel_matrix,
                                                   distance_result,
                                                   features = feature_columns(pixel_matrix),
                                                   analysis_mode = c("exploratory_pixel", "population_subject"),
                                                   aggregated_profiles = NULL,
                                                   k = 10,
                                                   subject_k = NULL,
                                                   min_subjects = 5,
                                                   p_adjust_method = "BH") {
  analysis_mode <- match.arg(analysis_mode)
  if (!inherits(distance_result, "reference_distance")) {
    stop("distance_result must come from compute_reference_distance().", call. = FALSE)
  }
  required_columns(pixel_matrix, features, "Pixel matrix")
  if (!length(features) || any(!vapply(pixel_matrix[features], is.numeric, logical(1)))) {
    stop("features must name at least one numeric column.", call. = FALSE)
  }
  if (analysis_mode == "population_subject") {
    if (is.null(aggregated_profiles) || !is.data.frame(aggregated_profiles)) {
      stop("aggregated_profiles is required for population_subject mode.", call. = FALSE)
    }
    required_columns(
      aggregated_profiles,
      c("subject_id", "section_id", "distance_midpoint", features),
      "Aggregated profiles"
    )
  }
  fits <- vector("list", length(features))
  rows <- lapply(seq_along(features), function(index) {
    feature <- features[index]
    calculation <- tryCatch({
      if (analysis_mode == "exploratory_pixel") {
        fit_distance_gam(
          pixel_matrix[[feature]], distance_result$matrix$distance,
          analysis_mode = analysis_mode, k = k, min_subjects = min_subjects
        )
      } else {
        fit_distance_gam(
          aggregated_profiles[[feature]], aggregated_profiles$distance_midpoint,
          analysis_mode = analysis_mode,
          subject_id = aggregated_profiles$subject_id,
          section_id = aggregated_profiles$section_id,
          data_are_aggregated = TRUE, k = k, subject_k = subject_k,
          min_subjects = min_subjects
        )
      }
    }, error = function(error) error)
    if (inherits(calculation, "error")) {
      return(data.frame(
        feature = feature, edf = NA_real_, p_value = NA_real_,
        r_squared = NA_real_, deviance_explained = NA_real_,
        n = NA_integer_, n_subjects = NA_integer_, k_used = NA_integer_,
        exploratory_only = analysis_mode == "exploratory_pixel",
        status = "failed", error_message = conditionMessage(calculation),
        stringsAsFactors = FALSE
      ))
    }
    fits[[index]] <<- calculation$fit
    data.frame(
      feature = feature, edf = calculation$edf, p_value = calculation$p_value,
      r_squared = calculation$r_squared,
      deviance_explained = calculation$deviance_explained,
      n = calculation$n, n_subjects = calculation$n_subjects,
      k_used = calculation$k_used,
      exploratory_only = calculation$exploratory_only,
      status = "fitted", error_message = NA_character_, stringsAsFactors = FALSE
    )
  })
  result <- do.call(rbind, rows)
  result$adj_p_value <- if (analysis_mode == "exploratory_pixel") {
    NA_real_
  } else stats::p.adjust(result$p_value, method = p_adjust_method)
  result <- result[order(
    if (analysis_mode == "exploratory_pixel") -result$deviance_explained else result$adj_p_value,
    na.last = TRUE
  ), , drop = FALSE]
  rownames(result) <- NULL
  names(fits) <- features
  list(
    table = result,
    fits = fits,
    settings = list(
      analysis_mode = analysis_mode, k = k, subject_k = subject_k,
      min_subjects = min_subjects, p_adjust_method = p_adjust_method
    )
  )
}

analyze_distance_gradient <- function(pixel_matrix,
                                      domain_column,
                                      target_domain,
                                      ring_method,
                                      x_resolution,
                                      y_resolution,
                                      topology_x_step,
                                      topology_y_step,
                                      distance_unit,
                                      domain_source_features = NULL,
                                      test_features = NULL,
                                      reference = c("boundary", "centroid"),
                                      boundary_type = c(
                                        "domain_interface",
                                        "observed_domain_perimeter",
                                        "domain_interface_and_perimeter"
                                      ),
                                      x_col = "x",
                                      y_col = "y",
                                      section_column = NULL,
                                      neighbor_method = c("queen", "rook"),
                                      component_action = c("error", "nearest", "largest"),
                                      n_rings = 5,
                                      ring_width = NULL,
                                      ring_breaks = NULL,
                                      separate_sides = TRUE,
                                      analysis_mode = c("exploratory_pixel", "population_subject"),
                                      subject_column = NULL,
                                      gam_k = 10,
                                      subject_k = NULL,
                                      min_subjects = 5,
                                      profile_summary = c("mean", "median"),
                                      min_pixels_per_subject_ring = 10,
                                      min_subjects_per_ring = 2,
                                      allow_circular_analysis = FALSE,
                                      p_adjust_method = "BH") {
  reference <- match.arg(reference)
  boundary_type <- match.arg(boundary_type)
  neighbor_method <- match.arg(neighbor_method)
  component_action <- match.arg(component_action)
  ring_method <- match.arg(ring_method, c("quantile", "fixed_width", "fixed_breaks"))
  analysis_mode <- match.arg(analysis_mode)
  profile_summary <- match.arg(profile_summary)
  all_features <- feature_columns(pixel_matrix)
  test_features <- if (is.null(test_features)) all_features else as.character(test_features)
  missing_features <- setdiff(test_features, names(pixel_matrix))
  non_mz_features <- setdiff(test_features, all_features)
  if (!length(test_features) || length(missing_features) || length(non_mz_features) ||
      any(!vapply(pixel_matrix[test_features], is.numeric, logical(1)))) {
    stop(
      "test_features must name existing numeric mz_ columns",
      if (length(missing_features)) paste0(": ", paste(missing_features, collapse = ", ")) else ".",
      call. = FALSE
    )
  }
  if (is.null(domain_source_features) && !isTRUE(allow_circular_analysis)) {
    stop(
      "domain_source_features must be declared. Use character(0) for an independent modality, ",
      "or allow_circular_analysis=TRUE for explicitly descriptive analysis.",
      call. = FALSE
    )
  }
  overlap <- intersect(as.character(domain_source_features), test_features)
  if (length(overlap) && !isTRUE(allow_circular_analysis)) {
    stop(
      "Circular analysis detected for feature(s): ", paste(overlap, collapse = ", "),
      ". Remove them from test_features or explicitly allow descriptive circular analysis.",
      call. = FALSE
    )
  }
  if (analysis_mode == "population_subject" && is.null(subject_column)) {
    stop("subject_column is required for population_subject mode.", call. = FALSE)
  }
  required_columns(pixel_matrix, c(domain_column, subject_column, section_column), "Pixel matrix")
  distance_result <- compute_reference_distance(
    pixel_matrix = pixel_matrix, domain_column = domain_column,
    target_domain = target_domain, x_resolution = x_resolution,
    y_resolution = y_resolution, topology_x_step = topology_x_step,
    topology_y_step = topology_y_step, distance_unit = distance_unit,
    reference = reference, boundary_type = boundary_type,
    x_col = x_col, y_col = y_col, section_column = section_column,
    neighbor_method = neighbor_method, component_action = component_action
  )
  rings <- bin_distance_to_rings(
    distance_result, method = ring_method, n_rings = n_rings,
    ring_width = ring_width, breaks = ring_breaks,
    separate_sides = separate_sides
  )
  aggregation <- NULL
  if (!is.null(subject_column)) {
    aggregation <- aggregate_distance_profiles(
      pixel_matrix, distance_result, rings,
      subject_column = subject_column, section_column = section_column,
      features = test_features, summary_method = profile_summary,
      min_pixels_per_subject_ring = min_pixels_per_subject_ring,
      min_subjects_per_ring = min_subjects_per_ring
    )
  }
  continuous <- compute_distance_variable_metabolites(
    pixel_matrix, distance_result, features = test_features,
    analysis_mode = analysis_mode,
    aggregated_profiles = if (is.null(aggregation)) NULL else aggregation$included_profiles,
    k = gam_k, subject_k = subject_k, min_subjects = min_subjects,
    p_adjust_method = p_adjust_method
  )
  discrete <- NULL
  if (!is.null(aggregation)) {
    included <- aggregation$included_profiles
    if (length(unique(included$distance_ring)) >= 2L) {
      discrete <- differential_region_analysis(
        included, group_column = "distance_ring",
        subject_column = "subject_id", section_column = "section_id",
        p_adjust_method = p_adjust_method
      )
    }
  }
  annotated <- pixel_matrix
  annotated$reference_distance <- distance_result$matrix$distance
  annotated$inside_target <- distance_result$matrix$inside_target
  annotated$distance_ring <- rings$labels
  annotated$nearest_reference_component <- distance_result$matrix$nearest_component_id
  annotated$reference_x <- distance_result$matrix$reference_x
  annotated$reference_y <- distance_result$matrix$reference_y
  list(
    distance_result = distance_result,
    rings = rings,
    aggregation = aggregation,
    continuous_result = continuous$table,
    continuous_fits = continuous$fits,
    discrete_result = discrete,
    annotated_pixel_matrix = annotated,
    provenance = list(
      domain_column = domain_column,
      domain_source_features = domain_source_features,
      test_features = test_features,
      circular_overlap = overlap,
      allow_circular_analysis = isTRUE(allow_circular_analysis),
      analysis_mode = analysis_mode
    ),
    interpretation = if (analysis_mode == "exploratory_pixel") {
      "Continuous curves are pixel-level descriptive summaries with no inferential p-values."
    } else {
      "Population inference uses aggregated subject-by-section ring profiles; inspect GAM diagnostics and subject coverage."
    }
  )
}
