# LC-MS/CCS and LCM multimodal integration -------------------------------

validate_ccs_evidence <- function(
    observed_table,
    reference_table,
    candidate_column = "candidate_id",
    observed_ccs_column = "observed_ccs",
    reference_ccs_column = "reference_ccs",
    observed_source_column = "observed_source",
    reference_source_column = "reference_source",
    observation_id_column = NULL,
    metadata_columns = c("adduct", "charge", "ion_mode"),
    ccs_tolerance_pct) {
  if (!is.data.frame(observed_table) || !is.data.frame(reference_table)) {
    stop("observed_table and reference_table must be data frames.", call. = FALSE)
  }
  required_columns(
    observed_table,
    c(candidate_column, observed_ccs_column, observed_source_column),
    "Observed CCS table"
  )
  if (!is.null(observation_id_column)) {
    required_columns(observed_table, observation_id_column, "Observed CCS table")
  }
  required_columns(
    reference_table,
    c(candidate_column, reference_ccs_column, reference_source_column),
    "Reference CCS table"
  )
  if (missing(ccs_tolerance_pct) || !is.numeric(ccs_tolerance_pct) ||
      length(ccs_tolerance_pct) != 1L || !is.finite(ccs_tolerance_pct) ||
      ccs_tolerance_pct <= 0) {
    stop("ccs_tolerance_pct must be one explicitly supplied positive finite number.", call. = FALSE)
  }
  observed_source <- as.character(observed_table[[observed_source_column]])
  allowed_observed <- c("msi_empirical", "lcms_empirical")
  if (any(is.na(observed_source) | !observed_source %in% allowed_observed)) {
    stop("observed_source must be 'msi_empirical' or 'lcms_empirical'.", call. = FALSE)
  }
  reference_source <- as.character(reference_table[[reference_source_column]])
  allowed_reference <- c("measured_library", "ml_predicted")
  if (any(is.na(reference_source) | !reference_source %in% allowed_reference)) {
    stop("reference_source must be 'measured_library' or 'ml_predicted'.", call. = FALSE)
  }

  observed <- observed_table
  reference <- reference_table
  observed$.__observed_row__ <- seq_len(nrow(observed))
  reference$.__reference_row__ <- seq_len(nrow(reference))
  common_metadata <- metadata_columns[
    metadata_columns %in% names(observed) & metadata_columns %in% names(reference)
  ]
  by_columns <- c(candidate_column, common_metadata)
  evidence <- merge(observed, reference, by = by_columns, all.x = TRUE, sort = FALSE)
  evidence <- evidence[order(evidence$.__observed_row__, evidence$.__reference_row__), , drop = FALSE]

  observed_name <- if (observed_ccs_column == reference_ccs_column) {
    paste0(observed_ccs_column, ".x")
  } else observed_ccs_column
  reference_name <- if (observed_ccs_column == reference_ccs_column) {
    paste0(reference_ccs_column, ".y")
  } else reference_ccs_column
  observed_source_name <- if (observed_source_column == reference_source_column) {
    paste0(observed_source_column, ".x")
  } else observed_source_column
  reference_source_name <- if (observed_source_column == reference_source_column) {
    paste0(reference_source_column, ".y")
  } else reference_source_column
  observed_ccs <- suppressWarnings(as.numeric(evidence[[observed_name]]))
  reference_ccs <- suppressWarnings(as.numeric(evidence[[reference_name]]))
  valid_observed <- is.finite(observed_ccs) & observed_ccs > 0
  valid_reference <- is.finite(reference_ccs) & reference_ccs > 0
  if (any(!valid_observed)) {
    stop("Observed CCS values must be finite and positive.", call. = FALSE)
  }
  relative_error <- rep(NA_real_, nrow(evidence))
  relative_error[valid_reference] <-
    abs(observed_ccs[valid_reference] - reference_ccs[valid_reference]) /
    reference_ccs[valid_reference] * 100
  within <- valid_reference & relative_error <= ccs_tolerance_pct + sqrt(.Machine$double.eps)
  scope <- ifelse(
    as.character(evidence[[observed_source_name]]) == "msi_empirical",
    "msi_to_reference_support",
    "lcms_candidate_characterization_only"
  )
  reference_kind <- as.character(evidence[[reference_source_name]])
  grade <- ifelse(
    !valid_reference, "ccs_not_available",
    ifelse(
      !within, "ccs_inconsistent",
      ifelse(reference_kind == "measured_library",
             "ccs_consistent_measured_reference",
             "ccs_consistent_predicted_reference")
    )
  )
  evidence$relative_error_pct <- relative_error
  evidence$within_tolerance <- within
  evidence$ccs_tolerance_pct <- ccs_tolerance_pct
  evidence$evidence_scope <- scope
  evidence$evidence_grade <- grade
  evidence$msi_identity_confirmed_by_ccs <- within & scope == "msi_to_reference_support"
  evidence$.__reference_row__ <- NULL

  pass_count <- stats::ave(as.integer(within), evidence$.__observed_row__, FUN = sum)
  evidence$passing_reference_count <- as.integer(pass_count)
  ambiguity_group <- if (is.null(observation_id_column)) {
    evidence$.__observed_row__
  } else {
    as.character(evidence[[observation_id_column]])
  }
  passing_candidates <- tapply(
    ifelse(within, as.character(evidence[[candidate_column]]), NA_character_),
    ambiguity_group,
    function(value) length(unique(value[!is.na(value)]))
  )
  evidence$passing_candidate_count <- as.integer(passing_candidates[as.character(ambiguity_group)])
  evidence$ambiguous_candidate_set <- evidence$passing_candidate_count > 1L
  evidence$.__observed_row__ <- NULL
  rownames(evidence) <- NULL
  attr(evidence, "interpretation") <- paste(
    "LC-MS-only CCS characterizes LC-MS candidates but cannot confirm which",
    "identity produced an MSI signal when MSI has no ion-mobility measurement."
  )
  evidence
}

roi_pixels_to_polygon <- function(x, y, x_resolution, y_resolution = x_resolution) {
  if (!is.numeric(x) || !is.numeric(y) || length(x) != length(y) || !length(x) ||
      any(!is.finite(x)) || any(!is.finite(y))) {
    stop("x and y must be equal-length, non-empty finite numeric vectors.", call. = FALSE)
  }
  resolution <- c(x_resolution, y_resolution)
  if (any(!is.finite(resolution) | resolution <= 0)) {
    stop("x_resolution and y_resolution must be positive finite numbers.", call. = FALSE)
  }
  centers <- unique(data.frame(x = x, y = y))
  tol <- sqrt(.Machine$double.eps) * max(1, abs(c(centers$x, centers$y)), resolution)
  aligned_x <- abs((centers$x - min(centers$x)) / x_resolution -
                     round((centers$x - min(centers$x)) / x_resolution)) <= tol
  aligned_y <- abs((centers$y - min(centers$y)) / y_resolution -
                     round((centers$y - min(centers$y)) / y_resolution)) <= tol
  if (!all(aligned_x) || !all(aligned_y)) {
    stop("Pixel centers are not aligned to the supplied regular-grid resolutions.", call. = FALSE)
  }

  half_x <- x_resolution / 2
  half_y <- y_resolution / 2
  edges <- do.call(rbind, lapply(seq_len(nrow(centers)), function(i) {
    left <- centers$x[i] - half_x; right <- centers$x[i] + half_x
    bottom <- centers$y[i] - half_y; top <- centers$y[i] + half_y
    data.frame(
      x1 = c(left, right, right, left), y1 = c(bottom, bottom, top, top),
      x2 = c(right, right, left, left), y2 = c(bottom, top, top, bottom)
    )
  }))
  point_key <- function(a, b) paste(sprintf("%.15g", a), sprintf("%.15g", b), sep = ",")
  start_key <- point_key(edges$x1, edges$y1)
  end_key <- point_key(edges$x2, edges$y2)
  undirected_key <- ifelse(start_key < end_key,
                           paste(start_key, end_key, sep = "|"),
                           paste(end_key, start_key, sep = "|"))
  boundary_edges <- edges[!duplicated(undirected_key) & !duplicated(undirected_key, fromLast = TRUE), , drop = FALSE]
  if (!nrow(boundary_edges)) stop("No exterior ROI edges could be constructed.", call. = FALSE)

  unused <- rep(TRUE, nrow(boundary_edges))
  rings <- list()
  while (any(unused)) {
    first <- which(unused)[1L]
    ring_indices <- first
    unused[first] <- FALSE
    first_key <- point_key(boundary_edges$x1[first], boundary_edges$y1[first])
    current_key <- point_key(boundary_edges$x2[first], boundary_edges$y2[first])
    while (current_key != first_key) {
      candidates <- which(unused & point_key(boundary_edges$x1, boundary_edges$y1) == current_key)
      if (!length(candidates)) {
        stop("ROI boundary is not a set of closed grid-aligned rings.", call. = FALSE)
      }
      next_edge <- candidates[1L]
      ring_indices <- c(ring_indices, next_edge)
      unused[next_edge] <- FALSE
      current_key <- point_key(boundary_edges$x2[next_edge], boundary_edges$y2[next_edge])
    }
    ring <- boundary_edges[ring_indices, , drop = FALSE]
    ring <- data.frame(
      vertex_order = seq_len(nrow(ring) + 1L),
      x = c(ring$x1, ring$x2[nrow(ring)]),
      y = c(ring$y1, ring$y2[nrow(ring)]),
      stringsAsFactors = FALSE
    )
    signed_area <- sum(
      ring$x[-nrow(ring)] * ring$y[-1L] -
        ring$x[-1L] * ring$y[-nrow(ring)]
    ) / 2
    ring$ring_type <- if (signed_area >= 0) "exterior" else "hole"
    rings[[length(rings) + 1L]] <- ring
  }
  out <- do.call(rbind, Map(function(ring, id) {
    ring$ring_id <- id
    ring[, c("ring_id", "ring_type", "vertex_order", "x", "y")]
  }, rings, seq_along(rings)))
  rownames(out) <- NULL
  out
}

export_lcm_targets <- function(roi_pixels,
                               registration,
                               x_resolution,
                               y_resolution = x_resolution,
                               roi_column = "roi_id",
                               section_column = NULL,
                               x_column = "x",
                               y_column = "y") {
  required_columns(
    roi_pixels, c(roi_column, section_column, x_column, y_column),
    "ROI pixel table"
  )
  if (!inherits(registration, "histology_msi_registration")) {
    stop("registration must come from fit_histology_msi_registration().", call. = FALSE)
  }
  group_columns <- c(section_column, roi_column)
  group_key <- do.call(interaction, c(roi_pixels[group_columns], list(drop = TRUE, lex.order = TRUE)))
  groups <- split(seq_len(nrow(roi_pixels)), group_key)
  polygons <- lapply(groups, function(index) {
    polygon <- roi_pixels_to_polygon(
      roi_pixels[[x_column]][index], roi_pixels[[y_column]][index],
      x_resolution, y_resolution
    )
    polygon[[roi_column]] <- as.character(roi_pixels[[roi_column]][index[1L]])
    if (!is.null(section_column)) {
      polygon[[section_column]] <- as.character(roi_pixels[[section_column]][index[1L]])
    }
    polygon
  })
  polygons <- do.call(rbind, polygons)
  rownames(polygons) <- NULL
  transformed <- transform_histology_coordinates(
    polygons, registration,
    x_column = "x", y_column = "y",
    output_columns = c("stage_x", "stage_y"),
    section_column = section_column
  )
  diagnostics <- registration_diagnostics(registration)
  if (is.null(section_column)) {
    transformed$registration_rmse <- diagnostics$rmse[1L]
    transformed$registration_max_error <- diagnostics$max_error[1L]
  } else {
    matched <- match(as.character(transformed[[section_column]]), as.character(diagnostics$section_id))
    transformed$registration_rmse <- diagnostics$rmse[matched]
    transformed$registration_max_error <- diagnostics$max_error[matched]
  }
  transformed$coordinate_unit <- "registration_target_unit"
  attr(transformed, "export_scope") <- paste(
    "Vendor-neutral coordinates only. Validate axes, units, fiducials, polygon",
    "topology, and instrument-specific constraints before physical cutting."
  )
  transformed
}

compare_msi_lcm_quantification <- function(
    msi_matrix,
    lcm_matrix,
    feature_mapping,
    id_columns = c("roi_id", "section_id"),
    msi_feature_column = "msi_feature",
    lcm_feature_column = "lcm_feature",
    method = c("spearman", "pearson"),
    subject_column = NULL,
    lcm_area_column = NULL,
    lcm_internal_standard_column = NULL,
    min_pairs = 3,
    p_adjust_method = "BH") {
  method <- match.arg(method)
  required_columns(msi_matrix, c(id_columns, subject_column), "MSI ROI matrix")
  required_columns(lcm_matrix, id_columns, "LCM quantification matrix")
  required_columns(feature_mapping, c(msi_feature_column, lcm_feature_column), "Feature mapping")
  if (!is.numeric(min_pairs) || length(min_pairs) != 1L || min_pairs < 3) {
    stop("min_pairs must be one integer of at least 3.", call. = FALSE)
  }
  if (!is.null(lcm_area_column)) required_columns(lcm_matrix, lcm_area_column, "LCM quantification matrix")
  if (!is.null(lcm_internal_standard_column)) {
    required_columns(lcm_matrix, lcm_internal_standard_column, "LCM quantification matrix")
  }
  duplicate_msi <- duplicated(msi_matrix[id_columns])
  duplicate_lcm <- duplicated(lcm_matrix[id_columns])
  if (any(duplicate_msi) || any(duplicate_lcm)) {
    stop("Each input must have at most one row per id_columns key.", call. = FALSE)
  }
  joined <- merge(msi_matrix, lcm_matrix, by = id_columns, suffixes = c(".msi", ".lcm"))
  if (!nrow(joined)) stop("MSI and LCM tables have no shared ROI keys.", call. = FALSE)
  joined_name <- function(column, source) {
    other <- if (source == "msi") lcm_matrix else msi_matrix
    if (column %in% names(other) && !column %in% id_columns) {
      paste0(column, if (source == "msi") ".msi" else ".lcm")
    } else column
  }
  area <- if (is.null(lcm_area_column)) rep(1, nrow(joined)) else {
    as.numeric(joined[[joined_name(lcm_area_column, "lcm")]])
  }
  internal_standard <- if (is.null(lcm_internal_standard_column)) rep(1, nrow(joined)) else {
    as.numeric(joined[[joined_name(lcm_internal_standard_column, "lcm")]])
  }
  denominator <- area * internal_standard
  if (any(!is.finite(denominator) | denominator <= 0)) {
    stop("LCM area and internal-standard normalization values must be positive and finite.", call. = FALSE)
  }

  rows <- lapply(seq_len(nrow(feature_mapping)), function(i) {
    msi_feature <- as.character(feature_mapping[[msi_feature_column]][i])
    lcm_feature <- as.character(feature_mapping[[lcm_feature_column]][i])
    if (!msi_feature %in% names(msi_matrix) || !lcm_feature %in% names(lcm_matrix)) {
      return(data.frame(msi_feature = msi_feature, lcm_feature = lcm_feature,
                        n_pairs = 0L, n_subjects = NA_integer_, correlation = NA_real_,
                        p_value = NA_real_, analysis_scale = "unavailable"))
    }
    msi_name <- joined_name(msi_feature, "msi")
    lcm_name <- joined_name(lcm_feature, "lcm")
    msi_value <- suppressWarnings(as.numeric(joined[[msi_name]]))
    lcm_value <- suppressWarnings(as.numeric(joined[[lcm_name]])) / denominator
    analysis_scale <- "across_roi"
    n_subjects <- NA_integer_
    if (!is.null(subject_column)) {
      subject_name <- joined_name(subject_column, "msi")
      subject <- as.character(joined[[subject_name]])
      n_subjects <- length(unique(subject[!is.na(subject)]))
      msi_value <- msi_value - stats::ave(msi_value, subject, FUN = function(z) mean(z, na.rm = TRUE))
      lcm_value <- lcm_value - stats::ave(lcm_value, subject, FUN = function(z) mean(z, na.rm = TRUE))
      analysis_scale <- "within_subject_centered"
    }
    keep <- is.finite(msi_value) & is.finite(lcm_value)
    n_pairs <- sum(keep)
    estimate <- p_value <- NA_real_
    if (n_pairs >= min_pairs && length(unique(msi_value[keep])) > 1L &&
        length(unique(lcm_value[keep])) > 1L) {
      test <- suppressWarnings(stats::cor.test(msi_value[keep], lcm_value[keep], method = method, exact = FALSE))
      estimate <- unname(test$estimate)
      p_value <- test$p.value
    }
    data.frame(
      msi_feature = msi_feature, lcm_feature = lcm_feature,
      n_pairs = n_pairs, n_subjects = n_subjects,
      correlation = estimate, p_value = p_value,
      analysis_scale = analysis_scale, stringsAsFactors = FALSE
    )
  })
  result <- do.call(rbind, rows)
  result$adj_p_value <- stats::p.adjust(result$p_value, method = p_adjust_method)
  result$method <- method
  result$lcm_area_normalized <- !is.null(lcm_area_column)
  result$lcm_internal_standard_normalized <- !is.null(lcm_internal_standard_column)
  result <- result[order(result$adj_p_value, -abs(result$correlation), na.last = TRUE), , drop = FALSE]
  rownames(result) <- NULL
  attr(result, "interpretation") <- paste(
    "Correlation assesses ROI-level concordance and does not establish identity.",
    "Within-subject centering controls subject offsets but is not a mixed-effects model."
  )
  result
}
