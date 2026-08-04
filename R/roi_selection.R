# Unified histology-guided ROI selection -----------------------------------

validate_roi_weights <- function(weights, names_expected, label) {
  if (is.null(weights)) weights <- rep(1, length(names_expected))
  if (!is.null(names(weights))) {
    missing_names <- setdiff(names_expected, names(weights))
    if (length(missing_names)) {
      stop(label, " is missing: ", paste(missing_names, collapse = ", "), call. = FALSE)
    }
    weights <- weights[names_expected]
  }
  weights <- as.numeric(weights)
  if (length(weights) != length(names_expected) || any(!is.finite(weights)) ||
      any(weights < 0) || sum(weights) <= 0) {
    stop(label, " must contain non-negative finite values with a positive sum.", call. = FALSE)
  }
  stats::setNames(weights / sum(weights), names_expected)
}

cluster_histology_superpixels <- function(superpixels,
                                          feature_columns,
                                          k,
                                          n_pcs = 80,
                                          cluster_column = "histology_cluster",
                                          center = TRUE,
                                          scale. = TRUE,
                                          nstart = 25,
                                          iter.max = 100,
                                          seed = 1) {
  required_columns(superpixels, c("x", "y", feature_columns), "Histology superpixel table")
  if (length(feature_columns) < 1) {
    stop("feature_columns must name at least one UNI/histology feature.", call. = FALSE)
  }
  x <- as.matrix(superpixels[, feature_columns, drop = FALSE])
  storage.mode(x) <- "double"
  for (j in seq_len(ncol(x))) {
    finite <- is.finite(x[, j])
    replacement <- if (any(finite)) stats::median(x[finite, j]) else 0
    x[!finite, j] <- replacement
  }
  variable <- apply(x, 2, stats::sd) > 0
  x <- x[, variable, drop = FALSE]
  if (ncol(x) < 1) stop("All histology features are constant.", call. = FALSE)

  max_pcs <- min(as.integer(n_pcs), ncol(x), nrow(x) - 1L)
  if (max_pcs < 1) stop("At least two superpixels are required.", call. = FALSE)
  pca <- stats::prcomp(x, center = center, scale. = scale., rank. = max_pcs)
  scores <- pca$x[, seq_len(min(max_pcs, ncol(pca$x))), drop = FALSE]
  k <- max(2L, min(as.integer(k), nrow(scores)))
  if (!is.null(seed)) set.seed(seed)
  fit <- stats::kmeans(
    scores,
    centers = k,
    nstart = nstart,
    iter.max = iter.max,
    algorithm = "Lloyd"
  )
  out <- superpixels
  out[[cluster_column]] <- fit$cluster
  list(
    matrix = out,
    fit = fit,
    pca = pca,
    pca_scores = scores,
    n_pcs = ncol(scores),
    feature_columns = feature_columns[variable]
  )
}

make_roi_target <- function(cluster_values, prior_weights = NULL) {
  cluster_names <- sort(unique(as.character(cluster_values[!is.na(cluster_values)])))
  if (!length(cluster_names)) stop("No non-missing histology clusters are available.", call. = FALSE)
  if (!is.null(prior_weights) && !is.null(names(prior_weights))) {
    unknown <- setdiff(names(prior_weights), cluster_names)
    if (length(unknown)) {
      stop("prior_weights contains unknown clusters: ", paste(unknown, collapse = ", "), call. = FALSE)
    }
    expanded <- stats::setNames(rep(1, length(cluster_names)), cluster_names)
    expanded[names(prior_weights)] <- prior_weights
    prior_weights <- expanded
  }
  validate_roi_weights(prior_weights, cluster_names, "prior_weights")
}

estimate_slide_area <- function(x, y) {
  positive_spacing <- function(v) {
    d <- diff(sort(unique(v[is.finite(v)])))
    d <- d[d > 0]
    if (length(d)) stats::median(d) else 1
  }
  (diff(range(x, na.rm = TRUE)) + positive_spacing(x)) *
    (diff(range(y, na.rm = TRUE)) + positive_spacing(y))
}

roi_membership <- function(x, y, center_x, center_y, roi_size, shape) {
  if (shape == "square") {
    half <- roi_size / 2
    abs(x - center_x) <= half & abs(y - center_y) <= half
  } else {
    (x - center_x)^2 + (y - center_y)^2 <= (roi_size / 2)^2
  }
}

generate_roi_candidates <- function(superpixels,
                                    roi_size,
                                    shape = c("square", "circle"),
                                    n_candidates = NULL,
                                    candidate_multiplier = 30,
                                    min_candidates = 200,
                                    max_candidates = 5000,
                                    seed = 1,
                                    section_column = NULL) {
  shape <- match.arg(shape)
  required_columns(superpixels, c("x", "y"), "Histology superpixel table")
  if (!is.numeric(roi_size) || length(roi_size) != 1 || !is.finite(roi_size) || roi_size <= 0) {
    stop("roi_size must be one positive value in the same coordinate unit as x and y.", call. = FALSE)
  }
  if (!is.null(section_column)) {
    required_columns(superpixels, section_column, "Histology superpixel table")
    if (length(unique(superpixels[[section_column]])) > 1) {
      stop("Automatic ROI selection currently supports one section at a time.", call. = FALSE)
    }
  }
  finite_rows <- is.finite(superpixels$x) & is.finite(superpixels$y)
  centers_source <- superpixels[finite_rows, c("x", "y"), drop = FALSE]
  if (!nrow(centers_source)) stop("No finite superpixel coordinates are available.", call. = FALSE)
  bounds_x <- range(centers_source$x)
  bounds_y <- range(centers_source$y)
  half <- roi_size / 2
  if (diff(bounds_x) >= roi_size && diff(bounds_y) >= roi_size) {
    inside_bounds <- centers_source$x >= bounds_x[1] + half &
      centers_source$x <= bounds_x[2] - half &
      centers_source$y >= bounds_y[1] + half &
      centers_source$y <= bounds_y[2] - half
    centers_source <- centers_source[inside_bounds, , drop = FALSE]
  }
  if (!nrow(centers_source)) {
    centers_source <- data.frame(x = mean(bounds_x), y = mean(bounds_y))
  }
  if (is.null(n_candidates)) {
    slide_area <- estimate_slide_area(centers_source$x, centers_source$y)
    roi_area <- if (shape == "square") roi_size^2 else pi * (roi_size / 2)^2
    n_candidates <- ceiling(candidate_multiplier * slide_area / roi_area)
    n_candidates <- max(min_candidates, min(max_candidates, n_candidates))
  }
  n_candidates <- max(1L, as.integer(n_candidates))
  if (!is.null(seed)) set.seed(seed)
  sampled_rows <- sample(seq_len(nrow(centers_source)), n_candidates, replace = n_candidates > nrow(centers_source))
  centers <- unique(centers_source[sampled_rows, , drop = FALSE])
  memberships <- lapply(seq_len(nrow(centers)), function(i) {
    which(roi_membership(
      superpixels$x,
      superpixels$y,
      centers$x[i],
      centers$y[i],
      roi_size,
      shape
    ))
  })
  keep <- lengths(memberships) > 0
  centers <- centers[keep, , drop = FALSE]
  memberships <- memberships[keep]
  centers$candidate_id <- sprintf("candidate_%05d", seq_len(nrow(centers)))
  centers$shape <- shape
  centers$roi_size <- roi_size
  centers$n_superpixels <- lengths(memberships)
  rownames(centers) <- NULL
  structure(list(candidates = centers, membership = memberships), class = "roi_candidates")
}

normalized_logistic <- function(x, midpoint = 0.1, steepness = 12) {
  lo <- stats::plogis(steepness * (0 - midpoint))
  hi <- stats::plogis(steepness * (1 - midpoint))
  pmax(0, pmin(1, (stats::plogis(steepness * (x - midpoint)) - lo) / (hi - lo)))
}

score_roi_indices <- function(superpixels,
                              candidates,
                              candidate_indices,
                              cluster_column = "histology_cluster",
                              valid_column = "qc_pass",
                              prior_weights = NULL,
                              score_weights = c(balance = 1, coverage = 1, size = 1),
                              size_midpoint = 0.1,
                              size_steepness = 12) {
  required_columns(superpixels, cluster_column, "Histology superpixel table")
  if (!is.null(valid_column)) required_columns(superpixels, valid_column, "Histology superpixel table")
  if (!inherits(candidates, "roi_candidates")) stop("candidates must come from generate_roi_candidates().", call. = FALSE)
  candidate_indices <- unique(as.integer(candidate_indices))
  if (!length(candidate_indices) || any(!candidate_indices %in% seq_along(candidates$membership))) {
    stop("candidate_indices must identify at least one ROI candidate.", call. = FALSE)
  }

  valid <- if (is.null(valid_column)) rep(TRUE, nrow(superpixels)) else {
    value <- as.logical(superpixels[[valid_column]])
    value[is.na(value)] <- FALSE
    value
  }
  selected_all <- unique(unlist(candidates$membership[candidate_indices], use.names = FALSE))
  selected_valid <- selected_all[valid[selected_all]]
  selected_clustered <- selected_valid[!is.na(superpixels[[cluster_column]][selected_valid])]
  denominator <- length(selected_all)
  coverage_raw <- if (denominator > 0) length(selected_valid) / denominator else 0
  coverage_score <- sqrt(max(0, min(1, coverage_raw)))

  target <- make_roi_target(superpixels[[cluster_column]], prior_weights)
  counts <- table(factor(
    as.character(superpixels[[cluster_column]][selected_clustered]),
    levels = names(target)
  ))
  observed <- as.numeric(counts)
  balance_score <- if (sum(observed) == 0) 0 else {
    sum(observed * target) / sqrt(sum(observed^2) * sum(target^2))
  }
  total_valid <- sum(valid)
  size_fraction <- if (total_valid > 0) length(selected_valid) / total_valid else 0
  size_score <- normalized_logistic(size_fraction, midpoint = size_midpoint, steepness = size_steepness)
  weights <- validate_roi_weights(score_weights, c("balance", "coverage", "size"), "score_weights")
  parts <- c(balance = balance_score, coverage = coverage_score, size = size_score)
  active <- weights > 0
  roi_score <- if (any(parts[active] <= 0)) 0 else exp(sum(weights[active] * log(parts[active])))

  c(
    roi_score = roi_score,
    balance_score = balance_score,
    coverage_score = coverage_score,
    size_score = size_score,
    coverage_fraction = coverage_raw,
    size_fraction = size_fraction,
    n_valid_superpixels = length(selected_valid),
    n_selected_superpixels = length(selected_all)
  )
}

candidate_overlap <- function(candidates, a, b) {
  ma <- candidates$membership[[a]]
  mb <- candidates$membership[[b]]
  if (!length(ma) || !length(mb)) return(0)
  length(intersect(ma, mb)) / min(length(ma), length(mb))
}

optimize_histology_rois <- function(superpixels,
                                    roi_size,
                                    shape = c("square", "circle"),
                                    cluster_column = "histology_cluster",
                                    valid_column = "qc_pass",
                                    prior_weights = NULL,
                                    score_weights = c(balance = 1, coverage = 1, size = 1),
                                    max_rois = 10,
                                    improvement_threshold = 0.03,
                                    n_candidates = NULL,
                                    candidate_multiplier = 30,
                                    beam_width = 25,
                                    candidate_pool = 250,
                                    max_overlap = 0,
                                    size_midpoint = 0.1,
                                    size_steepness = 12,
                                    seed = 1,
                                    section_column = NULL) {
  shape <- match.arg(shape)
  candidates <- generate_roi_candidates(
    superpixels,
    roi_size = roi_size,
    shape = shape,
    n_candidates = n_candidates,
    candidate_multiplier = candidate_multiplier,
    seed = seed,
    section_column = section_column
  )
  n <- nrow(candidates$candidates)
  if (!n) stop("No ROI candidates contain superpixels.", call. = FALSE)
  score_one <- function(indices) score_roi_indices(
    superpixels,
    candidates,
    indices,
    cluster_column = cluster_column,
    valid_column = valid_column,
    prior_weights = prior_weights,
    score_weights = score_weights,
    size_midpoint = size_midpoint,
    size_steepness = size_steepness
  )
  individual_scores <- vapply(seq_len(n), function(i) score_one(i)[["roi_score"]], numeric(1))
  pool <- order(individual_scores, decreasing = TRUE)[seq_len(min(n, max(1L, as.integer(candidate_pool))))]
  beams <- lapply(pool[seq_len(min(length(pool), max(1L, as.integer(beam_width))))], function(i) i)
  beam_scores <- vapply(beams, function(z) score_one(z)[["roi_score"]], numeric(1))
  best_at_r <- beams[[which.max(beam_scores)]]
  best_details <- score_one(best_at_r)
  history <- data.frame(
    n_rois = 1L,
    roi_score = best_details[["roi_score"]],
    improvement = NA_real_,
    accepted = TRUE,
    stringsAsFactors = FALSE
  )
  accepted <- best_at_r
  accepted_details <- best_details

  compatible <- function(current, addition) {
    all(vapply(current, function(old) candidate_overlap(candidates, old, addition) <= max_overlap, logical(1)))
  }
  if (max_rois > 1) {
    for (r in 2:max(1L, as.integer(max_rois))) {
      expanded <- list()
      expanded_scores <- numeric()
      seen <- character()
      for (current in beams) {
        additions <- pool[!pool %in% current]
        for (addition in additions) {
          if (!compatible(current, addition)) next
          proposal <- sort(c(current, addition))
          key <- paste(proposal, collapse = ",")
          if (key %in% seen) next
          seen <- c(seen, key)
          details <- score_one(proposal)
          expanded[[length(expanded) + 1L]] <- proposal
          expanded_scores <- c(expanded_scores, details[["roi_score"]])
        }
      }
      if (!length(expanded)) break
      keep <- order(expanded_scores, decreasing = TRUE)[seq_len(min(length(expanded), max(1L, as.integer(beam_width))))]
      beams <- expanded[keep]
      beam_scores <- expanded_scores[keep]
      proposed <- beams[[1]]
      proposed_details <- score_one(proposed)
      improvement <- proposed_details[["roi_score"]] - accepted_details[["roi_score"]]
      is_accepted <- is.finite(improvement) && improvement >= improvement_threshold
      history <- rbind(history, data.frame(
        n_rois = r,
        roi_score = proposed_details[["roi_score"]],
        improvement = improvement,
        accepted = is_accepted,
        stringsAsFactors = FALSE
      ))
      if (!is_accepted) break
      accepted <- proposed
      accepted_details <- proposed_details
    }
  }
  selected_candidates <- candidates$candidates[accepted, , drop = FALSE]
  selected_candidates$roi_id <- sprintf("roi_%02d", seq_len(nrow(selected_candidates)))
  selected_candidates$score <- accepted_details[["roi_score"]]
  list(
    selected_indices = accepted,
    selected_candidates = selected_candidates,
    score = accepted_details,
    history = history,
    candidates = candidates,
    target_proportions = make_roi_target(superpixels[[cluster_column]], prior_weights),
    settings = list(
      roi_size = roi_size,
      shape = shape,
      improvement_threshold = improvement_threshold,
      max_overlap = max_overlap,
      score_weights = score_weights,
      prior_weights = prior_weights
    )
  )
}

point_in_polygon <- function(x, y, polygon_x, polygon_y) {
  n <- length(polygon_x)
  if (n < 3) return(rep(FALSE, length(x)))
  inside <- rep(FALSE, length(x))
  j <- n
  for (i in seq_len(n)) {
    yi <- polygon_y[i]
    yj <- polygon_y[j]
    crosses <- (yi > y) != (yj > y)
    x_intersect <- (polygon_x[j] - polygon_x[i]) * (y - yi) / (yj - yi + .Machine$double.eps) + polygon_x[i]
    inside <- xor(inside, crosses & x <= x_intersect)
    j <- i
  }
  inside
}

manual_roi_masks <- function(superpixels,
                             method = c("cluster", "geometry", "polygon", "combined"),
                             cluster_column = "histology_cluster",
                             selected_clusters = NULL,
                             roi_table = NULL,
                             polygon_vertices = NULL,
                             section_column = NULL) {
  method <- match.arg(method)
  required_columns(superpixels, c("x", "y"), "Spatial pixel table")
  section_values <- if (is.null(section_column)) rep(".__single__", nrow(superpixels)) else {
    required_columns(superpixels, section_column, "Spatial pixel table")
    as.character(superpixels[[section_column]])
  }
  masks <- list()
  ids <- character()

  if (method == "cluster") {
    required_columns(superpixels, cluster_column, "Spatial pixel table")
    if (!length(selected_clusters)) stop("selected_clusters is required for manual cluster selection.", call. = FALSE)
    sections <- unique(section_values)
    for (section in sections) {
      for (cluster in selected_clusters) {
        prefix <- if (is.null(section_column)) "" else paste0(utils::URLencode(section, reserved = TRUE), "__")
        ids <- c(ids, paste0(prefix, "roi_cluster_", cluster))
        masks[[length(masks) + 1L]] <- section_values == section &
          !is.na(superpixels[[cluster_column]]) &
          as.character(superpixels[[cluster_column]]) == as.character(cluster)
      }
    }
    return(list(roi_id = make.unique(ids), masks = masks))
  }

  if (method == "combined") {
    geometry_masks <- if (!is.null(roi_table) && nrow(roi_table)) {
      manual_roi_masks(
        superpixels, method = "geometry", cluster_column = cluster_column,
        roi_table = roi_table, section_column = section_column
      )
    } else list(roi_id = character(), masks = list())
    polygon_masks <- if (!is.null(polygon_vertices) && nrow(polygon_vertices)) {
      manual_roi_masks(
        superpixels, method = "polygon", cluster_column = cluster_column,
        polygon_vertices = polygon_vertices, section_column = section_column
      )
    } else list(roi_id = character(), masks = list())
    if (!length(geometry_masks$masks) && !length(polygon_masks$masks)) {
      stop("At least one geometry or polygon ROI is required.", call. = FALSE)
    }
    return(list(
      roi_id = make.unique(c(geometry_masks$roi_id, polygon_masks$roi_id)),
      masks = c(geometry_masks$masks, polygon_masks$masks)
    ))
  }

  definitions <- if (method == "geometry") roi_table else polygon_vertices
  if (is.null(definitions)) {
    stop(if (method == "geometry") "roi_table is required." else "polygon_vertices is required.", call. = FALSE)
  }
  required_columns(definitions, "roi_id", if (method == "geometry") "ROI table" else "Polygon vertices")
  definition_sections <- if (is.null(section_column)) rep(".__single__", nrow(definitions)) else {
    required_columns(definitions, section_column, if (method == "geometry") "ROI table" else "Polygon vertices")
    as.character(definitions[[section_column]])
  }

  if (method == "geometry") {
    required_columns(definitions, "shape", "ROI table")
    for (i in seq_len(nrow(definitions))) {
      shape <- match.arg(as.character(definitions$shape[i]), c("square", "rectangle", "circle"))
      in_section <- section_values == definition_sections[i]
      if (shape == "rectangle") {
        required_columns(definitions, c("x_min", "x_max", "y_min", "y_max"), "Rectangle ROI table")
        mask <- in_section & superpixels$x >= min(definitions$x_min[i], definitions$x_max[i]) &
          superpixels$x <= max(definitions$x_min[i], definitions$x_max[i]) &
          superpixels$y >= min(definitions$y_min[i], definitions$y_max[i]) &
          superpixels$y <= max(definitions$y_min[i], definitions$y_max[i])
      } else {
        required_columns(definitions, c("center_x", "center_y"), "ROI table")
        size <- if (shape == "circle") {
          required_columns(definitions, "radius", "Circle ROI table")
          2 * definitions$radius[i]
        } else {
          required_columns(definitions, "size", "Square ROI table")
          definitions$size[i]
        }
        mask <- in_section & roi_membership(
          superpixels$x,
          superpixels$y,
          definitions$center_x[i],
          definitions$center_y[i],
          size,
          if (shape == "circle") "circle" else "square"
        )
      }
      ids <- c(ids, as.character(definitions$roi_id[i]))
      masks[[length(masks) + 1L]] <- mask
    }
  } else {
    required_columns(definitions, c("x", "y"), "Polygon vertices")
    group_key <- interaction(as.character(definitions$roi_id), definition_sections, drop = TRUE, lex.order = TRUE)
    groups <- split(seq_len(nrow(definitions)), group_key)
    for (idx in groups) {
      if ("vertex_order" %in% names(definitions)) idx <- idx[order(definitions$vertex_order[idx])]
      in_section <- section_values == definition_sections[idx[1]]
      mask <- in_section & point_in_polygon(
        superpixels$x,
        superpixels$y,
        definitions$x[idx],
        definitions$y[idx]
      )
      ids <- c(ids, as.character(definitions$roi_id[idx[1]]))
      masks[[length(masks) + 1L]] <- mask
    }
  }
  list(roi_id = make.unique(ids), masks = masks)
}

select_manual_rois <- function(superpixels,
                               method = c("cluster", "geometry", "polygon", "combined"),
                               cluster_column = "histology_cluster",
                               selected_clusters = NULL,
                               roi_table = NULL,
                               polygon_vertices = NULL,
                               section_column = NULL,
                               roi_column = "roi_id",
                               overlap = c("error", "first", "last")) {
  method <- match.arg(method)
  overlap <- match.arg(overlap)
  definitions <- manual_roi_masks(
    superpixels,
    method = method,
    cluster_column = cluster_column,
    selected_clusters = selected_clusters,
    roi_table = roi_table,
    polygon_vertices = polygon_vertices,
    section_column = section_column
  )
  membership_count <- Reduce(`+`, lapply(definitions$masks, as.integer), init = integer(nrow(superpixels)))
  if (overlap == "error" && any(membership_count > 1)) {
    stop("Manual ROI definitions overlap. Use overlap = 'first' or 'last' to set precedence.", call. = FALSE)
  }
  out <- superpixels
  out[[roi_column]] <- NA_character_
  order_indices <- if (overlap == "first") rev(seq_along(definitions$masks)) else seq_along(definitions$masks)
  for (i in order_indices) out[[roi_column]][definitions$masks[[i]]] <- definitions$roi_id[i]
  out$roi_selection_method <- paste0("manual_", method)
  selected <- out[!is.na(out[[roi_column]]) & out[[roi_column]] != "", , drop = FALSE]
  rownames(selected) <- NULL
  list(
    selected_pixels = selected,
    annotated_pixels = out,
    roi_summary = as.data.frame(table(selected[[roi_column]]), stringsAsFactors = FALSE)
  )
}

select_rois <- function(superpixels,
                        selection_mode = c("automatic", "manual"),
                        manual_method = c("cluster", "geometry", "polygon", "combined"),
                        roi_size = NULL,
                        shape = c("square", "circle"),
                        cluster_column = "histology_cluster",
                        valid_column = "qc_pass",
                        selected_clusters = NULL,
                        roi_table = NULL,
                        polygon_vertices = NULL,
                        section_column = NULL,
                        roi_column = "roi_id",
                        overlap = c("error", "first", "last"),
                        ...) {
  selection_mode <- match.arg(selection_mode)
  if (selection_mode == "manual") {
    return(select_manual_rois(
      superpixels,
      method = match.arg(manual_method),
      cluster_column = cluster_column,
      selected_clusters = selected_clusters,
      roi_table = roi_table,
      polygon_vertices = polygon_vertices,
      section_column = section_column,
      roi_column = roi_column,
      overlap = match.arg(overlap)
    ))
  }
  if (is.null(roi_size)) stop("roi_size is required for automatic ROI selection.", call. = FALSE)
  optimization <- optimize_histology_rois(
    superpixels,
    roi_size = roi_size,
    shape = match.arg(shape),
    cluster_column = cluster_column,
    valid_column = valid_column,
    section_column = section_column,
    ...
  )
  out <- superpixels
  out[[roi_column]] <- NA_character_
  for (i in seq_along(optimization$selected_indices)) {
    candidate_index <- optimization$selected_indices[i]
    out[[roi_column]][optimization$candidates$membership[[candidate_index]]] <-
      optimization$selected_candidates$roi_id[i]
  }
  out$roi_selection_method <- "automatic"
  selected <- out[!is.na(out[[roi_column]]), , drop = FALSE]
  rownames(selected) <- NULL
  list(
    selected_pixels = selected,
    annotated_pixels = out,
    roi_summary = optimization$selected_candidates,
    optimization = optimization
  )
}

apply_roi_labels <- function(pixel_matrix,
                             selected_pixel_coordinates,
                             roi_column = "roi_id",
                             section_column = NULL,
                             pixel_id_column = "pixel_id") {
  required_columns(selected_pixel_coordinates, roi_column, "Selected ROI coordinates")
  by <- character()
  if (!is.null(section_column)) {
    required_columns(pixel_matrix, section_column, "Pixel matrix")
    required_columns(selected_pixel_coordinates, section_column, "Selected ROI coordinates")
    by <- c(by, section_column)
  }
  if (pixel_id_column %in% names(pixel_matrix) && pixel_id_column %in% names(selected_pixel_coordinates)) {
    by <- c(by, pixel_id_column)
  } else {
    required_columns(pixel_matrix, c("x", "y"), "Pixel matrix")
    required_columns(selected_pixel_coordinates, c("x", "y"), "Selected ROI coordinates")
    by <- c(by, "x", "y")
  }
  labels <- unique(selected_pixel_coordinates[, c(by, roi_column), drop = FALSE])
  duplicate_keys <- duplicated(labels[, by, drop = FALSE])
  if (any(duplicate_keys)) stop("Selected coordinates assign more than one ROI to the same pixel key.", call. = FALSE)
  original_order <- seq_len(nrow(pixel_matrix))
  work <- pixel_matrix
  work$.__roi_order__ <- original_order
  if (roi_column %in% names(work)) work[[roi_column]] <- NULL
  out <- merge(work, labels, by = by, all.x = TRUE, sort = FALSE)
  out <- out[order(out$.__roi_order__), , drop = FALSE]
  out$.__roi_order__ <- NULL
  rownames(out) <- NULL
  out
}
