required_columns <- function(data, columns, label) {
  missing <- setdiff(columns, names(data))
  if (length(missing) > 0) {
    stop(label, " is missing required column(s): ", paste(missing, collapse = ", "), call. = FALSE)
  }
}

feature_columns <- function(data) {
  grep("^mz_", names(data), value = TRUE)
}

infer_feature_mapping <- function(pixel_matrix) {
  fcols <- feature_columns(pixel_matrix)
  mz_values <- as.numeric(sub("^mz_", "", fcols))
  data.frame(
    feature_id = seq_along(fcols),
    mzmine_id = seq_along(fcols),
    mz = mz_values,
    column_name = fcols,
    stringsAsFactors = FALSE
  )
}

normalize_mz_column_names <- function(columns) {
  out <- columns
  is_mz <- suppressWarnings(is.finite(as.numeric(columns)))
  out[is_mz] <- paste0("mz_", columns[is_mz])
  out <- gsub("-", "_", out, fixed = TRUE)
  out
}

import_peakpicked_msi_csv <- function(path,
                                      section_id = NULL,
                                      section_order = 1,
                                      sample_id = NULL) {
  data <- utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  required_columns(data, c("x", "y"), "Peak-picked MSI CSV")
  names(data) <- normalize_mz_column_names(names(data))

  if (is.null(section_id)) {
    section_id <- tools::file_path_sans_ext(basename(path))
  }
  if (is.null(sample_id)) {
    sample_id <- section_id
  }

  data$pixel_id <- seq_len(nrow(data))
  data$local_pixel_id <- data$pixel_id
  data$section_id <- section_id
  data$section_order <- section_order
  data$sample_id_source <- sample_id

  metadata_columns <- c("pixel_id", "local_pixel_id", "section_id", "section_order", "sample_id_source", "x", "y")
  data[, c(metadata_columns, setdiff(names(data), metadata_columns)), drop = FALSE]
}

load_peakpicked_msi_series <- function(paths,
                                       section_ids = NULL,
                                       section_order = NULL,
                                       sample_ids = NULL) {
  if (length(paths) == 0) {
    stop("paths must contain at least one peak-picked MSI CSV.", call. = FALSE)
  }

  n_sections <- length(paths)
  if (is.null(section_ids)) section_ids <- tools::file_path_sans_ext(basename(paths))
  if (is.null(section_order)) section_order <- seq_len(n_sections)
  if (is.null(sample_ids)) sample_ids <- section_ids
  if (length(section_ids) != n_sections || length(section_order) != n_sections || length(sample_ids) != n_sections) {
    stop("section_ids, section_order, and sample_ids must match length(paths).", call. = FALSE)
  }

  rows <- lapply(seq_along(paths), function(i) {
    import_peakpicked_msi_csv(
      path = paths[i],
      section_id = section_ids[i],
      section_order = section_order[i],
      sample_id = sample_ids[i]
    )
  })

  common_features <- Reduce(intersect, lapply(rows, feature_columns))
  if (length(common_features) == 0) {
    stop("No shared mz_* feature columns were found across the input CSV files.", call. = FALSE)
  }

  metadata_columns <- Reduce(union, lapply(rows, function(data) setdiff(names(data), feature_columns(data))))
  rows <- lapply(rows, function(data) {
    missing_metadata <- setdiff(metadata_columns, names(data))
    for (column in missing_metadata) data[[column]] <- NA
    data[, c(metadata_columns, common_features), drop = FALSE]
  })

  combined <- do.call(rbind, rows)
  combined$pixel_id <- seq_len(nrow(combined))
  rownames(combined) <- NULL

  list(
    pixel_matrix = combined,
    feature_mapping = infer_feature_mapping(combined),
    section_mapping = data.frame(
      section_id = section_ids,
      section_order = section_order,
      sample_id_source = sample_ids,
      path = paths,
      stringsAsFactors = FALSE
    )
  )
}

safe_cv <- function(x) {
  mu <- mean(x, na.rm = TRUE)
  if (!is.finite(mu) || mu == 0) return(0)
  stats::sd(x, na.rm = TRUE) / mu
}

otsu_threshold <- function(x, nbins = 256, return_details = FALSE) {
  x <- x[is.finite(x)]
  if (length(x) == 0 || length(unique(x)) < 2) {
    threshold <- if (length(x) == 0) NA_real_ else x[1]
    if (return_details) {
      return(list(
        threshold = threshold,
        actual_bins = 0L,
        hist = NULL,
        between_var = numeric()
      ))
    }
    return(threshold)
  }

  h <- hist(x, breaks = nbins, plot = FALSE)
  counts <- as.numeric(h$counts)
  mids <- as.numeric(h$mids)

  total <- sum(counts)
  sum_total <- sum(counts * mids)

  weight_bg <- cumsum(counts)
  weight_fg <- total - weight_bg

  sum_bg <- cumsum(counts * mids)
  sum_fg <- sum_total - sum_bg

  mean_bg <- sum_bg / weight_bg
  mean_fg <- sum_fg / weight_fg

  between_var <- weight_bg * weight_fg * (mean_bg - mean_fg)^2
  between_var[!is.finite(between_var)] <- NA

  threshold <- mids[which.max(between_var)]
  if (!is.finite(threshold)) threshold <- stats::median(x, na.rm = TRUE)

  if (return_details) {
    return(list(
      threshold = threshold,
      actual_bins = length(counts),
      hist = h,
      between_var = between_var
    ))
  }

  threshold
}

detect_ubiquitous_gap_features <- function(nonzero_rate,
                                           min_gap_multiplier = 6) {
  values <- nonzero_rate[is.finite(nonzero_rate)]
  if (length(values) < 3 || length(unique(values)) < 3) {
    return(list(
      threshold = NA_real_,
      gap_size = NA_real_,
      gap_cutoff = NA_real_,
      excluded = character(),
      n_excluded = 0L
    ))
  }

  xs <- sort(unique(values), decreasing = TRUE)
  gaps <- xs[-length(xs)] - xs[-1]
  gap_median <- stats::median(gaps, na.rm = TRUE)
  gap_mad <- stats::mad(gaps, constant = 1, na.rm = TRUE)
  gap_scale <- if (is.finite(gap_mad) && gap_mad > 0) gap_mad else stats::sd(gaps, na.rm = TRUE)
  if (!is.finite(gap_scale) || gap_scale == 0) {
    return(list(
      threshold = NA_real_,
      gap_size = NA_real_,
      gap_cutoff = NA_real_,
      excluded = character(),
      n_excluded = 0L
    ))
  }

  gap_cutoff <- gap_median + min_gap_multiplier * gap_scale
  candidate_idx <- which(gaps > gap_cutoff)
  if (length(candidate_idx) == 0) {
    return(list(
      threshold = NA_real_,
      gap_size = NA_real_,
      gap_cutoff = gap_cutoff,
      excluded = character(),
      n_excluded = 0L
    ))
  }

  cut_idx <- candidate_idx[1]
  threshold <- (xs[cut_idx] + xs[cut_idx + 1]) / 2
  excluded <- names(nonzero_rate)[is.finite(nonzero_rate) & nonzero_rate >= threshold]

  list(
    threshold = threshold,
    gap_size = gaps[cut_idx],
    gap_cutoff = gap_cutoff,
    excluded = excluded,
    n_excluded = length(excluded)
  )
}

detect_ubiquitous_mad_features <- function(nonzero_rate,
                                           z_cutoff = 3.5) {
  values <- nonzero_rate[is.finite(nonzero_rate)]
  med <- stats::median(values, na.rm = TRUE)
  mad_value <- stats::mad(values, constant = 1, na.rm = TRUE)
  if (!is.finite(mad_value) || mad_value == 0) {
    return(list(
      threshold = NA_real_,
      median = med,
      mad = mad_value,
      excluded = character(),
      n_excluded = 0L
    ))
  }

  modified_z <- 0.6745 * (nonzero_rate - med) / mad_value
  excluded <- names(nonzero_rate)[is.finite(modified_z) & modified_z > z_cutoff]

  list(
    threshold = med + z_cutoff * mad_value / 0.6745,
    median = med,
    mad = mad_value,
    excluded = excluded,
    n_excluded = length(excluded),
    modified_z = modified_z
  )
}

detect_ubiquitous_features <- function(pixel_matrix,
                                       fcols = feature_columns(pixel_matrix),
                                       method = c("gap", "mad", "consensus", "none"),
                                       gap_min_multiplier = 6,
                                       mad_z_cutoff = 3.5,
                                       min_agreement = 2) {
  method <- match.arg(method)
  if (length(fcols) == 0 || method == "none") {
    return(list(
      method = method,
      nonzero_rate = numeric(),
      excluded = character(),
      n_excluded = 0L,
      gap = NULL,
      mad = NULL,
      agreement = data.frame()
    ))
  }

  nonzero_rate <- vapply(fcols, function(feature) {
    mean(pixel_matrix[[feature]] > 0, na.rm = TRUE)
  }, numeric(1))

  gap <- detect_ubiquitous_gap_features(nonzero_rate, min_gap_multiplier = gap_min_multiplier)
  mad <- detect_ubiquitous_mad_features(nonzero_rate, z_cutoff = mad_z_cutoff)

  if (method == "gap") {
    excluded <- gap$excluded
  } else if (method == "mad") {
    excluded <- mad$excluded
  } else {
    feature_set <- unique(c(gap$excluded, mad$excluded))
    agreement_count <- vapply(feature_set, function(feature) {
      sum(feature %in% gap$excluded, feature %in% mad$excluded)
    }, integer(1))
    excluded <- names(agreement_count)[agreement_count >= min_agreement]
  }

  agreement <- data.frame(
    column_name = fcols,
    nonzero_rate = unname(nonzero_rate[fcols]),
    gap_outlier = fcols %in% gap$excluded,
    mad_outlier = fcols %in% mad$excluded,
    agreement_count = as.integer(fcols %in% gap$excluded) + as.integer(fcols %in% mad$excluded),
    stringsAsFactors = FALSE
  )

  list(
    method = method,
    nonzero_rate = nonzero_rate,
    excluded = excluded,
    n_excluded = length(excluded),
    gap = gap,
    mad = mad,
    agreement = agreement
  )
}

build_spatial_neighbors <- function(x,
                                    y,
                                    method = c("rook", "queen", "distance"),
                                    distance_threshold = NULL,
                                    weights = c("binary"),
                                    symmetric = TRUE) {
  method <- match.arg(method)
  weights <- match.arg(weights)
  x <- as.numeric(x)
  y <- as.numeric(y)
  if (length(x) != length(y) || length(x) == 0L) {
    stop("x and y must have the same positive length.", call. = FALSE)
  }
  if (any(!is.finite(x) | !is.finite(y))) {
    stop("Spatial coordinates must be finite and non-missing.", call. = FALSE)
  }
  if (anyDuplicated(paste(x, y, sep = "\r"))) {
    stop("Spatial coordinates must be unique.", call. = FALSE)
  }
  if (length(symmetric) != 1L || is.na(symmetric)) {
    stop("symmetric must be TRUE or FALSE.", call. = FALSE)
  }
  symmetric <- isTRUE(symmetric)
  default_threshold <- switch(method, rook = 1, queen = sqrt(2), distance = NA_real_)
  if (is.null(distance_threshold)) distance_threshold <- default_threshold
  if (length(distance_threshold) != 1L || !is.finite(distance_threshold) || distance_threshold <= 0) {
    stop("distance_threshold must be supplied as one positive finite value.", call. = FALSE)
  }

  n_nodes <- length(x)
  edge_from <- integer()
  edge_to <- integer()
  edge_distance <- numeric()
  if (method %in% c("rook", "queen")) {
    offsets <- if (method == "rook") {
      data.frame(dx = c(1, 0), dy = c(0, 1))
    } else {
      data.frame(dx = c(1, 0, 1, 1), dy = c(0, 1, 1, -1))
    }
    offsets$distance <- sqrt(offsets$dx^2 + offsets$dy^2)
    offsets <- offsets[offsets$distance <= distance_threshold + sqrt(.Machine$double.eps), , drop = FALSE]
    key_to_index <- stats::setNames(seq_len(n_nodes), paste(x, y, sep = "\r"))
    for (offset_index in seq_len(nrow(offsets))) {
      neighbor_index <- unname(key_to_index[paste(x + offsets$dx[offset_index], y + offsets$dy[offset_index], sep = "\r")])
      present <- !is.na(neighbor_index)
      edge_from <- c(edge_from, which(present))
      edge_to <- c(edge_to, neighbor_index[present])
      edge_distance <- c(edge_distance, rep(offsets$distance[offset_index], sum(present)))
    }
  } else if (n_nodes >= 2L) {
    from_rows <- vector("list", n_nodes - 1L)
    to_rows <- vector("list", n_nodes - 1L)
    distance_rows <- vector("list", n_nodes - 1L)
    for (i in seq_len(n_nodes - 1L)) {
      candidates <- seq.int(i + 1L, n_nodes)
      distances <- sqrt((x[candidates] - x[i])^2 + (y[candidates] - y[i])^2)
      keep <- distances <= distance_threshold + sqrt(.Machine$double.eps)
      from_rows[[i]] <- rep.int(i, sum(keep))
      to_rows[[i]] <- candidates[keep]
      distance_rows[[i]] <- distances[keep]
    }
    edge_from <- unlist(from_rows, use.names = FALSE)
    edge_to <- unlist(to_rows, use.names = FALSE)
    edge_distance <- unlist(distance_rows, use.names = FALSE)
  }
  if (any(edge_from == edge_to)) stop("Internal error: neighbor graph contains a self-loop.", call. = FALSE)
  undirected_edges <- data.frame(
    from = as.integer(edge_from), to = as.integer(edge_to),
    distance = as.numeric(edge_distance), weight = rep(1, length(edge_from)),
    stringsAsFactors = FALSE
  )
  directed_edges <- if (symmetric && nrow(undirected_edges) > 0L) {
    rbind(undirected_edges, data.frame(
      from = undirected_edges$to, to = undirected_edges$from,
      distance = undirected_edges$distance, weight = undirected_edges$weight,
      stringsAsFactors = FALSE
    ))
  } else undirected_edges
  degree <- tabulate(c(undirected_edges$from, undirected_edges$to), nbins = n_nodes)
  effective <- degree > 0L

  adjacency <- vector("list", n_nodes)
  if (nrow(undirected_edges) > 0L) {
    for (edge_index in seq_len(nrow(undirected_edges))) {
      a <- undirected_edges$from[edge_index]; b <- undirected_edges$to[edge_index]
      adjacency[[a]] <- c(adjacency[[a]], b)
      adjacency[[b]] <- c(adjacency[[b]], a)
    }
  }
  component_id <- integer(n_nodes)
  component_sizes <- integer()
  component_count <- 0L
  for (start in seq_len(n_nodes)) {
    if (component_id[start] != 0L) next
    component_count <- component_count + 1L
    queue <- integer(n_nodes); queue[1] <- start
    component_id[start] <- component_count
    head <- 1L; tail <- 1L; size <- 0L
    while (head <= tail) {
      node <- queue[head]; head <- head + 1L; size <- size + 1L
      for (neighbor in adjacency[[node]]) {
        if (component_id[neighbor] == 0L) {
          tail <- tail + 1L; queue[tail] <- neighbor
          component_id[neighbor] <- component_count
        }
      }
    }
    component_sizes[component_count] <- size
  }
  degree_distribution <- as.data.frame(table(degree), stringsAsFactors = FALSE)
  names(degree_distribution) <- c("degree", "node_count")
  degree_distribution$degree <- as.integer(as.character(degree_distribution$degree))
  component_table <- data.frame(
    component_id = seq_along(component_sizes), component_size = component_sizes,
    contains_edge = component_sizes > 1L, stringsAsFactors = FALSE
  )
  summary <- data.frame(
    method = method, weights = weights, symmetric = symmetric,
    n_nodes = n_nodes, n_effective = sum(effective), n_isolated = sum(!effective),
    isolated_proportion = mean(!effective), n_edges_undirected = nrow(undirected_edges),
    connected_component_count = length(component_sizes),
    effective_component_count = sum(component_sizes > 1L),
    distance_threshold = distance_threshold,
    includes_diagonal_neighbors = method == "queen" && distance_threshold >= sqrt(2) - sqrt(.Machine$double.eps),
    stringsAsFactors = FALSE
  )
  structure(list(
    edges = directed_edges, undirected_edges = undirected_edges,
    degree = degree, effective = effective, component_id = component_id,
    component_sizes = component_table, degree_distribution = degree_distribution,
    summary = summary, x = x, y = y, method = method, weights = weights,
    symmetric = symmetric, distance_threshold = distance_threshold
  ), class = "spatial_neighbors")
}

spatial_neighbor_diagnostics <- function(neighbors) {
  if (!inherits(neighbors, "spatial_neighbors")) stop("neighbors must come from build_spatial_neighbors().", call. = FALSE)
  summary_rows <- data.frame(
    section = "summary", metric = names(neighbors$summary),
    value = vapply(neighbors$summary, function(value) as.character(value[1]), character(1)),
    stringsAsFactors = FALSE
  )
  degree_rows <- data.frame(
    section = "degree_distribution",
    metric = paste0("degree_", neighbors$degree_distribution$degree),
    value = as.character(neighbors$degree_distribution$node_count), stringsAsFactors = FALSE
  )
  component_rows <- data.frame(
    section = "component_sizes",
    metric = sprintf("component_%04d", neighbors$component_sizes$component_id),
    value = as.character(neighbors$component_sizes$component_size), stringsAsFactors = FALSE
  )
  rbind(summary_rows, degree_rows, component_rows)
}

compute_morans_i_grid <- function(values,
                                  x = NULL,
                                  y = NULL,
                                  n_perm = 0,
                                  alternative = c("greater", "two.sided"),
                                  seed = NULL,
                                  neighbors = NULL,
                                  neighbor_method = c("rook", "queen", "distance"),
                                  distance_threshold = NULL,
                                  weights = "binary",
                                  symmetric = TRUE) {
  alternative <- match.arg(alternative)
  neighbor_method <- match.arg(neighbor_method)
  values <- as.numeric(values)
  if (is.null(neighbors)) {
    if (is.null(x) || is.null(y)) stop("x and y are required when neighbors is NULL.", call. = FALSE)
    neighbors <- build_spatial_neighbors(x, y, method = neighbor_method,
      distance_threshold = distance_threshold, weights = weights, symmetric = symmetric)
  }
  if (!inherits(neighbors, "spatial_neighbors")) stop("neighbors must come from build_spatial_neighbors().", call. = FALSE)
  n_total <- length(values)
  if (n_total != length(neighbors$x)) stop("values must contain one value per spatial-neighbor node.", call. = FALSE)
  edges <- neighbors$undirected_edges
  finite <- is.finite(values)
  if (nrow(edges) > 0L) edges <- edges[finite[edges$from] & finite[edges$to], , drop = FALSE]
  degree <- tabulate(c(edges$from, edges$to), nbins = n_total)
  effective <- finite & degree > 0L
  n_effective <- sum(effective); n_isolated <- n_total - n_effective; n_edges <- nrow(edges)
  if (n_edges == 0L) stop("Moran's I is undefined because the effective graph has no edges.", call. = FALSE)
  if (n_effective < 3L) stop("Moran's I requires at least three non-isolated finite nodes.", call. = FALSE)
  index_map <- integer(n_total); index_map[effective] <- seq_len(n_effective)
  edge_i <- index_map[edges$from]; edge_j <- index_map[edges$to]; edge_weight <- edges$weight
  effective_values <- values[effective]
  moran_stat <- function(v) {
    centered <- v - mean(v); denominator <- sum(centered^2)
    if (!is.finite(denominator) || denominator == 0) return(NA_real_)
    numerator <- 2 * sum(edge_weight * centered[edge_i] * centered[edge_j])
    (n_effective / (2 * sum(edge_weight))) * numerator / denominator
  }
  observed <- moran_stat(effective_values)
  p_value <- NA_real_
  if (n_perm > 0L && is.finite(observed)) {
    if (!is.null(seed)) set.seed(seed)
    permuted <- replicate(as.integer(n_perm), moran_stat(sample(effective_values, n_effective, replace = FALSE)))
    if (alternative == "greater") {
      p_value <- (sum(permuted >= observed, na.rm = TRUE) + 1) / (sum(is.finite(permuted)) + 1)
    } else {
      p_value <- (sum(abs(permuted) >= abs(observed), na.rm = TRUE) + 1) / (sum(is.finite(permuted)) + 1)
    }
  }
  list(I = observed, p_value = p_value, n = n_effective, n_total = n_total,
    n_effective = n_effective, n_isolated = n_isolated, n_edges = n_edges)
}

compute_gearys_c_grid <- function(values,
                                  x = NULL,
                                  y = NULL,
                                  n_perm = 0,
                                  alternative = c("less", "greater", "two.sided"),
                                  seed = NULL,
                                  neighbors = NULL,
                                  neighbor_method = c("rook", "queen", "distance"),
                                  distance_threshold = NULL,
                                  weights = "binary",
                                  symmetric = TRUE) {
  alternative <- match.arg(alternative)
  neighbor_method <- match.arg(neighbor_method)
  if (length(n_perm) != 1L || !is.finite(n_perm) || n_perm < 0 || n_perm != as.integer(n_perm)) {
    stop("n_perm must be one non-negative integer.", call. = FALSE)
  }
  n_perm <- as.integer(n_perm)
  values <- as.numeric(values)
  if (is.null(neighbors)) {
    if (is.null(x) || is.null(y)) stop("x and y are required when neighbors is NULL.", call. = FALSE)
    neighbors <- build_spatial_neighbors(
      x, y, method = neighbor_method, distance_threshold = distance_threshold,
      weights = weights, symmetric = symmetric
    )
  }
  if (!inherits(neighbors, "spatial_neighbors")) {
    stop("neighbors must come from build_spatial_neighbors().", call. = FALSE)
  }
  n_total <- length(values)
  if (n_total != length(neighbors$x)) {
    stop("values must contain one value per spatial-neighbor node.", call. = FALSE)
  }
  edges <- neighbors$undirected_edges
  finite <- is.finite(values)
  if (nrow(edges)) {
    edges <- edges[finite[edges$from] & finite[edges$to], , drop = FALSE]
  }
  if (any(!is.finite(edges$weight)) || any(edges$weight < 0)) {
    stop("Spatial neighbor weights must be finite and non-negative.", call. = FALSE)
  }
  degree <- tabulate(c(edges$from, edges$to), nbins = n_total)
  effective <- finite & degree > 0L
  n_effective <- sum(effective)
  n_isolated <- n_total - n_effective
  n_edges <- nrow(edges)
  if (!n_edges || sum(edges$weight) <= 0) {
    stop("Geary's C is undefined because the effective graph has no positive-weight edges.", call. = FALSE)
  }
  if (n_effective < 3L) {
    stop("Geary's C requires at least three non-isolated finite nodes.", call. = FALSE)
  }
  index_map <- integer(n_total)
  index_map[effective] <- seq_len(n_effective)
  edge_i <- index_map[edges$from]
  edge_j <- index_map[edges$to]
  edge_weight <- edges$weight
  effective_values <- values[effective]

  geary_stat <- function(v) {
    centered <- v - mean(v)
    denominator <- sum(centered^2)
    if (!is.finite(denominator) || denominator == 0) return(NA_real_)
    numerator <- sum(edge_weight * (v[edge_i] - v[edge_j])^2)
    ((n_effective - 1) / (2 * sum(edge_weight))) * numerator / denominator
  }
  observed <- geary_stat(effective_values)
  p_value <- NA_real_
  if (n_perm > 0L && is.finite(observed)) {
    if (!is.null(seed)) set.seed(seed)
    permuted <- replicate(
      n_perm,
      geary_stat(sample(effective_values, n_effective, replace = FALSE))
    )
    finite_permuted <- permuted[is.finite(permuted)]
    if (!length(finite_permuted)) {
      p_value <- NA_real_
    } else if (alternative == "less") {
      p_value <- (sum(finite_permuted <= observed) + 1) / (length(finite_permuted) + 1)
    } else if (alternative == "greater") {
      p_value <- (sum(finite_permuted >= observed) + 1) / (length(finite_permuted) + 1)
    } else {
      # Under random labeling, the exact randomization expectation of global
      # Geary's C is one for a fixed symmetric graph.
      p_value <- (sum(abs(finite_permuted - 1) >= abs(observed - 1)) + 1) /
        (length(finite_permuted) + 1)
    }
  }
  list(
    C = observed, expected_C = 1, p_value = p_value,
    alternative = alternative, n = n_effective, n_total = n_total,
    n_effective = n_effective, n_isolated = n_isolated, n_edges = n_edges,
    weight_sum_undirected = sum(edge_weight)
  )
}

compare_contamination_methods <- function(pixel_matrix,
                                          fcols = feature_columns(pixel_matrix),
                                          x_col = "x",
                                          y_col = "y",
                                          foreground_mask = NULL,
                                          gap_min_multiplier = 6,
                                          mad_z_cutoff = 3.5,
                                          enrichment_exclude_if = c("background_ge_tissue", "background_gt_tissue"),
                                          include_moran = FALSE,
                                          n_perm = 199,
                                          moran_alpha = 0.05,
                                          agreement_methods = c("gap", "mad", "enrichment")) {
  enrichment_exclude_if <- match.arg(enrichment_exclude_if)
  required_columns(pixel_matrix, c(x_col, y_col), "Pixel matrix")
  if (length(fcols) == 0) {
    return(list(summary = data.frame(), agreement = data.frame()))
  }

  detection <- detect_ubiquitous_features(
    pixel_matrix,
    fcols = fcols,
    method = "consensus",
    gap_min_multiplier = gap_min_multiplier,
    mad_z_cutoff = mad_z_cutoff,
    min_agreement = 2
  )

  if (is.null(foreground_mask)) {
    raw_tic <- rowSums(pixel_matrix[fcols], na.rm = TRUE)
    log_tic <- log10(raw_tic + 1)
    threshold <- otsu_threshold(log_tic)
    foreground_mask <- log_tic >= threshold
  }
  foreground_mask <- as.logical(foreground_mask)
  if (length(foreground_mask) != nrow(pixel_matrix)) {
    stop("foreground_mask must have one value per pixel.", call. = FALSE)
  }
  background_mask <- !foreground_mask

  tissue_mean <- vapply(fcols, function(feature) {
    mean(pixel_matrix[[feature]][foreground_mask], na.rm = TRUE)
  }, numeric(1))
  background_mean <- vapply(fcols, function(feature) {
    mean(pixel_matrix[[feature]][background_mask], na.rm = TRUE)
  }, numeric(1))
  enrichment_ratio <- background_mean / pmax(tissue_mean, .Machine$double.eps)
  enrichment_outlier <- if (enrichment_exclude_if == "background_ge_tissue") {
    enrichment_ratio >= 1
  } else {
    enrichment_ratio > 1
  }

  moran_i <- rep(NA_real_, length(fcols))
  moran_p <- rep(NA_real_, length(fcols))
  names(moran_i) <- fcols
  names(moran_p) <- fcols
  if (isTRUE(include_moran)) {
    for (feature in fcols) {
      stat <- compute_morans_i_grid(
        values = pixel_matrix[[feature]],
        x = pixel_matrix[[x_col]],
        y = pixel_matrix[[y_col]],
        n_perm = n_perm,
        alternative = "greater"
      )
      moran_i[[feature]] <- stat$I
      moran_p[[feature]] <- stat$p_value
    }
  }
  moran_low_structure <- is.finite(moran_p) & moran_p > moran_alpha

  summary <- detection$agreement
  summary$tissue_mean <- unname(tissue_mean[summary$column_name])
  summary$background_mean <- unname(background_mean[summary$column_name])
  summary$background_to_tissue_ratio <- unname(enrichment_ratio[summary$column_name])
  summary$enrichment_outlier <- unname(enrichment_outlier[summary$column_name])
  summary$moran_i <- unname(moran_i[summary$column_name])
  summary$moran_p <- unname(moran_p[summary$column_name])
  summary$moran_low_structure <- unname(moran_low_structure[summary$column_name])

  available_methods <- intersect(agreement_methods, c("gap", "mad", "enrichment", "moran"))
  method_matrix <- data.frame(
    gap = summary$gap_outlier,
    mad = summary$mad_outlier,
    enrichment = summary$enrichment_outlier,
    moran = summary$moran_low_structure
  )
  summary$contaminant_agreement <- rowSums(method_matrix[available_methods], na.rm = TRUE)

  agreement <- summary[summary$contaminant_agreement > 0, , drop = FALSE]
  agreement <- agreement[order(-agreement$contaminant_agreement, -agreement$nonzero_rate), , drop = FALSE]

  list(
    summary = summary,
    agreement = agreement,
    gap = detection$gap,
    mad = detection$mad,
    foreground_mask = foreground_mask,
    agreement_methods = available_methods
  )
}

compute_foreground_tic <- function(pixel_matrix,
                                   fcols = feature_columns(pixel_matrix),
                                   method = c("raw", "exclude_features", "trimmed", "median"),
                                   excluded_features = character(),
                                   trim = 0.1) {
  method <- match.arg(method)
  if (length(fcols) == 0) return(numeric(nrow(pixel_matrix)))

  if (method == "exclude_features") {
    tic_fcols <- setdiff(fcols, excluded_features)
    if (length(tic_fcols) == 0) tic_fcols <- fcols
    return(rowSums(pixel_matrix[tic_fcols], na.rm = TRUE))
  }

  values <- as.matrix(pixel_matrix[fcols])
  if (method == "trimmed") {
    return(apply(values, 1, mean, trim = trim, na.rm = TRUE) * ncol(values))
  }
  if (method == "median") {
    return(apply(values, 1, stats::median, na.rm = TRUE) * ncol(values))
  }

  rowSums(values, na.rm = TRUE)
}

cleanup_foreground_components <- function(pixel_matrix,
                                          keep_rows,
                                          x_col = "x",
                                          y_col = "y",
                                          method = c("none", "largest_component", "area_fraction"),
                                          connectivity = c("8", "4"),
                                          min_component_fraction = 0.05) {
  method <- match.arg(method)
  connectivity <- match.arg(connectivity)
  keep_rows <- as.logical(keep_rows)
  if (length(keep_rows) != nrow(pixel_matrix)) {
    stop("keep_rows must have one value per pixel.", call. = FALSE)
  }

  stats <- list(
    cleanup_method = method,
    cleanup_connectivity = connectivity,
    foreground_components = 0L,
    largest_component_pixels = sum(keep_rows, na.rm = TRUE),
    pixels_after_threshold = sum(keep_rows, na.rm = TRUE),
    pixels_after_cleanup = sum(keep_rows, na.rm = TRUE),
    pixels_removed_by_cleanup = 0L,
    min_component_fraction = min_component_fraction
  )

  if (method == "none" || sum(keep_rows, na.rm = TRUE) == 0) {
    return(list(keep_rows = keep_rows, stats = stats))
  }

  required_columns(pixel_matrix, c(x_col, y_col), "Pixel matrix")
  foreground_idx <- which(keep_rows)
  x <- pixel_matrix[[x_col]][foreground_idx]
  y <- pixel_matrix[[y_col]][foreground_idx]
  valid <- is.finite(x) & is.finite(y)
  if (!all(valid)) {
    foreground_idx <- foreground_idx[valid]
    x <- x[valid]
    y <- y[valid]
  }
  if (length(foreground_idx) == 0) {
    keep_rows[] <- FALSE
    stats$pixels_after_cleanup <- 0L
    stats$pixels_removed_by_cleanup <- stats$pixels_after_threshold
    return(list(keep_rows = keep_rows, stats = stats))
  }

  keys <- paste(x, y, sep = "\r")
  key_to_pos <- stats::setNames(seq_along(keys), keys)
  offsets <- if (connectivity == "8") {
    expand.grid(dx = -1:1, dy = -1:1)
  } else {
    data.frame(dx = c(-1, 1, 0, 0), dy = c(0, 0, -1, 1))
  }
  offsets <- offsets[!(offsets$dx == 0 & offsets$dy == 0), , drop = FALSE]

  visited <- rep(FALSE, length(foreground_idx))
  component_id <- integer(length(foreground_idx))
  component_sizes <- integer()
  current_component <- 0L

  for (start in seq_along(foreground_idx)) {
    if (visited[start]) next
    current_component <- current_component + 1L
    queue <- integer(length(foreground_idx))
    queue[[1]] <- start
    tail <- 1L
    visited[start] <- TRUE
    component_id[start] <- current_component
    head <- 1L
    size <- 0L

    while (head <= tail) {
      pos <- queue[[head]]
      head <- head + 1L
      size <- size + 1L

      neighbor_keys <- paste(x[pos] + offsets$dx, y[pos] + offsets$dy, sep = "\r")
      neighbor_pos <- unname(key_to_pos[neighbor_keys])
      neighbor_pos <- neighbor_pos[!is.na(neighbor_pos)]
      for (next_pos in neighbor_pos) {
        if (!visited[next_pos]) {
          visited[next_pos] <- TRUE
          component_id[next_pos] <- current_component
          tail <- tail + 1L
          queue[[tail]] <- next_pos
        }
      }
    }
    component_sizes[current_component] <- size
  }

  stats$foreground_components <- length(component_sizes)
  stats$largest_component_pixels <- if (length(component_sizes) > 0) max(component_sizes) else 0L

  if (length(component_sizes) == 0) {
    keep_rows[] <- FALSE
  } else {
    if (method == "largest_component") {
      keep_components <- which(component_sizes == max(component_sizes))[1]
    } else {
      keep_components <- which(component_sizes >= max(component_sizes) * min_component_fraction)
    }
    cleaned_foreground_idx <- foreground_idx[component_id %in% keep_components]
    keep_rows[] <- FALSE
    keep_rows[cleaned_foreground_idx] <- TRUE
  }

  stats$pixels_after_cleanup <- sum(keep_rows, na.rm = TRUE)
  stats$pixels_removed_by_cleanup <- stats$pixels_after_threshold - stats$pixels_after_cleanup
  list(keep_rows = keep_rows, stats = stats)
}

make_feature_mapping <- function(mzmine_features) {
  required_columns(mzmine_features, c("id", "mz"), "MZmine feature list")
  mz_values <- as.numeric(mzmine_features$mz)
  if (any(!is.finite(mz_values))) {
    stop("MZmine feature list contains non-numeric mz values.", call. = FALSE)
  }
  column_names <- paste0("mz_", format(mz_values, trim = TRUE, scientific = FALSE, digits = 15))
  if (any(duplicated(column_names))) {
    duplicated_names <- unique(column_names[duplicated(column_names)])
    stop(
      "Duplicate mz column names generated: ",
      paste(duplicated_names, collapse = ", "),
      ". Check for duplicate or insufficiently precise mz values.",
      call. = FALSE
    )
  }

  data.frame(
    feature_id = seq_along(mz_values),
    mzmine_id = mzmine_features$id,
    mz = mz_values,
    column_name = column_names,
    stringsAsFactors = FALSE
  )
}

nearest_intensity_ppm <- function(spectrum_mz, spectrum_intensity, targets, ppm = 10) {
  out <- numeric(length(targets))
  if (length(spectrum_mz) == 0 || length(spectrum_intensity) == 0) return(out)

  spectrum_mz <- as.numeric(spectrum_mz)
  spectrum_intensity <- as.numeric(spectrum_intensity)
  keep <- is.finite(spectrum_mz) & is.finite(spectrum_intensity)
  spectrum_mz <- spectrum_mz[keep]
  spectrum_intensity <- spectrum_intensity[keep]
  if (length(spectrum_mz) == 0) return(out)

  order_mz <- order(spectrum_mz)
  spectrum_mz <- spectrum_mz[order_mz]
  spectrum_intensity <- spectrum_intensity[order_mz]
  insertion <- findInterval(targets, spectrum_mz)
  left <- pmax(1L, pmin(length(spectrum_mz), insertion))
  right <- pmax(1L, pmin(length(spectrum_mz), insertion + 1L))
  left_distance <- abs(spectrum_mz[left] - targets)
  right_distance <- abs(spectrum_mz[right] - targets)
  closest <- ifelse(right_distance < left_distance, right, left)
  within <- abs(spectrum_mz[closest] - targets) <= targets * ppm / 1e6
  out[within] <- spectrum_intensity[closest[within]]

  out
}

as_numeric_vector <- function(x) {
  if (isS4(x)) {
    x <- as.vector(x)
  }
  if (is.list(x)) {
    x <- unlist(x, use.names = FALSE)
  }
  as.numeric(x)
}

is_list_like_spectra <- function(x) {
  is.list(x) || inherits(x, "matter_list")
}

spectra_element <- function(x, i) {
  x[[i]]
}

extract_target_intensity_list <- function(spectrum_mz, intensity_data, targets, ppm, n_pixels,
                                          progress_callback = NULL) {
  out <- matrix(0, nrow = n_pixels, ncol = length(targets))
  shared_mz <- NULL

  if (!is_list_like_spectra(spectrum_mz)) {
    shared_mz <- as_numeric_vector(spectrum_mz)
  }

  for (pixel_index in seq_len(n_pixels)) {
    mz_values <- if (is.null(shared_mz)) {
      as_numeric_vector(spectra_element(spectrum_mz, pixel_index))
    } else {
      shared_mz
    }
    intensity_values <- as_numeric_vector(spectra_element(intensity_data, pixel_index))
    out[pixel_index, ] <- nearest_intensity_ppm(
      spectrum_mz = mz_values,
      spectrum_intensity = intensity_values,
      targets = targets,
      ppm = ppm
    )
    if (!is.null(progress_callback) &&
        (pixel_index == 1L || pixel_index == n_pixels || pixel_index %% max(1L, floor(n_pixels / 100)) == 0L)) {
      progress_callback(pixel_index / n_pixels, paste0("Extracting pixel ", pixel_index, " / ", n_pixels))
    }
  }

  out
}

extract_target_intensity_matrix <- function(msi, targets, ppm = 10, progress_callback = NULL) {
  raw_mz <- Cardinal::mz(msi)
  intensity_data <- Cardinal::intensity(msi)
  n_pixels <- nrow(as.data.frame(Cardinal::coord(msi)))

  if (is_list_like_spectra(intensity_data) || is.null(dim(intensity_data))) {
    return(extract_target_intensity_list(
      raw_mz, intensity_data, targets, ppm, n_pixels,
      progress_callback = progress_callback
    ))
  }

  spectrum_mz <- as_numeric_vector(raw_mz)
  feature_by_pixel <- dim(intensity_data)[2] == n_pixels
  pixel_by_feature <- dim(intensity_data)[1] == n_pixels
  if (feature_by_pixel && pixel_by_feature) {
    warning(
      "Ambiguous Cardinal intensity matrix orientation: both dimensions equal the number of pixels. ",
      "Assuming features x pixels.",
      call. = FALSE
    )
  }
  if (!feature_by_pixel && !pixel_by_feature) {
    stop(
      "Cannot determine Cardinal intensity matrix orientation. Dimensions are ",
      paste(dim(intensity_data), collapse = " x "),
      " for ",
      n_pixels,
      " pixels.",
      call. = FALSE
    )
  }

  out <- matrix(0, nrow = n_pixels, ncol = length(targets))
  for (i in seq_along(targets)) {
    tolerance <- targets[i] * ppm / 1e6
    hits <- which(abs(spectrum_mz - targets[i]) <= tolerance)
    if (length(hits) == 0) {
      if (!is.null(progress_callback)) {
        progress_callback(i / length(targets), paste0("Extracting target ", i, " / ", length(targets)))
      }
      next
    }

    closest <- hits[which.min(abs(spectrum_mz[hits] - targets[i]))]
    values <- if (feature_by_pixel) {
      as_numeric_vector(intensity_data[closest, ])
    } else {
      as_numeric_vector(intensity_data[, closest])
    }
    out[, i] <- values
    if (!is.null(progress_callback)) {
      progress_callback(i / length(targets), paste0("Extracting target ", i, " / ", length(targets)))
    }
  }

  out
}

load_msi_target_features <- function(imzml_path,
                                     mzmine_features,
                                     ppm = 10,
                                     section_id = "Section01",
                                     section_order = 1,
                                     progress_callback = NULL) {
  if (!requireNamespace("Cardinal", quietly = TRUE)) {
    stop("Package 'Cardinal' is required for imzML loading.", call. = FALSE)
  }

  mapping <- make_feature_mapping(mzmine_features)
  if (!is.null(progress_callback)) progress_callback(0.02, "Opening imzML metadata")
  msi <- Cardinal::readMSIData(imzml_path)
  coordinates <- as.data.frame(Cardinal::coord(msi))
  required_columns(coordinates, c("x", "y"), "Cardinal coordinates")

  n_pixels <- nrow(coordinates)
  if (!is.null(progress_callback)) progress_callback(0.08, paste0("Loaded ", n_pixels, " pixel coordinates"))
  intensity_matrix <- extract_target_intensity_matrix(
    msi,
    mapping$mz,
    ppm = ppm,
    progress_callback = if (is.null(progress_callback)) NULL else function(value, detail) {
      progress_callback(0.08 + 0.9 * value, detail)
    }
  )
  colnames(intensity_matrix) <- mapping$column_name

  pixel_matrix <- data.frame(
    pixel_id = seq_len(n_pixels),
    local_pixel_id = seq_len(n_pixels),
    section_id = section_id,
    section_order = section_order,
    x = coordinates$x,
    y = coordinates$y,
    intensity_matrix,
    check.names = FALSE
  )

  if (!is.null(progress_callback)) progress_callback(1, "Extraction complete")

  list(pixel_matrix = pixel_matrix, feature_mapping = mapping)
}

stable_mz_column_names <- function(mz_values) {
  labels <- vapply(mz_values, function(value) {
    format(value, digits = 17, scientific = FALSE, trim = TRUE)
  }, character(1))
  make.unique(paste0("mz_", labels), sep = "_feature_")
}

orient_shared_intensity_matrix <- function(intensity_data, n_features, n_spectra) {
  intensity_dim <- dim(intensity_data)
  if (length(intensity_dim) != 2L) {
    stop("Shared-axis centroided intensity data must be a two-dimensional matrix.", call. = FALSE)
  }
  feature_by_spectrum <- identical(as.integer(intensity_dim), c(as.integer(n_features), as.integer(n_spectra)))
  spectrum_by_feature <- identical(as.integer(intensity_dim), c(as.integer(n_spectra), as.integer(n_features)))
  if (feature_by_spectrum && spectrum_by_feature) {
    warning(
      "Ambiguous square intensity matrix; assuming features x spectra.",
      call. = FALSE
    )
    spectrum_by_feature <- FALSE
  }
  if (!feature_by_spectrum && !spectrum_by_feature) {
    stop(
      "Intensity dimensions ", paste(intensity_dim, collapse = " x "),
      " do not match ", n_features, " shared m/z features and ", n_spectra, " spectra.",
      call. = FALSE
    )
  }
  values <- as.matrix(intensity_data)
  if (feature_by_spectrum) {
    values <- t(values)
    orientation <- "features_x_spectra"
  } else {
    orientation <- "spectra_x_features"
  }
  storage.mode(values) <- "double"
  list(values = values, source_orientation = orientation)
}

load_centroided_msi_features <- function(imzml_path,
                                         sample_id = NULL,
                                         section_id = NULL,
                                         ion_mode = NULL,
                                         ion_mode_source = NULL) {
  if (!is.character(imzml_path) || length(imzml_path) != 1L || !nzchar(imzml_path)) {
    stop("imzml_path must be one non-empty path.", call. = FALSE)
  }
  imzml_path <- normalizePath(imzml_path, mustWork = TRUE)
  extension <- tools::file_ext(imzml_path)
  if (tolower(extension) != "imzml") {
    stop("imzml_path must point to an .imzML file.", call. = FALSE)
  }
  ibd_path <- paste0(substr(imzml_path, 1L, nchar(imzml_path) - nchar(extension)), "ibd")
  if (!file.exists(ibd_path) || !isTRUE(file.info(ibd_path)$isdir == FALSE)) {
    stop("The companion .ibd file does not exist: ", ibd_path, call. = FALSE)
  }
  validate_scalar_text <- function(value, label) {
    if (missing(value) || is.null(value) || length(value) != 1L || is.na(value) || !nzchar(trimws(value))) {
      stop(label, " must be supplied explicitly as one non-empty value.", call. = FALSE)
    }
    as.character(value)
  }
  sample_id <- validate_scalar_text(sample_id, "sample_id")
  section_id <- validate_scalar_text(section_id, "section_id")
  ion_mode <- validate_scalar_text(ion_mode, "ion_mode")
  ion_mode_source <- validate_scalar_text(ion_mode_source, "ion_mode_source")
  if (!ion_mode %in% c("positive", "negative")) {
    stop("ion_mode must be either 'positive' or 'negative'.", call. = FALSE)
  }
  if (!requireNamespace("Cardinal", quietly = TRUE)) {
    stop("Package 'Cardinal' is required for imzML loading.", call. = FALSE)
  }

  msi <- Cardinal::readMSIData(imzml_path)
  experiment <- Cardinal::experimentData(msi)
  spectrum_type <- as.character(experiment$spectrumType)
  if (length(spectrum_type) != 1L || !identical(spectrum_type, "MS1 spectrum")) {
    stop(
      "Only MS1 imzML data are supported; observed spectrum type: ",
      paste(spectrum_type, collapse = ", "),
      call. = FALSE
    )
  }
  if (!isTRUE(Cardinal::isCentroided(msi))) {
    stop("The imzML data must be centroided.", call. = FALSE)
  }

  coordinates <- as.data.frame(Cardinal::coord(msi))
  required_columns(coordinates, c("x", "y"), "Cardinal coordinates")
  if (any(!is.finite(coordinates$x) | !is.finite(coordinates$y))) {
    stop("MSI coordinates must be finite and non-missing.", call. = FALSE)
  }
  if (any(duplicated(coordinates[, c("x", "y"), drop = FALSE]))) {
    stop("MSI coordinates must be unique within a sample.", call. = FALSE)
  }

  mz_data <- Cardinal::mz(msi)
  if (is_list_like_spectra(mz_data) || !is.null(dim(mz_data))) {
    stop("The imzML data do not have one shared m/z axis.", call. = FALSE)
  }
  mz_values <- as.numeric(mz_data)
  if (!length(mz_values) || any(!is.finite(mz_values))) {
    stop("The shared m/z axis is empty or contains non-finite values.", call. = FALSE)
  }
  oriented <- orient_shared_intensity_matrix(
    Cardinal::intensity(msi),
    n_features = length(mz_values),
    n_spectra = nrow(coordinates)
  )
  column_names <- stable_mz_column_names(mz_values)
  if (anyDuplicated(column_names)) {
    stop("Stable m/z column names are not unique.", call. = FALSE)
  }
  colnames(oriented$values) <- column_names

  pixel_matrix <- data.frame(
    pixel_id = seq_len(nrow(coordinates)),
    sample_id = rep(sample_id, nrow(coordinates)),
    section_id = rep(section_id, nrow(coordinates)),
    x = coordinates$x,
    y = coordinates$y,
    oriented$values,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  feature_metadata <- data.frame(
    feature_id = sprintf("feature_%04d", seq_along(mz_values)),
    column_name = column_names,
    mz = mz_values,
    ion_mode = rep(ion_mode, length(mz_values)),
    stringsAsFactors = FALSE
  )
  coordinate_table <- pixel_matrix[, c("pixel_id", "sample_id", "section_id", "x", "y"), drop = FALSE]
  qc_summary <- list(
    spectra_count = nrow(coordinates),
    feature_count = length(mz_values),
    spectrum_type = spectrum_type,
    spectrum_representation = "centroided",
    shared_mz_axis = TRUE,
    intensity_source_orientation = oriented$source_orientation,
    intensity_output_orientation = "spectra_x_features",
    duplicate_coordinates = 0L,
    missing_coordinates = 0L
  )
  parameters <- list(
    input_type = "centroided_imzML_shared_mz",
    sample_id = sample_id,
    section_id = section_id,
    ion_mode = ion_mode,
    ion_mode_source = ion_mode_source,
    polarity_confirmed_by_imzml_cv = FALSE
  )
  provenance <- make_pipeline_manifest(
    input_files = c(imzml_path, ibd_path),
    input_type = parameters$input_type,
    ion_mode = ion_mode,
    parameters = parameters
  )
  provenance <- rbind(
    provenance,
    data.frame(
      record_type = "polarity_statement",
      key = c("polarity_interpretation_source", "polarity_cv_metadata_status"),
      value = c(ion_mode_source, "not_confirmed"),
      size_bytes = NA_real_,
      md5 = NA_character_,
      stringsAsFactors = FALSE
    )
  )

  list(
    pixel_feature_matrix = pixel_matrix,
    coordinates = coordinate_table,
    feature_metadata = feature_metadata,
    qc_summary = qc_summary,
    parameters = parameters,
    provenance = provenance
  )
}

prepare_lipid_annotation_mapping <- function(annotation_path,
                                             feature_metadata,
                                             ppm) {
  if (!file.exists(annotation_path)) {
    stop("Annotation file does not exist: ", annotation_path, call. = FALSE)
  }
  required_columns(feature_metadata, c("feature_id", "column_name", "mz"), "Feature metadata")
  if (length(ppm) != 1L || !is.finite(ppm) || ppm <= 0) {
    stop("ppm must be one positive finite value.", call. = FALSE)
  }
  annotation <- utils::read.delim(
    annotation_path,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  annotation_columns <- c("Input Mass", "Matched Mass", "Delta", "Name", "Formula", "Ion", "LMSD Examples")
  required_columns(annotation, annotation_columns, "Lipid annotation table")
  annotation <- annotation[, annotation_columns, drop = FALSE]
  annotation[["Input Mass"]] <- suppressWarnings(as.numeric(annotation[["Input Mass"]]))
  if (any(!is.finite(annotation[["Input Mass"]]))) {
    stop("Input Mass must contain only finite numeric values.", call. = FALSE)
  }
  feature_mz <- suppressWarnings(as.numeric(feature_metadata$mz))
  if (any(!is.finite(feature_mz))) {
    stop("Feature metadata m/z values must be finite numeric values.", call. = FALSE)
  }

  nearest <- vapply(annotation[["Input Mass"]], function(query_mz) {
    which.min(abs(feature_mz - query_mz))
  }, integer(1))
  signed_ppm_error <- (feature_mz[nearest] - annotation[["Input Mass"]]) /
    annotation[["Input Mass"]] * 1e6
  matched <- abs(signed_ppm_error) <= ppm + (.Machine$double.eps * 1e6)
  annotation$msi_feature_id <- ifelse(matched, feature_metadata$feature_id[nearest], NA_character_)
  annotation$msi_column_name <- ifelse(matched, feature_metadata$column_name[nearest], NA_character_)
  annotation$msi_mz <- ifelse(matched, feature_mz[nearest], NA_real_)
  annotation$ppm_error <- ifelse(matched, abs(signed_ppm_error), NA_real_)
  annotation$match_status <- ifelse(matched, "matched", "unmatched")
  annotation$annotation_evidence <- "putative_accurate_mass_annotation"
  annotation
}

map_spatial_domain_labels <- function(pixel_matrix,
                                      umap_path,
                                      sample_id,
                                      sample_column = "sample",
                                      label_column = "label") {
  required_columns(pixel_matrix, c("sample_id", "x", "y"), "Pixel matrix")
  if (!file.exists(umap_path)) {
    stop("UMAP/domain file does not exist: ", umap_path, call. = FALSE)
  }
  if (length(sample_id) != 1L || is.na(sample_id) || !nzchar(sample_id)) {
    stop("sample_id must be supplied explicitly.", call. = FALSE)
  }
  domains <- utils::read.csv(umap_path, check.names = FALSE, stringsAsFactors = FALSE)
  required_columns(domains, c(sample_column, "x", "y", label_column), "UMAP/domain table")
  domains <- domains[as.character(domains[[sample_column]]) == as.character(sample_id), , drop = FALSE]
  if (nrow(domains) == 0L) {
    stop("No domain rows were found for sample_id '", sample_id, "'.", call. = FALSE)
  }
  domain_keys <- paste(domains[[sample_column]], domains$x, domains$y, sep = "\r")
  if (anyDuplicated(domain_keys)) {
    stop("Domain labels are not one-to-one by sample + x + y.", call. = FALSE)
  }
  pixel_keys <- paste(pixel_matrix$sample_id, pixel_matrix$x, pixel_matrix$y, sep = "\r")
  if (anyDuplicated(pixel_keys)) {
    stop("Pixel coordinates are not one-to-one by sample + x + y.", call. = FALSE)
  }
  if (length(pixel_keys) != length(domain_keys) || !setequal(pixel_keys, domain_keys)) {
    stop(
      "Pixel and domain coordinates must match one-to-one for the requested sample.",
      call. = FALSE
    )
  }
  match_index <- match(pixel_keys, domain_keys)
  raw_label <- as.character(domains[[label_column]][match_index])
  allowed <- c("-1", "0", "1", "2")
  if (any(!raw_label %in% allowed)) {
    stop("Domain labels must be one of -1, 0, 1, or 2.", call. = FALSE)
  }
  out <- pixel_matrix
  out$domain_id <- raw_label
  out$domain_label <- ifelse(
    raw_label == "-1",
    "unclassified/background",
    paste0("metabolic_domain_", raw_label)
  )
  out$domain_type <- ifelse(
    raw_label == "-1",
    "unclassified/background",
    "data-driven metabolic domain"
  )
  out
}

load_serial_msi_target_features <- function(imzml_paths,
                                            mzmine_features,
                                            ppm = 10,
                                            section_ids = NULL,
                                            section_order = NULL,
                                            progress_callback = NULL) {
  if (length(imzml_paths) == 0) {
    stop("imzml_paths must contain at least one file.", call. = FALSE)
  }

  n_sections <- length(imzml_paths)
  if (is.null(section_ids)) {
    section_ids <- sprintf("Section%02d", seq_len(n_sections))
  }
  if (is.null(section_order)) {
    section_order <- seq_len(n_sections)
  }
  if (length(section_ids) != n_sections || length(section_order) != n_sections) {
    stop("section_ids and section_order must match length(imzml_paths).", call. = FALSE)
  }

  section_results <- lapply(seq_along(imzml_paths), function(i) {
    load_msi_target_features(
      imzml_path = imzml_paths[i],
      mzmine_features = mzmine_features,
      ppm = ppm,
      section_id = section_ids[i],
      section_order = section_order[i],
      progress_callback = if (is.null(progress_callback)) NULL else function(value, detail) {
        progress_callback(
          ((i - 1) + value) / n_sections,
          paste0("Section ", i, " / ", n_sections, ": ", detail)
        )
      }
    )
  })

  mapping <- section_results[[1]]$feature_mapping
  pixel_rows <- list()
  next_pixel_id <- 1
  for (i in seq_along(section_results)) {
    section_matrix <- section_results[[i]]$pixel_matrix
    section_matrix$local_pixel_id <- section_matrix$pixel_id
    n_pixels <- nrow(section_matrix)
    section_matrix$pixel_id <- seq.int(next_pixel_id, length.out = n_pixels)
    next_pixel_id <- next_pixel_id + n_pixels
    pixel_rows[[i]] <- section_matrix
  }

  combined <- do.call(rbind, pixel_rows)
  rownames(combined) <- NULL

  list(
    pixel_matrix = combined,
    feature_mapping = mapping,
    section_mapping = data.frame(
      section_id = section_ids,
      section_order = section_order,
      imzml_path = imzml_paths,
      stringsAsFactors = FALSE
    )
  )
}

feature_stats <- function(pixel_matrix) {
  fcols <- feature_columns(pixel_matrix)
  if (length(fcols) == 0) return(data.frame())

  stats <- lapply(fcols, function(column_name) {
    values <- pixel_matrix[[column_name]]
    data.frame(
      column_name = column_name,
      mz = sub("^mz_", "", column_name),
      mean = mean(values, na.rm = TRUE),
      sd = stats::sd(values, na.rm = TRUE),
      cv = safe_cv(values),
      nonzero_rate = mean(values > 0, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  })

  do.call(rbind, stats)
}

select_features <- function(pixel_matrix,
                            feature_mapping = NULL,
                            cv_top_percent = 70,
                            mean_min = 0,
                            nonzero_min = 0.3,
                            manual_columns = character(),
                            combine_mode = c("union", "intersection")) {
  combine_mode <- match.arg(combine_mode)
  stats <- feature_stats(pixel_matrix)
  if (nrow(stats) == 0) return(data.frame())

  cutoff_prob <- max(0, min(1, 1 - cv_top_percent / 100))
  cv_cutoff <- as.numeric(stats::quantile(stats$cv, probs = cutoff_prob, na.rm = TRUE))
  rule_columns <- stats$column_name[
    stats$cv >= cv_cutoff &
      stats$mean > mean_min &
      stats$nonzero_rate >= nonzero_min
  ]
  manual_columns <- intersect(manual_columns, stats$column_name)

  if (combine_mode == "union") {
    selected <- union(rule_columns, manual_columns)
  } else {
    selected <- intersect(rule_columns, manual_columns)
  }

  out <- data.frame(
    column_name = selected,
    selection_method = vapply(selected, function(column_name) {
      in_rule <- column_name %in% rule_columns
      in_manual <- column_name %in% manual_columns
      if (in_rule && in_manual) "both" else if (in_manual) "manual" else "rule"
    }, character(1)),
    stringsAsFactors = FALSE
  )

  if (!is.null(feature_mapping)) {
    out <- merge(feature_mapping, out, by = "column_name", all.y = TRUE, sort = FALSE)
    out <- out[match(selected, out$column_name), ]
  } else {
    out$feature_id <- seq_len(nrow(out))
    out$mz <- as.numeric(sub("^mz_", "", out$column_name))
    out <- out[, c("feature_id", "mz", "column_name", "selection_method")]
  }

  rownames(out) <- NULL
  out
}

reduce_pixel_matrix <- function(pixel_matrix, selected_features) {
  required_columns(pixel_matrix, c("pixel_id", "x", "y"), "Pixel feature matrix")
  all_feature_columns <- feature_columns(pixel_matrix)
  selected_columns <- intersect(selected_features$column_name, names(pixel_matrix))
  metadata_columns <- setdiff(names(pixel_matrix), all_feature_columns)
  pixel_matrix[, c(metadata_columns, selected_columns), drop = FALSE]
}

select_shared_features <- function(pixel_matrix,
                                   feature_mapping = NULL,
                                   section_column = "section_id",
                                   min_section_fraction = 1,
                                   cv_top_percent = 70,
                                   mean_min = 0,
                                   nonzero_min = 0.3,
                                   manual_columns = character(),
                                   combine_mode = c("union", "intersection")) {
  combine_mode <- match.arg(combine_mode)
  required_columns(pixel_matrix, section_column, "Serial pixel matrix")

  sections <- unique(pixel_matrix[[section_column]])
  if (length(sections) == 0) return(data.frame())

  selections <- lapply(sections, function(section) {
    section_matrix <- pixel_matrix[pixel_matrix[[section_column]] == section, , drop = FALSE]
    selected <- select_features(
      pixel_matrix = section_matrix,
      feature_mapping = feature_mapping,
      cv_top_percent = cv_top_percent,
      mean_min = mean_min,
      nonzero_min = nonzero_min,
      manual_columns = manual_columns,
      combine_mode = combine_mode
    )
    if (nrow(selected) == 0) {
      warning("Section ", section, " yielded no selected features and will be excluded.", call. = FALSE)
      return(NULL)
    }
    selected$section_id <- section
    selected
  })
  selections <- Filter(Negate(is.null), selections)
  if (length(selections) == 0) return(data.frame())
  selection_rows <- do.call(rbind, selections)
  if (is.null(selection_rows) || nrow(selection_rows) == 0) return(data.frame())

  section_counts <- stats::aggregate(
    selection_rows$section_id,
    by = list(column_name = selection_rows$column_name),
    FUN = function(x) length(unique(x))
  )
  names(section_counts)[2] <- "selected_section_count"
  section_counts$selected_section_fraction <- section_counts$selected_section_count / length(sections)
  keep <- section_counts$column_name[
    section_counts$selected_section_fraction >= max(0, min(1, min_section_fraction))
  ]

  out <- unique(selection_rows[selection_rows$column_name %in% keep, , drop = FALSE])
  out <- out[!duplicated(out$column_name), , drop = FALSE]
  out <- merge(out, section_counts, by = "column_name", all.x = TRUE, sort = FALSE)
  out <- out[match(keep, out$column_name), , drop = FALSE]
  out$section_id <- NULL
  rownames(out) <- NULL
  out
}

preprocess_matrix <- function(pixel_matrix,
                              do_background = TRUE,
                              background_method = c("log_tic_kmeans", "tic_percentile", "median_fraction", "gap_otsu_log_tic", "otsu_log_tic", "trimmed_otsu_log_tic", "median_otsu_log_tic"),
                              background_percent = 1,
                              background_percentile = 10,
                              min_nonzero_features = 1,
                              ubiquitous_method = c("gap", "mad", "consensus", "none"),
                              do_ubiquitous_exclusion = FALSE,
                              gap_min_multiplier = 6,
                              mad_z_cutoff = 3.5,
                              min_ubiquitous_agreement = 2,
                              otsu_bins = 256,
                              trimmed_tic_trim = 0.1,
                              kmeans_background_k = 2,
                              kmeans_background_seed = 42,
                              kmeans_background_nstart = 25,
                              kmeans_background_iter_max = 100,
                              kmeans_background_algorithm = "Lloyd",
                              foreground_cleanup = c("none", "largest_component", "area_fraction"),
                              cleanup_connectivity = c("8", "4"),
                              min_component_fraction = 0.05,
                              do_tic = TRUE,
                              do_log = TRUE,
                              do_scale = FALSE) {
  background_method <- match.arg(background_method)
  ubiquitous_method <- match.arg(ubiquitous_method)
  foreground_cleanup <- match.arg(foreground_cleanup)
  cleanup_connectivity <- match.arg(cleanup_connectivity)
  required_columns(pixel_matrix, c("pixel_id", "x", "y"), "Pixel matrix")
  fcols <- feature_columns(pixel_matrix)
  out <- pixel_matrix
  distributions <- list(raw = unlist(out[fcols], use.names = FALSE))
  background_stats <- data.frame(
    method = background_method,
    pixels_before = nrow(out),
    pixels_after = nrow(out),
    pixels_removed = 0,
    removed_fraction = 0,
    threshold = NA_real_,
    tissue_cluster = NA_integer_,
    cluster1_mean_log_tic = NA_real_,
    cluster2_mean_log_tic = NA_real_,
    foreground_tic_method = NA_character_,
    ubiquitous_method = ubiquitous_method,
    ubiquitous_threshold = NA_real_,
    ubiquitous_gap_size = NA_real_,
    ubiquitous_gap_cutoff = NA_real_,
    excluded_ubiquitous_features = 0L,
    removed_ubiquitous_features = 0L,
    background_tic_features = length(fcols),
    excluded_ubiquitous_columns = "",
    removed_ubiquitous_columns = "",
    otsu_bins = NA_integer_,
    kmeans_background_k = kmeans_background_k,
    kmeans_background_seed = kmeans_background_seed,
    kmeans_background_nstart = kmeans_background_nstart,
    kmeans_background_iter_max = kmeans_background_iter_max,
    kmeans_background_algorithm = kmeans_background_algorithm,
    foreground_cleanup = foreground_cleanup,
    cleanup_connectivity = cleanup_connectivity,
    foreground_components = NA_integer_,
    largest_component_pixels = NA_integer_,
    pixels_after_threshold = nrow(out),
    pixels_removed_by_cleanup = 0L,
    min_component_fraction = min_component_fraction,
    stringsAsFactors = FALSE
  )

  if (length(fcols) == 0) return(list(matrix = out, distributions = distributions, background_stats = background_stats))

  ubiquitous <- NULL
  excluded_ubiquitous <- character()

  if (do_background) {
    nonzero_features <- rowSums(out[fcols] > 0, na.rm = TRUE)
    threshold <- NA_real_
    row_tic <- rowSums(out[fcols], na.rm = TRUE)

    if (background_method %in% c("gap_otsu_log_tic", "otsu_log_tic")) {
      detector_method <- if (background_method == "gap_otsu_log_tic") ubiquitous_method else "none"
      if (background_method == "gap_otsu_log_tic" && detector_method == "none") detector_method <- "gap"
      ubiquitous <- detect_ubiquitous_features(
        out,
        fcols = fcols,
        method = detector_method,
        gap_min_multiplier = gap_min_multiplier,
        mad_z_cutoff = mad_z_cutoff,
        min_agreement = min_ubiquitous_agreement
      )
      excluded_ubiquitous <- ubiquitous$excluded
      if (length(setdiff(fcols, excluded_ubiquitous)) == 0) {
        excluded_ubiquitous <- character()
      }
      row_tic <- compute_foreground_tic(
        out,
        fcols = fcols,
        method = "exclude_features",
        excluded_features = excluded_ubiquitous
      )
      log_tic <- log10(row_tic + 1)
      otsu <- otsu_threshold(log_tic, nbins = otsu_bins, return_details = TRUE)
      threshold <- otsu$threshold
      keep_rows <- log_tic >= threshold & nonzero_features >= min_nonzero_features
      background_stats$foreground_tic_method <- if (length(excluded_ubiquitous) > 0) "clean_tic_excluding_ubiquitous" else "raw_tic"
      background_stats$ubiquitous_method <- detector_method
      background_stats$ubiquitous_threshold <- if (detector_method == "mad" && !is.null(ubiquitous$mad)) {
        ubiquitous$mad$threshold
      } else if (!is.null(ubiquitous$gap)) {
        ubiquitous$gap$threshold
      } else {
        NA_real_
      }
      background_stats$ubiquitous_gap_size <- if (!is.null(ubiquitous$gap)) ubiquitous$gap$gap_size else NA_real_
      background_stats$ubiquitous_gap_cutoff <- if (!is.null(ubiquitous$gap)) ubiquitous$gap$gap_cutoff else NA_real_
      background_stats$excluded_ubiquitous_features <- length(excluded_ubiquitous)
      background_stats$background_tic_features <- length(setdiff(fcols, excluded_ubiquitous))
      background_stats$excluded_ubiquitous_columns <- paste(excluded_ubiquitous, collapse = ";")
      background_stats$otsu_bins <- otsu$actual_bins
    } else if (background_method == "trimmed_otsu_log_tic") {
      row_tic <- compute_foreground_tic(out, fcols = fcols, method = "trimmed", trim = trimmed_tic_trim)
      log_tic <- log10(row_tic + 1)
      otsu <- otsu_threshold(log_tic, nbins = otsu_bins, return_details = TRUE)
      threshold <- otsu$threshold
      keep_rows <- log_tic >= threshold & nonzero_features >= min_nonzero_features
      background_stats$foreground_tic_method <- paste0("trimmed_tic_", trimmed_tic_trim)
      background_stats$otsu_bins <- otsu$actual_bins
    } else if (background_method == "median_otsu_log_tic") {
      row_tic <- compute_foreground_tic(out, fcols = fcols, method = "median")
      log_tic <- log10(row_tic + 1)
      otsu <- otsu_threshold(log_tic, nbins = otsu_bins, return_details = TRUE)
      threshold <- otsu$threshold
      keep_rows <- log_tic >= threshold & nonzero_features >= min_nonzero_features
      background_stats$foreground_tic_method <- "median_tic"
      background_stats$otsu_bins <- otsu$actual_bins
    } else if (background_method == "log_tic_kmeans" && length(unique(row_tic)) > 1) {
      kmeans_background_k <- max(2, min(as.integer(kmeans_background_k), length(unique(row_tic))))
      if (!is.null(kmeans_background_seed)) set.seed(kmeans_background_seed)
      log_tic <- log10(row_tic + 1)
      fit <- stats::kmeans(
        log_tic,
        centers = kmeans_background_k,
        nstart = kmeans_background_nstart,
        iter.max = kmeans_background_iter_max,
        algorithm = kmeans_background_algorithm
      )
      cluster_means <- tapply(log_tic, fit$cluster, mean, na.rm = TRUE)
      tissue_cluster <- as.integer(names(cluster_means)[which.max(cluster_means)])
      keep_rows <- fit$cluster == tissue_cluster & nonzero_features >= min_nonzero_features
      background_stats$tissue_cluster <- tissue_cluster
      background_stats$cluster1_mean_log_tic <- unname(cluster_means[["1"]])
      background_stats$cluster2_mean_log_tic <- unname(cluster_means[["2"]])
      if (kmeans_background_k != 2) {
        background_stats$cluster1_mean_log_tic <- NA_real_
        background_stats$cluster2_mean_log_tic <- NA_real_
      }
      background_stats$foreground_tic_method <- "raw_tic"
    } else {
      positive_tic <- row_tic[row_tic > 0]
      threshold <- if (background_method == "tic_percentile" && length(positive_tic) > 0) {
        stats::quantile(positive_tic, probs = background_percentile / 100, na.rm = TRUE)
      } else {
        stats::median(row_tic, na.rm = TRUE) * background_percent / 100
      }
      keep_rows <- row_tic >= threshold & nonzero_features >= min_nonzero_features
      background_stats$foreground_tic_method <- "raw_tic"
    }

    cleanup <- cleanup_foreground_components(
      out,
      keep_rows = keep_rows,
      x_col = "x",
      y_col = "y",
      method = foreground_cleanup,
      connectivity = cleanup_connectivity,
      min_component_fraction = min_component_fraction
    )
    keep_rows <- cleanup$keep_rows

    out <- out[keep_rows, , drop = FALSE]
    background_stats$pixels_after <- nrow(out)
    background_stats$pixels_removed <- background_stats$pixels_before - background_stats$pixels_after
    background_stats$removed_fraction <- background_stats$pixels_removed / background_stats$pixels_before
    background_stats$threshold <- threshold
    background_stats$foreground_components <- cleanup$stats$foreground_components
    background_stats$largest_component_pixels <- cleanup$stats$largest_component_pixels
    background_stats$pixels_after_threshold <- cleanup$stats$pixels_after_threshold
    background_stats$pixels_removed_by_cleanup <- cleanup$stats$pixels_removed_by_cleanup
  }

  if (isTRUE(do_ubiquitous_exclusion)) {
    if (is.null(ubiquitous)) {
      detector_method <- ubiquitous_method
      if (detector_method == "none") detector_method <- "gap"
      ubiquitous <- detect_ubiquitous_features(
        out,
        fcols = fcols,
        method = detector_method,
        gap_min_multiplier = gap_min_multiplier,
        mad_z_cutoff = mad_z_cutoff,
        min_agreement = min_ubiquitous_agreement
      )
      excluded_ubiquitous <- ubiquitous$excluded
      background_stats$ubiquitous_method <- detector_method
      background_stats$ubiquitous_threshold <- if (detector_method == "mad" && !is.null(ubiquitous$mad)) {
        ubiquitous$mad$threshold
      } else if (!is.null(ubiquitous$gap)) {
        ubiquitous$gap$threshold
      } else {
        NA_real_
      }
      background_stats$ubiquitous_gap_size <- if (!is.null(ubiquitous$gap)) ubiquitous$gap$gap_size else NA_real_
      background_stats$ubiquitous_gap_cutoff <- if (!is.null(ubiquitous$gap)) ubiquitous$gap$gap_cutoff else NA_real_
      background_stats$excluded_ubiquitous_features <- length(excluded_ubiquitous)
      background_stats$excluded_ubiquitous_columns <- paste(excluded_ubiquitous, collapse = ";")
    }

    removable <- intersect(excluded_ubiquitous, names(out))
    if (length(removable) > 0 && length(setdiff(fcols, removable)) > 0) {
      out <- out[, setdiff(names(out), removable), drop = FALSE]
      fcols <- setdiff(fcols, removable)
      background_stats$removed_ubiquitous_features <- length(removable)
      background_stats$removed_ubiquitous_columns <- paste(removable, collapse = ";")
      background_stats$background_tic_features <- length(fcols)
    }
  }

  distributions$background <- unlist(out[fcols], use.names = FALSE)
  if (nrow(out) == 0) {
    warning(
      "All pixels removed by background filtering. Check background_method and threshold.",
      call. = FALSE
    )
    return(list(matrix = out, distributions = distributions, background_stats = background_stats))
  }

  if (do_tic) {
    row_tic <- rowSums(out[fcols], na.rm = TRUE)
    scale_factor <- stats::median(row_tic[row_tic > 0], na.rm = TRUE)
    if (!is.finite(scale_factor)) scale_factor <- 1
    normalized <- sweep(as.matrix(out[fcols]), 1, ifelse(row_tic == 0, NA, row_tic), "/") * scale_factor
    normalized[!is.finite(normalized)] <- 0
    out[fcols] <- normalized
  }
  distributions$tic <- unlist(out[fcols], use.names = FALSE)

  if (do_log) {
    out[fcols] <- log10(as.matrix(out[fcols]) + 1)
  }
  distributions$log <- unlist(out[fcols], use.names = FALSE)

  if (do_scale) {
    scaled <- scale(as.matrix(out[fcols]))
    scaled[!is.finite(scaled)] <- 0
    out[fcols] <- scaled
  }
  distributions$scaled <- unlist(out[fcols], use.names = FALSE)

  rownames(out) <- NULL
  list(matrix = out, distributions = distributions, background_stats = background_stats)
}

filter_background_cardinal_spatial_kmeans <- function(pixel_matrix,
                                                      k = 2,
                                                      r = 1,
                                                      ncomp = NULL,
                                                      weights = c("adaptive", "gaussian"),
                                                      transform = c("log10", "none"),
                                                      tissue_cluster_method = c("low_border", "high_tic"),
                                                      foreground_cleanup = c("largest_component", "area_fraction", "none"),
                                                      cleanup_connectivity = c("8", "4"),
                                                      min_component_fraction = 0.05,
                                                      seed = 42,
                                                      verbose = FALSE) {
  if (!requireNamespace("Cardinal", quietly = TRUE)) {
    stop("Package 'Cardinal' is required for Cardinal spatial background filtering.", call. = FALSE)
  }
  required_columns(pixel_matrix, c("pixel_id", "x", "y"), "Pixel matrix")
  fcols <- feature_columns(pixel_matrix)
  if (length(fcols) == 0) {
    stop("No mz feature columns available for Cardinal spatial background filtering.", call. = FALSE)
  }

  weights <- match.arg(weights)
  transform <- match.arg(transform)
  tissue_cluster_method <- match.arg(tissue_cluster_method)
  foreground_cleanup <- match.arg(foreground_cleanup)
  cleanup_connectivity <- match.arg(cleanup_connectivity)
  k <- max(2, min(as.integer(k), nrow(pixel_matrix)))
  if (is.null(ncomp)) ncomp <- min(10L, length(fcols), nrow(pixel_matrix) - 1L)
  ncomp <- max(k, min(as.integer(ncomp), length(fcols), nrow(pixel_matrix) - 1L))

  matrix_data <- as.matrix(pixel_matrix[fcols])
  if (transform == "log10") {
    matrix_data[matrix_data < 0] <- 0
    matrix_data <- log10(matrix_data + 1)
  }
  matrix_data[!is.finite(matrix_data)] <- 0
  matrix_data <- unname(matrix_data)
  coord <- unname(as.matrix(pixel_matrix[, c("x", "y"), drop = FALSE]))

  if (!is.null(seed)) set.seed(seed)
  fit <- Cardinal::spatialKMeans(
    t(matrix_data),
    coord = coord,
    r = r,
    k = k,
    ncomp = ncomp,
    weights = weights,
    verbose = verbose
  )

  raw_tic <- rowSums(pixel_matrix[fcols], na.rm = TRUE)
  log_tic <- log10(raw_tic + 1)
  cluster_medians <- tapply(raw_tic, fit$cluster, stats::median, na.rm = TRUE)
  cluster_median_log_tic <- tapply(log_tic, fit$cluster, stats::median, na.rm = TRUE)
  border_pixel <- pixel_matrix$x %in% range(pixel_matrix$x, na.rm = TRUE) |
    pixel_matrix$y %in% range(pixel_matrix$y, na.rm = TRUE)
  cluster_border_fraction <- tapply(border_pixel, fit$cluster, mean, na.rm = TRUE)
  tissue_cluster <- if (tissue_cluster_method == "low_border") {
    candidates <- names(cluster_border_fraction)[
      cluster_border_fraction == min(cluster_border_fraction, na.rm = TRUE)
    ]
    if (length(candidates) > 1) {
      candidates[which.max(cluster_median_log_tic[candidates])]
    } else {
      candidates[[1]]
    }
  } else {
    names(cluster_medians)[which.max(cluster_medians)]
  }
  tissue_cluster <- as.integer(tissue_cluster)
  keep_rows <- fit$cluster == tissue_cluster

  cleanup <- cleanup_foreground_components(
    pixel_matrix,
    keep_rows = keep_rows,
    x_col = "x",
    y_col = "y",
    method = foreground_cleanup,
    connectivity = cleanup_connectivity,
    min_component_fraction = min_component_fraction
  )
  keep_rows <- cleanup$keep_rows
  out <- pixel_matrix[keep_rows, , drop = FALSE]
  rownames(out) <- NULL

  total_tic <- sum(raw_tic, na.rm = TRUE)
  retained_tic <- sum(raw_tic[keep_rows], na.rm = TRUE)
  background_stats <- data.frame(
    method = "cardinal_spatial_kmeans",
    pixels_before = nrow(pixel_matrix),
    pixels_after = nrow(out),
    pixels_removed = nrow(pixel_matrix) - nrow(out),
    removed_fraction = (nrow(pixel_matrix) - nrow(out)) / nrow(pixel_matrix),
    retained_tic_fraction = retained_tic / total_tic,
    tissue_cluster = tissue_cluster,
    tissue_cluster_method = tissue_cluster_method,
    cluster_median_log_tic = paste(
      names(cluster_median_log_tic),
      signif(cluster_median_log_tic, 4),
      sep = ":",
      collapse = ";"
    ),
    cluster_border_fraction = paste(
      names(cluster_border_fraction),
      signif(cluster_border_fraction, 4),
      sep = ":",
      collapse = ";"
    ),
    cardinal_k = k,
    cardinal_r = r,
    cardinal_ncomp = ncomp,
    cardinal_weights = weights,
    cardinal_transform = transform,
    foreground_cleanup = foreground_cleanup,
    cleanup_connectivity = cleanup_connectivity,
    foreground_components = cleanup$stats$foreground_components,
    largest_component_pixels = cleanup$stats$largest_component_pixels,
    pixels_after_threshold = cleanup$stats$pixels_after_threshold,
    pixels_removed_by_cleanup = cleanup$stats$pixels_removed_by_cleanup,
    min_component_fraction = min_component_fraction,
    stringsAsFactors = FALSE
  )

  list(
    matrix = out,
    keep = keep_rows,
    cluster = fit$cluster,
    fit = fit,
    background_stats = background_stats
  )
}

preprocess_sections <- function(pixel_matrix,
                                section_column = "section_id",
                                do_background = TRUE,
                                background_method = c("log_tic_kmeans", "tic_percentile", "median_fraction", "gap_otsu_log_tic", "otsu_log_tic", "trimmed_otsu_log_tic", "median_otsu_log_tic"),
                                background_percent = 1,
                                background_percentile = 10,
                                min_nonzero_features = 1,
                                ubiquitous_method = c("gap", "mad", "consensus", "none"),
                                do_ubiquitous_exclusion = FALSE,
                                gap_min_multiplier = 6,
                                mad_z_cutoff = 3.5,
                                min_ubiquitous_agreement = 2,
                                otsu_bins = 256,
                                trimmed_tic_trim = 0.1,
                                kmeans_background_k = 2,
                                kmeans_background_seed = 42,
                                kmeans_background_nstart = 25,
                                kmeans_background_iter_max = 100,
                                kmeans_background_algorithm = "Lloyd",
                                foreground_cleanup = c("none", "largest_component", "area_fraction"),
                                cleanup_connectivity = c("8", "4"),
                                min_component_fraction = 0.05,
                                do_tic = TRUE,
                                do_log = TRUE,
                                do_scale = FALSE,
                                scale_scope = c("global", "section")) {
  background_method <- match.arg(background_method)
  ubiquitous_method <- match.arg(ubiquitous_method)
  foreground_cleanup <- match.arg(foreground_cleanup)
  cleanup_connectivity <- match.arg(cleanup_connectivity)
  scale_scope <- match.arg(scale_scope)
  required_columns(pixel_matrix, section_column, "Serial pixel matrix")

  sections <- unique(pixel_matrix[[section_column]])
  results <- lapply(sections, function(section) {
    section_matrix <- pixel_matrix[pixel_matrix[[section_column]] == section, , drop = FALSE]
    result <- preprocess_matrix(
      section_matrix,
      do_background = do_background,
      background_method = background_method,
      background_percent = background_percent,
      background_percentile = background_percentile,
      min_nonzero_features = min_nonzero_features,
      ubiquitous_method = ubiquitous_method,
      do_ubiquitous_exclusion = do_ubiquitous_exclusion,
      gap_min_multiplier = gap_min_multiplier,
      mad_z_cutoff = mad_z_cutoff,
      min_ubiquitous_agreement = min_ubiquitous_agreement,
      otsu_bins = otsu_bins,
      trimmed_tic_trim = trimmed_tic_trim,
      kmeans_background_k = kmeans_background_k,
      kmeans_background_seed = kmeans_background_seed,
      kmeans_background_nstart = kmeans_background_nstart,
      kmeans_background_iter_max = kmeans_background_iter_max,
      kmeans_background_algorithm = kmeans_background_algorithm,
      foreground_cleanup = foreground_cleanup,
      cleanup_connectivity = cleanup_connectivity,
      min_component_fraction = min_component_fraction,
      do_tic = do_tic,
      do_log = do_log,
      do_scale = do_scale && identical(scale_scope, "section")
    )
    result$section_id <- section
    if (!is.null(result$background_stats)) {
      result$background_stats$section_id <- section
    }
    result
  })

  combined <- do.call(rbind, lapply(results, function(result) result$matrix))
  rownames(combined) <- NULL

  fcols <- feature_columns(combined)
  if (do_scale && identical(scale_scope, "global") && length(fcols) > 0) {
    scaled <- scale(as.matrix(combined[fcols]))
    scaled[!is.finite(scaled)] <- 0
    combined[fcols] <- scaled
  }

  distribution_names <- unique(unlist(lapply(results, function(result) names(result$distributions)), use.names = FALSE))
  distributions <- lapply(distribution_names, function(name) {
    unlist(lapply(results, function(result) result$distributions[[name]]), use.names = FALSE)
  })
  names(distributions) <- distribution_names
  if (do_scale && identical(scale_scope, "global") && length(fcols) > 0) {
    distributions$scaled <- unlist(combined[fcols], use.names = FALSE)
  }
  background_stats <- do.call(rbind, lapply(results, function(result) result$background_stats))
  rownames(background_stats) <- NULL

  list(
    matrix = combined,
    distributions = distributions,
    section_results = results,
    background_stats = background_stats,
    scale_scope = if (do_scale) scale_scope else "none"
  )
}

preprocess_select_features <- function(pixel_matrix,
                                       feature_mapping = NULL,
                                       serial = NULL,
                                       section_column = "section_id",
                                       min_section_fraction = 1,
                                       cv_top_percent = 70,
                                       mean_min = 0,
                                       nonzero_min = 0.3,
                                       manual_columns = character(),
                                       combine_mode = c("union", "intersection"),
                                       do_background = TRUE,
                                       background_method = c("log_tic_kmeans", "tic_percentile", "median_fraction", "gap_otsu_log_tic", "otsu_log_tic", "trimmed_otsu_log_tic", "median_otsu_log_tic"),
                                       background_percent = 1,
                                       background_percentile = 10,
                                       min_nonzero_features = 1,
                                       ubiquitous_method = c("gap", "mad", "consensus", "none"),
                                       do_ubiquitous_exclusion = FALSE,
                                       gap_min_multiplier = 6,
                                       mad_z_cutoff = 3.5,
                                       min_ubiquitous_agreement = 2,
                                       otsu_bins = 256,
                                       trimmed_tic_trim = 0.1,
                                       kmeans_background_k = 2,
                                       kmeans_background_seed = 42,
                                       kmeans_background_nstart = 25,
                                       kmeans_background_iter_max = 100,
                                       kmeans_background_algorithm = "Lloyd",
                                       foreground_cleanup = c("none", "largest_component", "area_fraction"),
                                       cleanup_connectivity = c("8", "4"),
                                       min_component_fraction = 0.05,
                                       do_tic = TRUE,
                                       do_log = TRUE) {
  combine_mode <- match.arg(combine_mode)
  background_method <- match.arg(background_method)
  ubiquitous_method <- match.arg(ubiquitous_method)
  foreground_cleanup <- match.arg(foreground_cleanup)
  cleanup_connectivity <- match.arg(cleanup_connectivity)
  if (is.null(serial)) {
    serial <- section_column %in% names(pixel_matrix) &&
      length(unique(pixel_matrix[[section_column]])) > 1
  }

  preprocessing <- if (isTRUE(serial)) {
    preprocess_sections(
      pixel_matrix,
      section_column = section_column,
      do_background = do_background,
      background_method = background_method,
      background_percent = background_percent,
      background_percentile = background_percentile,
      min_nonzero_features = min_nonzero_features,
      ubiquitous_method = ubiquitous_method,
      do_ubiquitous_exclusion = do_ubiquitous_exclusion,
      gap_min_multiplier = gap_min_multiplier,
      mad_z_cutoff = mad_z_cutoff,
      min_ubiquitous_agreement = min_ubiquitous_agreement,
      otsu_bins = otsu_bins,
      trimmed_tic_trim = trimmed_tic_trim,
      kmeans_background_k = kmeans_background_k,
      kmeans_background_seed = kmeans_background_seed,
      kmeans_background_nstart = kmeans_background_nstart,
      kmeans_background_iter_max = kmeans_background_iter_max,
      kmeans_background_algorithm = kmeans_background_algorithm,
      foreground_cleanup = foreground_cleanup,
      cleanup_connectivity = cleanup_connectivity,
      min_component_fraction = min_component_fraction,
      do_tic = do_tic,
      do_log = do_log,
      do_scale = FALSE,
      scale_scope = "global"
    )
  } else {
    preprocess_matrix(
      pixel_matrix,
      do_background = do_background,
      background_method = background_method,
      background_percent = background_percent,
      background_percentile = background_percentile,
      min_nonzero_features = min_nonzero_features,
      ubiquitous_method = ubiquitous_method,
      do_ubiquitous_exclusion = do_ubiquitous_exclusion,
      gap_min_multiplier = gap_min_multiplier,
      mad_z_cutoff = mad_z_cutoff,
      min_ubiquitous_agreement = min_ubiquitous_agreement,
      otsu_bins = otsu_bins,
      trimmed_tic_trim = trimmed_tic_trim,
      kmeans_background_k = kmeans_background_k,
      kmeans_background_seed = kmeans_background_seed,
      kmeans_background_nstart = kmeans_background_nstart,
      kmeans_background_iter_max = kmeans_background_iter_max,
      kmeans_background_algorithm = kmeans_background_algorithm,
      foreground_cleanup = foreground_cleanup,
      cleanup_connectivity = cleanup_connectivity,
      min_component_fraction = min_component_fraction,
      do_tic = do_tic,
      do_log = do_log,
      do_scale = FALSE
    )
  }

  selected <- if (isTRUE(serial)) {
    select_shared_features(
      preprocessing$matrix,
      feature_mapping = feature_mapping,
      section_column = section_column,
      min_section_fraction = min_section_fraction,
      cv_top_percent = cv_top_percent,
      mean_min = mean_min,
      nonzero_min = nonzero_min,
      manual_columns = manual_columns,
      combine_mode = combine_mode
    )
  } else {
    select_features(
      preprocessing$matrix,
      feature_mapping = feature_mapping,
      cv_top_percent = cv_top_percent,
      mean_min = mean_min,
      nonzero_min = nonzero_min,
      manual_columns = manual_columns,
      combine_mode = combine_mode
    )
  }

  reduced <- reduce_pixel_matrix(preprocessing$matrix, selected)
  list(
    normalized_matrix = preprocessing$matrix,
    selected_features = selected,
    reduced_matrix = reduced,
    background_stats = preprocessing$background_stats,
    distributions = preprocessing$distributions,
    preprocessing = preprocessing,
    serial = isTRUE(serial)
  )
}

apply_matched_region_labels <- function(pixel_matrix,
                                        region_labels,
                                        region_column = "matched_region_label") {
  required_columns(region_labels, region_column, "Region label table")
  out <- pixel_matrix

  if ("pixel_id" %in% names(region_labels)) {
    label_data <- region_labels[, c("pixel_id", region_column), drop = FALSE]
    out <- merge(out, label_data, by = "pixel_id", all.x = TRUE, sort = FALSE)
  } else if (all(c("section_id", "local_pixel_id") %in% names(region_labels)) &&
             all(c("section_id", "local_pixel_id") %in% names(pixel_matrix))) {
    label_data <- region_labels[, c("section_id", "local_pixel_id", region_column), drop = FALSE]
    out <- merge(out, label_data, by = c("section_id", "local_pixel_id"), all.x = TRUE, sort = FALSE)
  } else if (all(c("section_id", "x", "y") %in% names(region_labels)) &&
             all(c("section_id", "x", "y") %in% names(pixel_matrix))) {
    label_data <- region_labels[, c("section_id", "x", "y", region_column), drop = FALSE]
    duplicate_keys <- duplicated(label_data[, c("section_id", "x", "y"), drop = FALSE])
    if (any(duplicate_keys)) {
      warning(
        "Duplicate (section_id, x, y) entries in region_labels. Only the first match will be used.",
        call. = FALSE
      )
      label_data <- label_data[!duplicate_keys, , drop = FALSE]
    }
    out <- merge(out, label_data, by = c("section_id", "x", "y"), all.x = TRUE, sort = FALSE)
  } else {
    stop(
      "Region label table must contain matched_region_label plus pixel_id, ",
      "section_id + local_pixel_id, or section_id + x + y.",
      call. = FALSE
    )
  }

  out[order(out$pixel_id), , drop = FALSE]
}

cluster_pixels <- function(pixel_matrix,
                           k = 3,
                           nstart = 25,
                           iter.max = 100,
                           algorithm = "Lloyd",
                           seed = 1,
                           pca_components = NULL) {
  fcols <- feature_columns(pixel_matrix)
  if (length(fcols) == 0) stop("No mz feature columns available for clustering.", call. = FALSE)
  k <- max(2, min(as.integer(k), nrow(pixel_matrix)))
  matrix_data <- as.matrix(pixel_matrix[fcols])
  pca <- NULL
  if (!is.null(pca_components) && is.finite(pca_components) && pca_components > 0) {
    n_components <- min(as.integer(pca_components), ncol(matrix_data), nrow(matrix_data) - 1)
    if (n_components >= 1) {
      pca <- stats::prcomp(matrix_data, center = TRUE, scale. = FALSE, rank. = n_components)
      matrix_data <- pca$x[, seq_len(n_components), drop = FALSE]
    }
  }
  if (!is.null(seed)) set.seed(seed)
  fit <- stats::kmeans(
    matrix_data,
    centers = k,
    nstart = nstart,
    iter.max = iter.max,
    algorithm = algorithm
  )
  out <- pixel_matrix
  out$cluster <- fit$cluster
  list(matrix = out, fit = fit, pca = pca, pca_components = if (is.null(pca)) 0 else ncol(matrix_data))
}

# Compute the H0 (mean-neighborhood) feature block used by the
# BANKSY-inspired clustering helper below. This is not the full BANKSY
# representation: no azimuthal Gabor feature (AGF) block is computed here.
compute_neighborhood_average <- function(feature_matrix,
                                         neighbors,
                                         isolate_action = c("self", "error")) {
  isolate_action <- match.arg(isolate_action)
  if (!inherits(neighbors, "spatial_neighbors")) {
    stop("neighbors must come from build_spatial_neighbors().", call. = FALSE)
  }
  feature_matrix <- as.matrix(feature_matrix)
  storage.mode(feature_matrix) <- "double"
  n_nodes <- nrow(feature_matrix)
  if (n_nodes != length(neighbors$x)) {
    stop("feature_matrix must have one row per spatial-neighbor node.", call. = FALSE)
  }
  if (!ncol(feature_matrix)) {
    stop("feature_matrix must contain at least one feature.", call. = FALSE)
  }
  if (any(!is.finite(feature_matrix))) {
    stop("feature_matrix must contain only finite values.", call. = FALSE)
  }

  edges <- neighbors$edges
  if (!nrow(edges)) {
    stop("Spatial neighbor graph has no directed edges.", call. = FALSE)
  }
  if (!isTRUE(neighbors$symmetric)) {
    stop("Neighborhood averaging requires a symmetric neighbor graph.", call. = FALSE)
  }
  if (any(!is.finite(edges$weight)) || any(edges$weight < 0)) {
    stop("Spatial neighbor weights must be finite and non-negative.", call. = FALSE)
  }

  neighbor_sum <- matrix(0, nrow = n_nodes, ncol = ncol(feature_matrix))
  weight_sum <- numeric(n_nodes)
  weighted_values <- feature_matrix[edges$to, , drop = FALSE] * edges$weight
  grouped_values <- rowsum(weighted_values, group = edges$from, reorder = FALSE)
  grouped_weights <- rowsum(matrix(edges$weight, ncol = 1L), group = edges$from, reorder = FALSE)
  grouped_nodes <- as.integer(rownames(grouped_values))
  neighbor_sum[grouped_nodes, ] <- grouped_values
  weight_sum[grouped_nodes] <- grouped_weights[, 1]
  isolated <- weight_sum <= 0
  if (any(isolated) && isolate_action == "error") {
    stop(sum(isolated), " spatial node(s) have no positive-weight neighbors.", call. = FALSE)
  }
  connected <- !isolated
  neighbor_sum[connected, ] <- neighbor_sum[connected, , drop = FALSE] / weight_sum[connected]
  if (any(isolated)) {
    neighbor_sum[isolated, ] <- feature_matrix[isolated, , drop = FALSE]
  }
  dimnames(neighbor_sum) <- dimnames(feature_matrix)
  attr(neighbor_sum, "isolated") <- isolated
  neighbor_sum
}

scale_feature_block <- function(values) {
  scaled <- scale(values, center = TRUE, scale = TRUE)
  scaled[!is.finite(scaled)] <- 0
  unname(scaled)
}

# BANKSY-inspired H0-only spatial feature augmentation. It concatenates
# z-scaled own and mean-neighborhood features with sqrt(1-lambda) and
# sqrt(lambda) weights. Unlike full BANKSY, this helper does not compute AGF,
# does not construct the paper's default kNN Gaussian kernel, and uses k-means
# rather than Leiden community detection.
cluster_pixels_spatial <- function(pixel_matrix,
                                   k = 3,
                                   lambda = 0.5,
                                   neighbor_method = c("queen", "rook", "distance"),
                                   distance_threshold = NULL,
                                   neighbors = NULL,
                                   isolate_action = c("self", "error"),
                                   pca_components = NULL,
                                   nstart = 25,
                                   iter.max = 100,
                                   algorithm = "Lloyd",
                                   seed = 1) {
  neighbor_method <- match.arg(neighbor_method)
  isolate_action <- match.arg(isolate_action)
  if (length(lambda) != 1L || !is.finite(lambda) || lambda < 0 || lambda > 1) {
    stop("lambda must be one finite value between 0 and 1.", call. = FALSE)
  }
  required_columns(pixel_matrix, c("x", "y"), "Pixel matrix")
  if (nrow(pixel_matrix) < 2L) stop("At least two pixels are required for clustering.", call. = FALSE)
  fcols <- feature_columns(pixel_matrix)
  if (!length(fcols)) stop("No mz feature columns available for clustering.", call. = FALSE)
  own_features <- as.matrix(pixel_matrix[fcols])
  storage.mode(own_features) <- "double"
  if (any(!is.finite(own_features))) {
    stop("Feature values must be finite before spatial clustering.", call. = FALSE)
  }

  if (is.null(neighbors)) {
    neighbors <- build_spatial_neighbors(
      pixel_matrix$x, pixel_matrix$y,
      method = neighbor_method,
      distance_threshold = distance_threshold,
      symmetric = TRUE
    )
  } else {
    if (!inherits(neighbors, "spatial_neighbors")) {
      stop("neighbors must come from build_spatial_neighbors().", call. = FALSE)
    }
    if (!isTRUE(all.equal(as.numeric(pixel_matrix$x), neighbors$x)) ||
        !isTRUE(all.equal(as.numeric(pixel_matrix$y), neighbors$y))) {
      stop("neighbors coordinates and pixel_matrix row order do not match.", call. = FALSE)
    }
  }

  neighbor_features <- compute_neighborhood_average(
    own_features, neighbors, isolate_action = isolate_action
  )
  own_scaled <- scale_feature_block(own_features)
  neighbor_scaled <- scale_feature_block(neighbor_features)
  augmented <- cbind(
    sqrt(1 - lambda) * own_scaled,
    sqrt(lambda) * neighbor_scaled
  )
  colnames(augmented) <- c(paste0("own__", fcols), paste0("neighbor__", fcols))

  matrix_data <- augmented
  pca <- NULL
  if (!is.null(pca_components) && length(pca_components) == 1L &&
      is.finite(pca_components) && pca_components > 0) {
    n_components <- min(as.integer(pca_components), ncol(matrix_data), nrow(matrix_data) - 1L)
    if (n_components >= 1L) {
      pca <- stats::prcomp(matrix_data, center = TRUE, scale. = FALSE, rank. = n_components)
      matrix_data <- pca$x[, seq_len(n_components), drop = FALSE]
    }
  }

  k <- as.integer(k)
  if (length(k) != 1L || is.na(k) || k < 2L || k > nrow(pixel_matrix)) {
    stop("k must be an integer from 2 through the number of pixels.", call. = FALSE)
  }
  if (!is.null(seed)) set.seed(seed)
  fit <- stats::kmeans(
    matrix_data, centers = k, nstart = nstart, iter.max = iter.max,
    algorithm = algorithm
  )
  out <- pixel_matrix
  out$cluster <- fit$cluster
  list(
    matrix = out,
    fit = fit,
    pca = pca,
    neighbors = neighbors,
    neighborhood_features = neighbor_features,
    augmented_features = augmented,
    lambda = lambda,
    neighbor_method = neighbors$method,
    isolate_action = isolate_action,
    n_isolated = sum(attr(neighbor_features, "isolated")),
    pca_components = if (is.null(pca)) 0L else ncol(matrix_data),
    method = "BANKSY-inspired H0-only spatial augmentation followed by k-means",
    full_banksy = FALSE
  )
}

cluster_diagnostics_spatial <- function(pixel_matrix,
                                        k = 3,
                                        lambda_grid = seq(0, 1, by = 0.25),
                                        neighbor_method = c("queen", "rook", "distance"),
                                        distance_threshold = NULL,
                                        isolate_action = c("self", "error"),
                                        pca_components = NULL,
                                        nstart = 25,
                                        iter.max = 100,
                                        algorithm = "Lloyd",
                                        seed = 1) {
  neighbor_method <- match.arg(neighbor_method)
  isolate_action <- match.arg(isolate_action)
  if (!length(lambda_grid) || any(!is.finite(lambda_grid)) || any(lambda_grid < 0 | lambda_grid > 1)) {
    stop("lambda_grid must contain finite values between 0 and 1.", call. = FALSE)
  }
  neighbors <- build_spatial_neighbors(
    pixel_matrix$x, pixel_matrix$y, method = neighbor_method,
    distance_threshold = distance_threshold, symmetric = TRUE
  )
  edges <- neighbors$undirected_edges
  if (!nrow(edges)) stop("Spatial neighbor graph has no edges.", call. = FALSE)

  rows <- lapply(lambda_grid, function(lambda) {
    result <- cluster_pixels_spatial(
      pixel_matrix, k = k, lambda = lambda, neighbors = neighbors,
      isolate_action = isolate_action, pca_components = pca_components,
      nstart = nstart, iter.max = iter.max, algorithm = algorithm, seed = seed
    )
    labels <- result$matrix$cluster
    same_cluster <- labels[edges$from] == labels[edges$to]
    data.frame(
      lambda = lambda,
      k = k,
      tot_withinss = result$fit$tot.withinss,
      adjacent_pair_agreement = mean(same_cluster),
      boundary_edge_fraction = mean(!same_cluster),
      n_edges_undirected = nrow(edges),
      n_isolated = result$n_isolated,
      stringsAsFactors = FALSE
    )
  })
  output <- do.call(rbind, rows)
  attr(output, "interpretation") <- paste(
    "Adjacent-pair agreement is descriptive and generally favors smoother",
    "solutions; it is not an accuracy metric and should not select lambda alone."
  )
  output
}

# Gaussian/Potts-inspired spatially regularized k-means optimized by
# sequential iterated conditional modes (ICM). This is a deterministic local
# optimizer for the stated objective, not full Bayesian HMRF inference.
cluster_pixels_hmrf <- function(pixel_matrix,
                                k = 3,
                                beta = 1,
                                neighbor_method = c("queen", "rook", "distance"),
                                distance_threshold = NULL,
                                neighbors = NULL,
                                scale_features = TRUE,
                                data_term_scale = c("per_feature", "raw"),
                                pca_components = NULL,
                                max_iter = 20,
                                update_order = c("random", "fixed"),
                                init_nstart = 25,
                                init_iter_max = 100,
                                init_algorithm = "Lloyd",
                                seed = 1,
                                energy_tolerance = 1e-10) {
  neighbor_method <- match.arg(neighbor_method)
  data_term_scale <- match.arg(data_term_scale)
  update_order <- match.arg(update_order)
  if (length(beta) != 1L || !is.finite(beta) || beta < 0) {
    stop("beta must be one finite non-negative value.", call. = FALSE)
  }
  if (length(max_iter) != 1L || !is.finite(max_iter) || max_iter < 1) {
    stop("max_iter must be one positive integer.", call. = FALSE)
  }
  max_iter <- as.integer(max_iter)
  if (length(energy_tolerance) != 1L || !is.finite(energy_tolerance) || energy_tolerance < 0) {
    stop("energy_tolerance must be one finite non-negative value.", call. = FALSE)
  }
  required_columns(pixel_matrix, c("x", "y"), "Pixel matrix")
  if (nrow(pixel_matrix) < 2L) stop("At least two pixels are required for clustering.", call. = FALSE)
  fcols <- feature_columns(pixel_matrix)
  if (!length(fcols)) stop("No mz feature columns available for clustering.", call. = FALSE)
  matrix_data <- as.matrix(pixel_matrix[fcols])
  storage.mode(matrix_data) <- "double"
  if (any(!is.finite(matrix_data))) {
    stop("Feature values must be finite before HMRF/Potts clustering.", call. = FALSE)
  }
  if (isTRUE(scale_features)) matrix_data <- scale_feature_block(matrix_data)

  pca <- NULL
  if (!is.null(pca_components) && length(pca_components) == 1L &&
      is.finite(pca_components) && pca_components > 0) {
    n_components <- min(as.integer(pca_components), ncol(matrix_data), nrow(matrix_data) - 1L)
    if (n_components >= 1L) {
      pca <- stats::prcomp(matrix_data, center = TRUE, scale. = FALSE, rank. = n_components)
      matrix_data <- pca$x[, seq_len(n_components), drop = FALSE]
    }
  }
  n <- nrow(matrix_data)
  k <- as.integer(k)
  if (length(k) != 1L || is.na(k) || k < 2L || k > n) {
    stop("k must be an integer from 2 through the number of pixels.", call. = FALSE)
  }

  if (is.null(neighbors)) {
    neighbors <- build_spatial_neighbors(
      pixel_matrix$x, pixel_matrix$y, method = neighbor_method,
      distance_threshold = distance_threshold, symmetric = TRUE
    )
  } else {
    if (!inherits(neighbors, "spatial_neighbors")) {
      stop("neighbors must come from build_spatial_neighbors().", call. = FALSE)
    }
    if (!isTRUE(neighbors$symmetric)) stop("HMRF/Potts clustering requires a symmetric graph.", call. = FALSE)
    if (!isTRUE(all.equal(as.numeric(pixel_matrix$x), neighbors$x)) ||
        !isTRUE(all.equal(as.numeric(pixel_matrix$y), neighbors$y))) {
      stop("neighbors coordinates and pixel_matrix row order do not match.", call. = FALSE)
    }
  }
  edges <- neighbors$edges
  undirected_edges <- neighbors$undirected_edges
  if (!nrow(undirected_edges)) stop("Spatial neighbor graph has no edges.", call. = FALSE)
  if (any(!is.finite(edges$weight)) || any(edges$weight < 0)) {
    stop("Spatial neighbor weights must be finite and non-negative.", call. = FALSE)
  }
  adjacency_to <- split(edges$to, factor(edges$from, levels = seq_len(n)))
  adjacency_weight <- split(edges$weight, factor(edges$from, levels = seq_len(n)))
  data_divisor <- if (data_term_scale == "per_feature") ncol(matrix_data) else 1

  if (!is.null(seed)) set.seed(seed)
  init_fit <- stats::kmeans(
    matrix_data, centers = k, nstart = init_nstart,
    iter.max = init_iter_max, algorithm = init_algorithm
  )
  labels <- init_fit$cluster
  centers <- init_fit$centers

  objective_parts <- function(current_labels, current_centers) {
    residual <- matrix_data - current_centers[current_labels, , drop = FALSE]
    data_energy <- sum(residual^2) / data_divisor
    spatial_energy <- beta * sum(
      undirected_edges$weight *
        (current_labels[undirected_edges$from] != current_labels[undirected_edges$to])
    )
    c(data_energy = data_energy, spatial_energy = spatial_energy,
      total_energy = data_energy + spatial_energy)
  }
  initial_energy <- objective_parts(labels, centers)
  iteration_log <- data.frame(
    iteration = 0L, n_changed = NA_integer_,
    data_energy = initial_energy[["data_energy"]],
    spatial_energy = initial_energy[["spatial_energy"]],
    total_energy = initial_energy[["total_energy"]],
    stringsAsFactors = FALSE
  )
  converged <- FALSE

  for (iteration in seq_len(max_iter)) {
    previous_labels <- labels
    node_order <- if (update_order == "random") sample.int(n) else seq_len(n)
    class_sizes <- tabulate(labels, nbins = k)
    for (node in node_order) {
      current_class <- labels[node]
      data_energy <- rowSums(
        sweep(centers, 2, matrix_data[node, ], "-")^2
      ) / data_divisor
      neighbor_nodes <- adjacency_to[[node]]
      neighbor_weights <- adjacency_weight[[node]]
      same_weight <- numeric(k)
      if (length(neighbor_nodes)) {
        same_weight <- vapply(seq_len(k), function(class_id) {
          sum(neighbor_weights[labels[neighbor_nodes] == class_id])
        }, numeric(1))
      }
      conditional_energy <- data_energy - beta * same_weight
      minimum <- min(conditional_energy)
      candidates <- which(abs(conditional_energy - minimum) <= energy_tolerance)
      proposed <- if (current_class %in% candidates) current_class else candidates[1L]
      # Preserve all k states; an empty Gaussian component has no defined mean.
      if (proposed != current_class && class_sizes[current_class] > 1L) {
        labels[node] <- proposed
        class_sizes[current_class] <- class_sizes[current_class] - 1L
        class_sizes[proposed] <- class_sizes[proposed] + 1L
      }
    }

    centers <- t(vapply(seq_len(k), function(class_id) {
      colMeans(matrix_data[labels == class_id, , drop = FALSE])
    }, numeric(ncol(matrix_data))))
    energy <- objective_parts(labels, centers)
    n_changed <- sum(labels != previous_labels)
    iteration_log <- rbind(iteration_log, data.frame(
      iteration = iteration, n_changed = n_changed,
      data_energy = energy[["data_energy"]],
      spatial_energy = energy[["spatial_energy"]],
      total_energy = energy[["total_energy"]],
      stringsAsFactors = FALSE
    ))
    previous_energy <- iteration_log$total_energy[nrow(iteration_log) - 1L]
    audit_tolerance <- energy_tolerance * max(1, abs(previous_energy))
    if (energy[["total_energy"]] > previous_energy + audit_tolerance) {
      stop("Internal error: sequential ICM increased the Potts objective.", call. = FALSE)
    }
    if (n_changed == 0L) {
      converged <- TRUE
      break
    }
  }

  withinss <- vapply(seq_len(k), function(class_id) {
    members <- matrix_data[labels == class_id, , drop = FALSE]
    sum(sweep(members, 2, centers[class_id, ], "-")^2)
  }, numeric(1))
  totss <- sum(scale(matrix_data, center = TRUE, scale = FALSE)^2)
  fit <- list(
    cluster = labels, centers = centers, totss = totss,
    withinss = withinss, tot.withinss = sum(withinss),
    betweenss = totss - sum(withinss), size = tabulate(labels, nbins = k),
    iter = max(iteration_log$iteration), ifault = if (converged) 0L else 2L
  )
  out <- pixel_matrix
  out$cluster <- labels
  list(
    matrix = out, labels = labels, centers = centers, fit = fit,
    beta = beta, neighbor_method = neighbors$method, neighbors = neighbors,
    iteration_log = iteration_log, converged = converged,
    update_order = update_order, scale_features = isTRUE(scale_features),
    data_term_scale = data_term_scale, pca = pca,
    pca_components = if (is.null(pca)) 0L else ncol(matrix_data),
    method = "Gaussian/Potts-inspired spatially regularized k-means with sequential ICM",
    full_bayesian_hmrf = FALSE
  )
}

cluster_diagnostics_hmrf <- function(pixel_matrix,
                                     k = 3,
                                     beta_grid = c(0, 0.5, 1, 2, 4),
                                     neighbor_method = c("queen", "rook", "distance"),
                                     distance_threshold = NULL,
                                     scale_features = TRUE,
                                     data_term_scale = c("per_feature", "raw"),
                                     pca_components = NULL,
                                     max_iter = 20,
                                     update_order = c("random", "fixed"),
                                     seed = 1) {
  neighbor_method <- match.arg(neighbor_method)
  data_term_scale <- match.arg(data_term_scale)
  update_order <- match.arg(update_order)
  if (!length(beta_grid) || any(!is.finite(beta_grid)) || any(beta_grid < 0)) {
    stop("beta_grid must contain finite non-negative values.", call. = FALSE)
  }
  required_columns(pixel_matrix, c("x", "y"), "Pixel matrix")
  neighbors <- build_spatial_neighbors(
    pixel_matrix$x, pixel_matrix$y, method = neighbor_method,
    distance_threshold = distance_threshold, symmetric = TRUE
  )
  edges <- neighbors$undirected_edges
  if (!nrow(edges)) stop("Spatial neighbor graph has no edges.", call. = FALSE)
  rows <- lapply(beta_grid, function(beta) {
    result <- cluster_pixels_hmrf(
      pixel_matrix, k = k, beta = beta, neighbors = neighbors,
      scale_features = scale_features, data_term_scale = data_term_scale,
      pca_components = pca_components, max_iter = max_iter,
      update_order = update_order, seed = seed
    )
    labels <- result$labels
    same_cluster <- labels[edges$from] == labels[edges$to]
    final <- result$iteration_log[nrow(result$iteration_log), , drop = FALSE]
    data.frame(
      beta = beta, k = k,
      data_energy = final$data_energy,
      spatial_energy = final$spatial_energy,
      total_energy = final$total_energy,
      tot_withinss = result$fit$tot.withinss,
      n_iterations = final$iteration,
      converged = result$converged,
      adjacent_pair_agreement = mean(same_cluster),
      boundary_edge_fraction = mean(!same_cluster),
      n_edges_undirected = nrow(edges),
      stringsAsFactors = FALSE
    )
  })
  output <- do.call(rbind, rows)
  attr(output, "interpretation") <- paste(
    "Beta is scale-dependent. Adjacent-pair agreement measures smoothness,",
    "not accuracy, and must not be the sole beta-selection criterion."
  )
  output
}

cluster_diagnostics <- function(pixel_matrix,
                                max_k = 10,
                                pca_components = NULL,
                                nstart = 25,
                                iter.max = 100,
                                algorithm = "Lloyd",
                                seed = 1,
                                max_subsample = 2000) {
  fcols <- feature_columns(pixel_matrix)
  max_k <- max(2, min(max_k, nrow(pixel_matrix) - 1))
  matrix_data_full <- as.matrix(pixel_matrix[fcols])
  n_pixels_total <- nrow(matrix_data_full)
  pca_components_used <- 0L
  if (!is.null(pca_components) && is.finite(pca_components) && pca_components > 0) {
    n_components <- min(as.integer(pca_components), ncol(matrix_data_full), nrow(matrix_data_full) - 1)
    if (n_components >= 1) {
      pca <- stats::prcomp(matrix_data_full, center = TRUE, scale. = FALSE, rank. = n_components)
      matrix_data_full <- pca$x[, seq_len(n_components), drop = FALSE]
      pca_components_used <- n_components
    }
  }
  sample_index <- seq_len(nrow(matrix_data_full))
  subsampled <- FALSE
  if (!is.null(max_subsample) && is.finite(max_subsample) && nrow(matrix_data_full) > max_subsample) {
    if (!is.null(seed)) set.seed(seed)
    sample_index <- sample(seq_len(nrow(matrix_data_full)), as.integer(max_subsample))
    subsampled <- TRUE
  }
  matrix_data_used <- matrix_data_full[sample_index, , drop = FALSE]

  rows <- lapply(2:max_k, function(k) {
    if (!is.null(seed)) set.seed(seed)
    matrix_data <- matrix_data_used
    fit <- stats::kmeans(matrix_data, centers = k, nstart = nstart, iter.max = iter.max, algorithm = algorithm)
    sil_mean <- NA_real_
    silhouette_limit <- if (is.null(max_subsample) || !is.finite(max_subsample)) Inf else max_subsample
    if (nrow(matrix_data) <= silhouette_limit) {
      if (!requireNamespace("cluster", quietly = TRUE)) {
        stop("Package 'cluster' is required for silhouette computation.", call. = FALSE)
      }
      dist_obj <- stats::dist(matrix_data)
      sil <- cluster::silhouette(fit$cluster, dist_obj)
      sil_mean <- mean(sil[, "sil_width"])
    }
    data.frame(
      k = k,
      tot_withinss = fit$tot.withinss,
      silhouette = sil_mean,
      n_pixels_total = n_pixels_total,
      n_pixels_used = nrow(matrix_data),
      subsampled = subsampled,
      pca_components = pca_components_used
    )
  })
  do.call(rbind, rows)
}

compute_neighborhood_composition <- function(cell_type_labels,
                                             x,
                                             y,
                                             k = 10,
                                             cell_types = NULL,
                                             pixel_id = seq_along(cell_type_labels),
                                             region_id = rep("region_1", length(cell_type_labels)),
                                             include_self = TRUE) {
  n <- length(cell_type_labels)
  if (n == 0L || length(x) != n || length(y) != n ||
      length(pixel_id) != n || length(region_id) != n) {
    stop(
      "cell_type_labels, x, y, pixel_id, and region_id must have the same non-zero length.",
      call. = FALSE
    )
  }
  if (anyNA(cell_type_labels) || any(!nzchar(as.character(cell_type_labels)))) {
    stop("cell_type_labels must not contain missing or empty values.", call. = FALSE)
  }
  if (!is.numeric(x) || !is.numeric(y) || any(!is.finite(x)) || any(!is.finite(y))) {
    stop("x and y must contain only finite numeric coordinates.", call. = FALSE)
  }
  if (anyNA(pixel_id) || anyDuplicated(pixel_id)) {
    stop("pixel_id must be non-missing and unique.", call. = FALSE)
  }
  if (anyNA(region_id) || any(!nzchar(as.character(region_id)))) {
    stop("region_id must not contain missing or empty values.", call. = FALSE)
  }
  if (length(k) != 1L || !is.finite(k) || k < 1 || k != as.integer(k)) {
    stop("k must be one positive integer.", call. = FALSE)
  }
  if (length(include_self) != 1L || is.na(include_self)) {
    stop("include_self must be TRUE or FALSE.", call. = FALSE)
  }
  k <- as.integer(k)
  labels <- as.character(cell_type_labels)
  if (is.null(cell_types)) {
    cell_types <- sort(unique(labels))
  } else {
    cell_types <- as.character(cell_types)
    if (anyNA(cell_types) || any(!nzchar(cell_types)) || anyDuplicated(cell_types)) {
      stop("cell_types must contain unique, non-missing, non-empty labels.", call. = FALSE)
    }
    missing_types <- setdiff(unique(labels), cell_types)
    if (length(missing_types) > 0L) {
      stop("cell_types is missing observed label(s): ", paste(missing_types, collapse = ", "), call. = FALSE)
    }
  }

  region <- as.character(region_id)
  region_sizes <- table(region)
  minimum_size <- k + if (isTRUE(include_self)) 0L else 1L
  too_small <- names(region_sizes)[region_sizes < minimum_size]
  if (length(too_small) > 0L) {
    stop(
      "Each region must contain at least ", minimum_size,
      " positions for this k/include_self setting. Too small: ",
      paste(too_small, collapse = ", "),
      call. = FALSE
    )
  }

  safe_type_names <- make.unique(make.names(cell_types))
  composition_columns <- paste0("composition__", safe_type_names)
  count_matrix <- matrix(0L, nrow = n, ncol = length(cell_types))
  colnames(count_matrix) <- composition_columns
  neighbor_indices <- matrix(NA_integer_, nrow = n, ncol = k)
  x_numeric <- as.numeric(x)
  y_numeric <- as.numeric(y)

  for (region_name in unique(region)) {
    region_indices <- which(region == region_name)
    for (i in region_indices) {
      candidates <- setdiff(region_indices, i)
      distance_squared <- (x_numeric[candidates] - x_numeric[i])^2 +
        (y_numeric[candidates] - y_numeric[i])^2
      ordered_candidates <- candidates[order(distance_squared, candidates)]
      selected <- if (isTRUE(include_self)) {
        c(i, ordered_candidates[seq_len(k - 1L)])
      } else {
        ordered_candidates[seq_len(k)]
      }
      neighbor_indices[i, ] <- selected
      count_matrix[i, ] <- tabulate(match(labels[selected], cell_types), nbins = length(cell_types))
    }
  }

  composition_matrix <- count_matrix / k
  matrix_output <- data.frame(
    pixel_id = pixel_id,
    region_id = region,
    focal_cell_type = labels,
    window_size = rep(k, n),
    composition_matrix,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  mapping <- data.frame(
    cell_type = cell_types,
    composition_column = composition_columns,
    stringsAsFactors = FALSE
  )
  structure(
    list(
      matrix = matrix_output,
      counts = count_matrix,
      neighbor_indices = neighbor_indices,
      cell_type_mapping = mapping,
      settings = list(k = k, include_self = isTRUE(include_self), region_column = "region_id")
    ),
    class = "neighborhood_composition"
  )
}

define_niches <- function(composition,
                          k_niches,
                          nstart = 25,
                          iter.max = 100,
                          algorithm = "Lloyd",
                          seed = 1) {
  if (inherits(composition, "neighborhood_composition")) {
    composition_matrix <- composition$matrix
    composition_columns <- composition$cell_type_mapping$composition_column
  } else if (is.data.frame(composition)) {
    composition_matrix <- composition
    composition_columns <- grep("^composition__", names(composition_matrix), value = TRUE)
  } else {
    stop("composition must come from compute_neighborhood_composition() or be a data frame.", call. = FALSE)
  }
  if (length(composition_columns) == 0L ||
      any(!composition_columns %in% names(composition_matrix))) {
    stop("No valid composition__ columns were found.", call. = FALSE)
  }
  values <- as.matrix(composition_matrix[, composition_columns, drop = FALSE])
  storage.mode(values) <- "double"
  if (nrow(values) < 2L || any(!is.finite(values)) || any(values < 0)) {
    stop("Composition values must be a finite, non-negative matrix with at least two rows.", call. = FALSE)
  }
  row_totals <- rowSums(values)
  if (any(abs(row_totals - 1) > sqrt(.Machine$double.eps))) {
    stop("Each neighborhood composition row must sum to one.", call. = FALSE)
  }
  if (length(k_niches) != 1L || !is.finite(k_niches) ||
      k_niches < 2 || k_niches != as.integer(k_niches)) {
    stop("k_niches must be one integer of at least two.", call. = FALSE)
  }
  k_niches <- as.integer(k_niches)
  unique_profiles <- nrow(unique(as.data.frame(values, check.names = FALSE)))
  if (k_niches > unique_profiles) {
    stop(
      "k_niches cannot exceed the number of distinct composition profiles (",
      unique_profiles, ").",
      call. = FALSE
    )
  }
  if (!is.null(seed)) set.seed(seed)
  fit <- stats::kmeans(
    values,
    centers = k_niches,
    nstart = nstart,
    iter.max = iter.max,
    algorithm = algorithm
  )
  output <- composition_matrix
  output$niche_id <- as.integer(fit$cluster)
  output$niche_label <- paste0("niche_", fit$cluster)
  list(
    matrix = output,
    fit = fit,
    centers = fit$centers,
    composition_columns = composition_columns,
    settings = list(
      k_niches = k_niches,
      nstart = nstart,
      iter.max = iter.max,
      algorithm = algorithm,
      seed = seed
    ),
    interpretation = paste(
      "Niche numbers are arbitrary k-means labels, not ordered biological identities;",
      "annotate them from their cell-type composition centers."
    )
  )
}

define_rois <- function(pixel_matrix,
                        mode = c("cluster", "coordinate"),
                        cluster_mapping = NULL,
                        roi_table = NULL,
                        roi_column = "roi_id") {
  mode <- match.arg(mode)
  required_columns(pixel_matrix, c("pixel_id", "x", "y"), "Pixel matrix")
  out <- pixel_matrix
  out[[roi_column]] <- NA_character_

  if (mode == "cluster") {
    required_columns(pixel_matrix, "cluster", "Clustered matrix")
    clusters <- sort(unique(pixel_matrix$cluster[!is.na(pixel_matrix$cluster)]))
    if (is.null(cluster_mapping)) {
      cluster_mapping <- data.frame(
        cluster = clusters,
        roi_id = paste0("roi_cluster_", clusters),
        stringsAsFactors = FALSE
      )
    }
    required_columns(cluster_mapping, c("cluster", "roi_id"), "Cluster ROI mapping")
    for (i in seq_len(nrow(cluster_mapping))) {
      cluster_value <- cluster_mapping$cluster[i]
      roi_id <- as.character(cluster_mapping$roi_id[i])
      out[[roi_column]][pixel_matrix$cluster == cluster_value] <- roi_id
    }
    return(out)
  }

  if (is.null(roi_table)) {
    stop("roi_table is required when mode = 'coordinate'.", call. = FALSE)
  }
  required_columns(roi_table, c("roi_id", "x_min", "x_max", "y_min", "y_max"), "Coordinate ROI table")
  roi_table <- roi_table[!is.na(roi_table$roi_id) & roi_table$roi_id != "", , drop = FALSE]
  if (nrow(roi_table) == 0) {
    stop("roi_table does not contain any non-empty roi_id values.", call. = FALSE)
  }

  for (i in seq_len(nrow(roi_table))) {
    roi <- roi_table[i, , drop = FALSE]
    x_min <- min(roi$x_min, roi$x_max, na.rm = TRUE)
    x_max <- max(roi$x_min, roi$x_max, na.rm = TRUE)
    y_min <- min(roi$y_min, roi$y_max, na.rm = TRUE)
    y_max <- max(roi$y_min, roi$y_max, na.rm = TRUE)
    if (!all(is.finite(c(x_min, x_max, y_min, y_max)))) {
      warning("Skipping ROI with non-finite coordinates: ", roi$roi_id[1], call. = FALSE)
      next
    }
    in_roi <- pixel_matrix$x >= x_min &
      pixel_matrix$x <= x_max &
      pixel_matrix$y >= y_min &
      pixel_matrix$y <= y_max
    out[[roi_column]][in_roi] <- as.character(roi$roi_id[1])
  }

  out
}

sample_subregions <- function(pixel_matrix,
                              grid_size = 5,
                              min_pixels = 30,
                              roi_column = "roi_id",
                              grid_scope = c("global", "roi")) {
  grid_scope <- match.arg(grid_scope)
  required_columns(pixel_matrix, c("pixel_id", "x", "y", roi_column), "ROI pixel matrix")
  fcols <- feature_columns(pixel_matrix)
  if (length(fcols) == 0) {
    stop("No mz feature columns available for ROI subregion sampling.", call. = FALSE)
  }
  grid_size <- max(2L, as.integer(grid_size))
  min_pixels <- max(1L, as.integer(min_pixels))
  work <- pixel_matrix[!is.na(pixel_matrix[[roi_column]]) & pixel_matrix[[roi_column]] != "", , drop = FALSE]
  if (nrow(work) == 0) {
    stop("No pixels have ROI labels. Run define_rois() before sample_subregions().", call. = FALSE)
  }
  if (grid_scope == "global") {
    work$grid_x <- as.integer(cut(work$x, breaks = grid_size, labels = FALSE, include.lowest = TRUE))
    work$grid_y <- as.integer(cut(work$y, breaks = grid_size, labels = FALSE, include.lowest = TRUE))
  } else {
    work$grid_x <- NA_integer_
    work$grid_y <- NA_integer_
    for (roi_id in unique(work[[roi_column]])) {
      in_roi <- work[[roi_column]] == roi_id
      work$grid_x[in_roi] <- as.integer(cut(work$x[in_roi], breaks = grid_size, labels = FALSE, include.lowest = TRUE))
      work$grid_y[in_roi] <- as.integer(cut(work$y[in_roi], breaks = grid_size, labels = FALSE, include.lowest = TRUE))
    }
  }
  work$grid_cell <- paste(work$grid_x, work$grid_y, sep = "_")

  split_groups <- split(work, interaction(work[[roi_column]], work$grid_cell, drop = TRUE), drop = TRUE)
  sample_rows <- list()
  mapping_rows <- list()

  for (group_name in names(split_groups)) {
    group_data <- split_groups[[group_name]]
    n_pixels <- nrow(group_data)
    if (n_pixels < min_pixels) next

    roi_id <- as.character(group_data[[roi_column]][1])
    grid_cell <- group_data$grid_cell[1]
    sample_id <- paste0(utils::URLencode(roi_id, reserved = TRUE), "_", grid_cell)
    means <- colMeans(group_data[fcols], na.rm = TRUE)

    sample_rows[[sample_id]] <- data.frame(
      sample_id = sample_id,
      roi_id = roi_id,
      grid_cell = grid_cell,
      grid_scope = grid_scope,
      grid_x = group_data$grid_x[1],
      grid_y = group_data$grid_y[1],
      x_min = min(group_data$x, na.rm = TRUE),
      x_max = max(group_data$x, na.rm = TRUE),
      y_min = min(group_data$y, na.rm = TRUE),
      y_max = max(group_data$y, na.rm = TRUE),
      x_span = diff(range(group_data$x, na.rm = TRUE)),
      y_span = diff(range(group_data$y, na.rm = TRUE)),
      n_pixels = n_pixels,
      t(means),
      check.names = FALSE
    )
    mapping_rows[[sample_id]] <- data.frame(
      sample_id = sample_id,
      roi_id = roi_id,
      grid_cell = grid_cell,
      grid_scope = grid_scope,
      grid_x = group_data$grid_x[1],
      grid_y = group_data$grid_y[1],
      x_min = min(group_data$x, na.rm = TRUE),
      x_max = max(group_data$x, na.rm = TRUE),
      y_min = min(group_data$y, na.rm = TRUE),
      y_max = max(group_data$y, na.rm = TRUE),
      x_span = diff(range(group_data$x, na.rm = TRUE)),
      y_span = diff(range(group_data$y, na.rm = TRUE)),
      n_pixels = n_pixels,
      pixel_ids = paste(group_data$pixel_id, collapse = ","),
      stringsAsFactors = FALSE
    )
  }

  sample_matrix <- if (length(sample_rows)) do.call(rbind, sample_rows) else data.frame()
  sample_mapping <- if (length(mapping_rows)) do.call(rbind, mapping_rows) else data.frame()
  rownames(sample_matrix) <- NULL
  rownames(sample_mapping) <- NULL
  list(sample_matrix = sample_matrix, sample_mapping = sample_mapping, annotated_pixels = work)
}

sample_coordinate_rois <- function(pixel_matrix,
                                   roi_table,
                                   label_column = "roi_label",
                                   min_pixels = 30) {
  required_columns(pixel_matrix, c("pixel_id", "x", "y"), "Pixel matrix")
  required_columns(roi_table, c(label_column, "x_min", "x_max", "y_min", "y_max"), "ROI table")
  warning(
    "sample_coordinate_rois() is deprecated. Use define_rois(mode = 'coordinate') ",
    "followed by sample_subregions(), or sample by roi_id directly.",
    call. = FALSE
  )
  fcols <- feature_columns(pixel_matrix)
  if (length(fcols) == 0) {
    stop("No mz feature columns available for coordinate ROI sampling.", call. = FALSE)
  }

  min_pixels <- max(1L, as.integer(min_pixels))
  names(roi_table)[names(roi_table) == label_column] <- "roi_id"
  labeled <- define_rois(pixel_matrix, mode = "coordinate", roi_table = roi_table)
  roi_table <- roi_table[!is.na(roi_table$roi_id) & roi_table$roi_id != "", , drop = FALSE]
  if (nrow(roi_table) == 0) {
    stop("ROI table does not contain any non-empty ROI labels.", call. = FALSE)
  }

  annotated_rows <- list()
  sample_rows <- list()
  mapping_rows <- list()

  for (i in seq_len(nrow(roi_table))) {
    roi <- roi_table[i, , drop = FALSE]
    roi_label <- as.character(roi$roi_id[1])
    x_min <- min(roi$x_min, roi$x_max, na.rm = TRUE)
    x_max <- max(roi$x_min, roi$x_max, na.rm = TRUE)
    y_min <- min(roi$y_min, roi$y_max, na.rm = TRUE)
    y_max <- max(roi$y_min, roi$y_max, na.rm = TRUE)
    if (!all(is.finite(c(x_min, x_max, y_min, y_max)))) {
      warning("Skipping ROI with non-finite coordinates: ", roi_label, call. = FALSE)
      next
    }

    group_data <- labeled[labeled$roi_id == roi_label, , drop = FALSE]
    n_pixels <- nrow(group_data)
    if (n_pixels < min_pixels) next

    sample_id <- paste0("roi_", utils::URLencode(roi_label, reserved = TRUE))
    means <- colMeans(group_data[fcols], na.rm = TRUE)

    sample_rows[[sample_id]] <- data.frame(
      sample_id = sample_id,
      roi_label = roi_label,
      x_min = x_min,
      x_max = x_max,
      y_min = y_min,
      y_max = y_max,
      x_span = x_max - x_min,
      y_span = y_max - y_min,
      n_pixels = n_pixels,
      t(means),
      check.names = FALSE
    )
    mapping_rows[[sample_id]] <- data.frame(
      sample_id = sample_id,
      roi_label = roi_label,
      x_min = x_min,
      x_max = x_max,
      y_min = y_min,
      y_max = y_max,
      x_span = x_max - x_min,
      y_span = y_max - y_min,
      n_pixels = n_pixels,
      pixel_ids = paste(group_data$pixel_id, collapse = ","),
      stringsAsFactors = FALSE
    )

    annotated <- group_data[, c("pixel_id", "x", "y"), drop = FALSE]
    annotated$roi_label <- roi_label
    annotated$sample_id <- sample_id
    annotated_rows[[sample_id]] <- annotated
  }

  sample_matrix <- if (length(sample_rows)) do.call(rbind, sample_rows) else data.frame()
  sample_mapping <- if (length(mapping_rows)) do.call(rbind, mapping_rows) else data.frame()
  annotated_pixels <- if (length(annotated_rows)) do.call(rbind, annotated_rows) else data.frame()
  rownames(sample_matrix) <- NULL
  rownames(sample_mapping) <- NULL
  rownames(annotated_pixels) <- NULL

  if (nrow(sample_matrix) == 0) {
    stop("No coordinate ROI contained enough pixels. Check roi_table and min_pixels.", call. = FALSE)
  }

  list(sample_matrix = sample_matrix, sample_mapping = sample_mapping, annotated_pixels = annotated_pixels)
}

make_matched_sample_id <- function(section_id, matched_region_label) {
  # URL-encode both parts so sample_id can safely use "__" as a separator.
  # Use utils::URLdecode() before parsing these labels back downstream.
  paste(
    utils::URLencode(as.character(section_id), reserved = TRUE),
    utils::URLencode(as.character(matched_region_label), reserved = TRUE),
    sep = "__"
  )
}

sample_matched_regions <- function(pixel_matrix,
                                   region_column = "matched_region_label",
                                   section_column = "section_id",
                                   section_order_column = "section_order",
                                   min_pixels = 30) {
  required_columns(pixel_matrix, c("pixel_id", section_column, region_column), "Matched-region pixel matrix")
  fcols <- feature_columns(pixel_matrix)
  if (length(fcols) == 0) {
    stop("No mz feature columns available for matched-region sampling.", call. = FALSE)
  }

  work <- pixel_matrix[!is.na(pixel_matrix[[region_column]]) & pixel_matrix[[region_column]] != "", , drop = FALSE]
  if (nrow(work) == 0) {
    stop("No pixels have matched region labels.", call. = FALSE)
  }

  split_groups <- split(
    work,
    interaction(work[[section_column]], work[[region_column]], drop = TRUE),
    drop = TRUE
  )
  sample_rows <- list()
  mapping_rows <- list()

  for (group_name in names(split_groups)) {
    group_data <- split_groups[[group_name]]
    n_pixels <- nrow(group_data)
    if (n_pixels < min_pixels) next

    section_id <- as.character(group_data[[section_column]][1])
    matched_region_label <- as.character(group_data[[region_column]][1])
    section_order <- if (section_order_column %in% names(group_data)) {
      group_data[[section_order_column]][1]
    } else {
      NA
    }
    sample_id <- make_matched_sample_id(section_id, matched_region_label)
    means <- colMeans(group_data[fcols], na.rm = TRUE)

    sample_rows[[sample_id]] <- data.frame(
      sample_id = sample_id,
      section_id = section_id,
      section_order = section_order,
      matched_region_label = matched_region_label,
      n_pixels = n_pixels,
      t(means),
      check.names = FALSE
    )
    mapping_rows[[sample_id]] <- data.frame(
      sample_id = sample_id,
      section_id = section_id,
      section_order = section_order,
      matched_region_label = matched_region_label,
      n_pixels = n_pixels,
      pixel_ids = paste(group_data$pixel_id, collapse = ","),
      stringsAsFactors = FALSE
    )
  }

  sample_matrix <- if (length(sample_rows)) do.call(rbind, sample_rows) else data.frame()
  sample_mapping <- if (length(mapping_rows)) do.call(rbind, mapping_rows) else data.frame()
  rownames(sample_matrix) <- NULL
  rownames(sample_mapping) <- NULL
  list(sample_matrix = sample_matrix, sample_mapping = sample_mapping, annotated_pixels = work)
}

make_metaboanalyst_data <- function(sample_matrix,
                                    group_column = NULL,
                                    group_prefix = NULL) {
  required_columns(sample_matrix, "sample_id", "Sample matrix")
  fcols <- feature_columns(sample_matrix)
  if (length(fcols) == 0) {
    stop("No mz_ feature columns found in sample_matrix.", call. = FALSE)
  }
  if (is.null(group_column)) {
    group_column <- if ("cluster" %in% names(sample_matrix)) {
      "cluster"
    } else if ("roi_id" %in% names(sample_matrix)) {
      "roi_id"
    } else if ("roi_label" %in% names(sample_matrix)) {
      "roi_label"
    } else if ("matched_region_label" %in% names(sample_matrix)) {
      "matched_region_label"
    } else {
      NULL
    }
  }
  if (is.null(group_column) || !group_column %in% names(sample_matrix)) {
    stop("Sample matrix must contain a group column such as cluster or matched_region_label.", call. = FALSE)
  }

  group_values <- sample_matrix[[group_column]]
  if (is.null(group_prefix)) {
    group_prefix <- if (identical(group_column, "cluster")) "Region_" else ""
  }
  data.frame(
    Sample = sample_matrix$sample_id,
    Group = paste0(group_prefix, group_values),
    sample_matrix[fcols],
    check.names = FALSE
  )
}

differential_region_analysis <- function(sample_matrix,
                                         group_column = "roi_id",
                                         subject_column = NULL,
                                         section_column = NULL,
                                         reference_group = NULL,
                                         p_adjust_method = "BH") {
  required_columns(sample_matrix, c("sample_id", group_column), "Sample matrix")
  fcols <- feature_columns(sample_matrix)
  if (length(fcols) == 0) {
    stop("No mz_ feature columns found in sample_matrix.", call. = FALSE)
  }

  optional_columns <- c(subject_column, section_column)
  optional_columns <- optional_columns[!is.na(optional_columns) & nzchar(optional_columns)]
  required_columns(sample_matrix, optional_columns, "Sample matrix")
  groups <- unique(as.character(sample_matrix[[group_column]]))
  groups <- groups[!is.na(groups) & nzchar(groups)]
  if (length(groups) < 2L) {
    stop("At least two non-empty regions are required for differential analysis.", call. = FALSE)
  }

  if (!is.null(reference_group)) {
    reference_group <- as.character(reference_group)[1]
    if (!reference_group %in% groups) {
      stop("reference_group is not present in ", group_column, ".", call. = FALSE)
    }
    comparison_groups <- setdiff(groups, reference_group)
    comparisons <- rbind(rep(reference_group, length(comparison_groups)), comparison_groups)
  } else {
    comparisons <- utils::combn(groups, 2L)
  }

  if (!is.null(subject_column)) {
    replicate_column <- subject_column
    inference_unit <- "biological_subject"
    pseudoreplication_warning <- FALSE
  } else if (!is.null(section_column)) {
    replicate_column <- section_column
    inference_unit <- "section"
    pseudoreplication_warning <- TRUE
    warning(
      "No biological subject column was supplied. Tiles are aggregated within section and ",
      "section-level results must not be described as independent biological replication.",
      call. = FALSE
    )
  } else {
    replicate_column <- "sample_id"
    inference_unit <- "sample_or_tile"
    pseudoreplication_warning <- TRUE
    warning(
      "No biological subject or section column was supplied. P-values treat sample rows as ",
      "independent and may be pseudoreplicated; use results as exploratory only.",
      call. = FALSE
    )
  }

  work_columns <- unique(c(group_column, replicate_column, section_column, fcols))
  work <- sample_matrix[, work_columns, drop = FALSE]
  work[[group_column]] <- as.character(work[[group_column]])
  work[[replicate_column]] <- as.character(work[[replicate_column]])
  work <- work[
    !is.na(work[[group_column]]) & nzchar(work[[group_column]]) &
      !is.na(work[[replicate_column]]) & nzchar(work[[replicate_column]]),
    ,
    drop = FALSE
  ]
  if (nrow(work) == 0L) {
    stop("No rows remain after removing missing group or replicate identifiers.", call. = FALSE)
  }

  aggregate_data <- work
  aggregate_data$.__group__ <- aggregate_data[[group_column]]
  aggregate_data$.__replicate__ <- aggregate_data[[replicate_column]]
  aggregated <- stats::aggregate(
    aggregate_data[, fcols, drop = FALSE],
    by = list(
      .__group__ = aggregate_data$.__group__,
      .__replicate__ = aggregate_data$.__replicate__
    ),
    FUN = function(x) mean(as.numeric(x), na.rm = TRUE)
  )

  standardized_effect <- function(left, right) {
    n_left <- sum(is.finite(left))
    n_right <- sum(is.finite(right))
    if (n_left < 2L || n_right < 2L) return(NA_real_)
    pooled_df <- n_left + n_right - 2
    pooled_sd <- sqrt(((n_left - 1) * stats::var(left) + (n_right - 1) * stats::var(right)) / pooled_df)
    if (!is.finite(pooled_sd) || pooled_sd == 0) return(NA_real_)
    d <- (mean(right) - mean(left)) / pooled_sd
    correction <- if (pooled_df > 1) 1 - 3 / (4 * pooled_df - 1) else NA_real_
    d * correction
  }

  rows <- list()
  row_index <- 1L
  for (comparison_index in seq_len(ncol(comparisons))) {
    group_a <- comparisons[1, comparison_index]
    group_b <- comparisons[2, comparison_index]
    left_rows <- aggregated$.__group__ == group_a
    right_rows <- aggregated$.__group__ == group_b
    for (feature in fcols) {
      left <- as.numeric(aggregated[[feature]][left_rows])
      right <- as.numeric(aggregated[[feature]][right_rows])
      names(left) <- aggregated$.__replicate__[left_rows]
      names(right) <- aggregated$.__replicate__[right_rows]
      left <- left[is.finite(left)]
      right <- right[is.finite(right)]
      shared_replicates <- intersect(names(left), names(right))
      paired <- inference_unit != "sample_or_tile" && length(shared_replicates) >= 2L
      p_value <- NA_real_
      if (paired) {
        left_test <- left[shared_replicates]
        right_test <- right[shared_replicates]
      } else {
        left_test <- left
        right_test <- right
      }
      if (length(left_test) >= 2L && length(right_test) >= 2L) {
        if (paired && stats::sd(right_test - left_test) == 0) {
          p_value <- if (mean(right_test - left_test) == 0) 1 else 0
        } else if (!paired && stats::sd(left_test) == 0 && stats::sd(right_test) == 0) {
          p_value <- if (mean(right_test) == mean(left_test)) 1 else 0
        } else {
          p_value <- tryCatch(
            stats::t.test(right_test, left_test, paired = paired)$p.value,
            error = function(e) NA_real_
          )
        }
      }
      paired_effect <- if (paired) {
        differences <- right_test - left_test
        difference_sd <- stats::sd(differences)
        if (is.finite(difference_sd) && difference_sd > 0) mean(differences) / difference_sd else NA_real_
      } else {
        standardized_effect(left_test, right_test)
      }
      rows[[row_index]] <- data.frame(
        feature = feature,
        mz = suppressWarnings(as.numeric(sub("^mz_", "", feature))),
        group_a = group_a,
        group_b = group_b,
        mean_group_a = if (length(left)) mean(left) else NA_real_,
        mean_group_b = if (length(right)) mean(right) else NA_real_,
        effect_size = if (length(left) && length(right)) mean(right) - mean(left) else NA_real_,
        standardized_effect = paired_effect,
        p_value = p_value,
        n_group_a = length(left),
        n_group_b = length(right),
        n_pairs = if (paired) length(shared_replicates) else NA_integer_,
        test_type = if (paired) "paired_t_test" else "welch_t_test",
        inference_unit = inference_unit,
        subject_column = if (is.null(subject_column)) NA_character_ else subject_column,
        section_column = if (is.null(section_column)) NA_character_ else section_column,
        pseudoreplication_warning = pseudoreplication_warning,
        stringsAsFactors = FALSE
      )
      row_index <- row_index + 1L
    }
  }
  result <- do.call(rbind, rows)
  result$fdr <- stats::ave(
    result$p_value,
    interaction(result$group_a, result$group_b, drop = TRUE),
    FUN = function(x) stats::p.adjust(x, method = p_adjust_method)
  )
  result <- result[, c(
    "feature", "mz", "group_a", "group_b", "mean_group_a", "mean_group_b",
    "effect_size", "standardized_effect", "p_value", "fdr", "n_group_a", "n_group_b", "n_pairs", "test_type",
    "inference_unit", "subject_column", "section_column", "pseudoreplication_warning"
  )]
  result[order(result$fdr, -abs(result$effect_size), na.last = TRUE), , drop = FALSE]
}

differential_region_analysis_wilcoxon <- function(
    sample_matrix,
    group_column = "roi_id",
    subject_column = NULL,
    section_column = NULL,
    reference_group = NULL,
    aggregate_fun = c("mean", "median"),
    exact_method = c("auto", "exact", "asymptotic"),
    min_replicates = 5L,
    confidence_level = 0.95,
    p_adjust_method = "BH") {
  aggregate_fun <- match.arg(aggregate_fun)
  exact_method <- match.arg(exact_method)
  required_columns(sample_matrix, c("sample_id", group_column), "Sample matrix")
  fcols <- feature_columns(sample_matrix)
  if (!length(fcols)) {
    stop("No mz_ feature columns found in sample_matrix.", call. = FALSE)
  }
  if (any(!vapply(sample_matrix[fcols], is.numeric, logical(1)))) {
    stop("All mz_ feature columns must be numeric.", call. = FALSE)
  }
  if (length(min_replicates) != 1L || !is.finite(min_replicates) ||
      min_replicates < 3L || min_replicates != as.integer(min_replicates)) {
    stop("min_replicates must be one integer of at least three.", call. = FALSE)
  }
  min_replicates <- as.integer(min_replicates)
  if (length(confidence_level) != 1L || !is.finite(confidence_level) ||
      confidence_level <= 0 || confidence_level >= 1) {
    stop("confidence_level must lie strictly between zero and one.", call. = FALSE)
  }
  optional_columns <- c(subject_column, section_column)
  optional_columns <- optional_columns[!is.na(optional_columns) & nzchar(optional_columns)]
  required_columns(sample_matrix, optional_columns, "Sample matrix")

  groups <- sort(unique(as.character(sample_matrix[[group_column]])))
  groups <- groups[!is.na(groups) & nzchar(groups)]
  if (length(groups) < 2L) {
    stop("At least two non-empty regions are required for differential analysis.", call. = FALSE)
  }
  if (!is.null(reference_group)) {
    reference_group <- as.character(reference_group)[1]
    if (!reference_group %in% groups) {
      stop("reference_group is not present in ", group_column, ".", call. = FALSE)
    }
    comparison_groups <- setdiff(groups, reference_group)
    comparisons <- rbind(rep(reference_group, length(comparison_groups)), comparison_groups)
  } else {
    comparisons <- utils::combn(groups, 2L)
  }

  if (!is.null(subject_column)) {
    replicate_column <- subject_column
    inference_unit <- "biological_subject"
    paired_design <- TRUE
    pseudoreplication_warning <- FALSE
  } else if (!is.null(section_column)) {
    replicate_column <- section_column
    inference_unit <- "section"
    paired_design <- TRUE
    pseudoreplication_warning <- TRUE
    warning(
      "No biological subject column was supplied. Wilcoxon results are aggregated by section ",
      "and must not be described as independent biological replication.",
      call. = FALSE
    )
  } else {
    replicate_column <- "sample_id"
    inference_unit <- "sample_or_tile"
    paired_design <- FALSE
    pseudoreplication_warning <- TRUE
    warning(
      "No biological subject or section column was supplied. Wilcoxon results treat sample ",
      "rows as independent and may be pseudoreplicated; use them as exploratory only.",
      call. = FALSE
    )
  }

  work_columns <- unique(c(group_column, replicate_column, fcols))
  work <- sample_matrix[, work_columns, drop = FALSE]
  work[[group_column]] <- as.character(work[[group_column]])
  work[[replicate_column]] <- as.character(work[[replicate_column]])
  work <- work[
    !is.na(work[[group_column]]) & nzchar(work[[group_column]]) &
      !is.na(work[[replicate_column]]) & nzchar(work[[replicate_column]]),
    , drop = FALSE
  ]
  if (!nrow(work)) {
    stop("No rows remain after removing missing group or replicate identifiers.", call. = FALSE)
  }
  summarize_replicate <- switch(
    aggregate_fun,
    mean = function(x) {
      x <- as.numeric(x)
      if (all(!is.finite(x))) NA_real_ else mean(x[is.finite(x)])
    },
    median = function(x) {
      x <- as.numeric(x)
      if (all(!is.finite(x))) NA_real_ else stats::median(x[is.finite(x)])
    }
  )
  aggregate_data <- work
  aggregate_data$.__group__ <- aggregate_data[[group_column]]
  aggregate_data$.__replicate__ <- aggregate_data[[replicate_column]]
  aggregated <- stats::aggregate(
    aggregate_data[, fcols, drop = FALSE],
    by = list(
      .__group__ = aggregate_data$.__group__,
      .__replicate__ = aggregate_data$.__replicate__
    ),
    FUN = summarize_replicate
  )

  rows <- list()
  row_index <- 1L
  for (comparison_index in seq_len(ncol(comparisons))) {
    group_a <- comparisons[1, comparison_index]
    group_b <- comparisons[2, comparison_index]
    left_rows <- aggregated$.__group__ == group_a
    right_rows <- aggregated$.__group__ == group_b
    left_ids <- as.character(aggregated$.__replicate__[left_rows])
    right_ids <- as.character(aggregated$.__replicate__[right_rows])
    shared_replicates <- sort(intersect(left_ids, right_ids))
    contrast_eligible <- if (paired_design) {
      length(shared_replicates) >= min_replicates
    } else {
      length(left_ids) >= min_replicates && length(right_ids) >= min_replicates
    }
    contrast_skip_reason <- if (contrast_eligible) {
      NA_character_
    } else if (paired_design) {
      "too_few_shared_replicates"
    } else {
      "too_few_replicates_per_group"
    }
    pair_coverage <- if (paired_design) {
      length(shared_replicates) / length(unique(c(left_ids, right_ids)))
    } else {
      NA_real_
    }

    for (feature in fcols) {
      left <- as.numeric(aggregated[[feature]][left_rows])
      right <- as.numeric(aggregated[[feature]][right_rows])
      names(left) <- left_ids
      names(right) <- right_ids
      if (paired_design) {
        left_test <- left[shared_replicates]
        right_test <- right[shared_replicates]
        complete <- is.finite(left_test) & is.finite(right_test)
        left_test <- left_test[complete]
        right_test <- right_test[complete]
      } else {
        left_test <- left[is.finite(left)]
        right_test <- right[is.finite(right)]
      }
      differences <- if (paired_design) right_test - left_test else numeric()
      n_zero_differences <- if (paired_design) sum(differences == 0) else NA_integer_
      n_ties <- if (paired_design) {
        sum(duplicated(abs(differences[differences != 0])))
      } else {
        sum(duplicated(c(left_test, right_test)))
      }
      enough_feature_data <- if (paired_design) {
        length(left_test) >= min_replicates
      } else {
        length(left_test) >= min_replicates && length(right_test) >= min_replicates
      }
      nonzero_information <- if (paired_design) {
        any(differences != 0)
      } else {
        length(unique(c(left_test, right_test))) > 1L
      }
      exact_eligible <- if (paired_design) {
        length(differences) < 50L && n_zero_differences == 0L && n_ties == 0L
      } else {
        length(left_test) + length(right_test) < 50L && n_ties == 0L
      }
      use_exact <- switch(
        exact_method,
        auto = exact_eligible,
        exact = exact_eligible,
        asymptotic = FALSE
      )
      status <- "fitted"
      skip_reason <- NA_character_
      if (!contrast_eligible) {
        status <- "skipped"
        skip_reason <- contrast_skip_reason
      } else if (!enough_feature_data) {
        status <- "skipped"
        skip_reason <- "feature_missingness_reduces_replication"
      } else if (!nonzero_information) {
        status <- "constant"
        skip_reason <- "no_rank_information"
      } else if (exact_method == "exact" && !exact_eligible) {
        status <- "skipped"
        skip_reason <- "exact_test_unavailable_with_ties_zeros_or_large_sample"
      }

      test <- NULL
      warning_messages <- character()
      if (status == "fitted") {
        test <- withCallingHandlers(
          tryCatch(
            stats::wilcox.test(
              right_test, left_test,
              paired = paired_design,
              exact = use_exact,
              correct = !use_exact,
              conf.int = TRUE,
              conf.level = confidence_level
            ),
            error = function(error) error
          ),
          warning = function(warning) {
            warning_messages <<- unique(c(warning_messages, conditionMessage(warning)))
            invokeRestart("muffleWarning")
          }
        )
        if (inherits(test, "error")) {
          status <- "failed"
          skip_reason <- conditionMessage(test)
          test <- NULL
        }
      }
      if (status == "constant") {
        p_value <- 1
        statistic <- 0
      } else {
        p_value <- if (is.null(test)) NA_real_ else unname(test$p.value)
        statistic <- if (is.null(test)) NA_real_ else unname(test$statistic)
      }
      rank_biserial <- NA_real_
      if (!is.null(test)) {
        if (paired_design) {
          nonzero <- differences[differences != 0]
          ranks <- rank(abs(nonzero), ties.method = "average")
          denominator <- sum(ranks)
          rank_biserial <- if (denominator > 0) {
            (sum(ranks[nonzero > 0]) - sum(ranks[nonzero < 0])) / denominator
          } else NA_real_
        } else {
          rank_biserial <- 2 * statistic / (length(right_test) * length(left_test)) - 1
        }
      } else if (status == "constant") {
        rank_biserial <- 0
      }
      estimate <- if (is.null(test) || is.null(test$estimate)) NA_real_ else unname(test$estimate[1])
      confidence_interval <- if (is.null(test) || is.null(test$conf.int)) {
        c(NA_real_, NA_real_)
      } else {
        unname(test$conf.int[1:2])
      }
      median_left <- if (length(left_test)) stats::median(left_test) else NA_real_
      median_right <- if (length(right_test)) stats::median(right_test) else NA_real_

      rows[[row_index]] <- data.frame(
        feature = feature,
        mz = suppressWarnings(as.numeric(sub("^mz_", "", feature))),
        group_a = group_a,
        group_b = group_b,
        median_group_a = median_left,
        median_group_b = median_right,
        median_difference = median_right - median_left,
        median_paired_difference = if (paired_design && length(differences)) stats::median(differences) else NA_real_,
        hodges_lehmann_shift = estimate,
        confidence_lower = confidence_interval[1],
        confidence_upper = confidence_interval[2],
        rank_biserial_correlation = rank_biserial,
        statistic = statistic,
        p_value = p_value,
        n_group_a = length(left_test),
        n_group_b = length(right_test),
        n_pairs = if (paired_design) length(left_test) else NA_integer_,
        n_shared_replicates_design = if (paired_design) length(shared_replicates) else NA_integer_,
        pair_coverage = pair_coverage,
        n_ties = n_ties,
        n_zero_differences = n_zero_differences,
        test_type = if (paired_design) "paired_wilcoxon_signed_rank" else "wilcoxon_rank_sum",
        p_value_method = if (status %in% c("fitted", "constant")) {
          if (use_exact) "exact" else "asymptotic"
        } else NA_character_,
        status = status,
        skip_reason = skip_reason,
        warning_message = paste(warning_messages, collapse = " | "),
        inference_unit = inference_unit,
        subject_column = if (is.null(subject_column)) NA_character_ else subject_column,
        section_column = if (is.null(section_column)) NA_character_ else section_column,
        aggregation = aggregate_fun,
        pseudoreplication_warning = pseudoreplication_warning,
        stringsAsFactors = FALSE
      )
      row_index <- row_index + 1L
    }
  }
  result <- do.call(rbind, rows)
  result$fdr <- stats::ave(
    result$p_value,
    interaction(result$group_a, result$group_b, drop = TRUE),
    FUN = function(p) stats::p.adjust(p, method = p_adjust_method)
  )
  result$fdr_scope <- "within_contrast_across_features"
  result <- result[order(result$fdr, -abs(result$rank_biserial_correlation), na.last = TRUE), , drop = FALSE]
  rownames(result) <- NULL
  attr(result, "interpretation") <- paste(
    "Wilcoxon results are sensitivity analyses, not independent validation of parametric or spatial models.",
    "median_difference is descriptive; rank_biserial_correlation is the rank-based effect size,",
    "and hodges_lehmann_shift assumes a location-shift interpretation."
  )
  result
}

detect_metaboanalyst_result <- function(result_data) {
  columns <- names(result_data)
  if (all(c("Feature", "VIP") %in% columns)) return("vip")
  if (all(c("Feature", "p.value") %in% columns) || all(c("Feature", "FDR") %in% columns)) return("differential")
  if ("Sample" %in% columns && any(grepl("^PC[0-9]+$", columns))) return("pca_scores")
  if ("Sample" %in% columns && any(grepl("^Comp[0-9]+$", columns))) return("plsda_scores")
  "unknown"
}

normalize_metaboanalyst_result <- function(result_data, source_name = "") {
  if (is.null(result_data) || ncol(result_data) == 0) return(result_data)

  original_names <- names(result_data)
  clean_names <- tolower(gsub("[^a-z0-9]+", "", original_names))
  if (is.null(source_name) || length(source_name) == 0) source_name <- ""
  source_clean <- tolower(source_name)

  sample_candidates <- which(clean_names %in% c("sample", "samples", "sampleid", "samplename"))
  if (length(sample_candidates) > 0) {
    names(result_data)[sample_candidates[1]] <- "Sample"
  }

  has_score_columns <- any(grepl("^(PC|Comp)[0-9]+$", names(result_data)))
  feature_candidates <- which(clean_names %in% c(
    "feature", "features", "metabolite", "metabolites", "name", "variable", "variables", "compound", "compounds"
  ))
  if (length(feature_candidates) == 0 && !has_score_columns) {
    character_cols <- which(vapply(result_data, function(x) is.character(x) || is.factor(x), logical(1)))
    character_cols <- setdiff(character_cols, match("Sample", names(result_data)))
    feature_candidates <- character_cols
  }
  if (length(feature_candidates) > 0) {
    names(result_data)[feature_candidates[1]] <- "Feature"
  }

  vip_candidates <- which(grepl("vip", clean_names))
  if (length(vip_candidates) == 0 && grepl("vip", source_clean) && "Feature" %in% names(result_data)) {
    numeric_cols <- which(vapply(result_data, is.numeric, logical(1)))
    numeric_cols <- setdiff(numeric_cols, match("Feature", names(result_data)))
    vip_candidates <- numeric_cols
  }
  if (length(vip_candidates) > 0) {
    names(result_data)[vip_candidates[1]] <- "VIP"
  }

  p_candidates <- which(clean_names %in% c("pvalue", "pval", "p"))
  if (length(p_candidates) > 0) names(result_data)[p_candidates[1]] <- "p.value"

  fdr_candidates <- which(clean_names %in% c("fdr", "adjp", "padj", "qvalue"))
  if (length(fdr_candidates) > 0) names(result_data)[fdr_candidates[1]] <- "FDR"

  result_data
}

first_existing_column <- function(data, candidates, label) {
  hit <- candidates[candidates %in% names(data)]
  if (length(hit) == 0) {
    stop(
      label, " is missing. Tried column(s): ",
      paste(candidates, collapse = ", "),
      call. = FALSE
    )
  }
  hit[1]
}

collapse_unique_values <- function(x) {
  x <- unique(as.character(x[!is.na(x) & x != ""]))
  if (length(x) == 0) return(NA_character_)
  paste(x, collapse = "; ")
}

match_rows_by_ppm <- function(query_mz, reference_mz, ppm) {
  if (!is.finite(query_mz)) return(integer())
  diff_ppm <- abs(query_mz - reference_mz) / reference_mz * 1e6
  which(is.finite(diff_ppm) & diff_ppm <= ppm)
}

parse_formula_counts <- function(formula) {
  if (is.na(formula) || formula == "") return(NULL)
  matches <- regmatches(formula, gregexpr("([A-Z][a-z]?)([0-9]*)", formula, perl = TRUE))[[1]]
  if (length(matches) == 0) return(NULL)

  elements <- sub("^([A-Z][a-z]?).*", "\\1", matches, perl = TRUE)
  counts <- sub("^[A-Z][a-z]?([0-9]*)$", "\\1", matches, perl = TRUE)
  counts <- ifelse(counts == "", "1", counts)
  counts <- suppressWarnings(as.integer(counts))
  if (any(!is.finite(counts))) return(NULL)
  stats::setNames(as.integer(tapply(counts, elements, sum)), sort(unique(elements)))
}

format_formula_counts <- function(counts) {
  if (is.null(counts) || length(counts) == 0) return(NA_character_)
  counts <- counts[counts > 0]
  order_names <- c("C", "H", setdiff(sort(names(counts)), c("C", "H")))
  order_names <- order_names[order_names %in% names(counts)]
  paste0(order_names, ifelse(counts[order_names] == 1L, "", counts[order_names]), collapse = "")
}

add_formula_hydrogen <- function(formula, n = 1L) {
  counts <- parse_formula_counts(formula)
  if (is.null(counts)) return(NA_character_)
  if (!"H" %in% names(counts)) counts <- c(counts, H = 0L)
  counts[["H"]] <- counts[["H"]] + as.integer(n)
  if (counts[["H"]] < 0L) return(NA_character_)
  format_formula_counts(counts)
}

neutralize_mminush_formulas <- function(formulas, enabled = TRUE) {
  formulas <- as.character(formulas)
  if (!isTRUE(enabled)) return(formulas)
  vapply(formulas, add_formula_hydrogen, character(1), n = 1L, USE.NAMES = FALSE)
}

classify_formula_overlap <- function(left_formula, right_formula) {
  left_formula <- unique(as.character(left_formula[!is.na(left_formula) & left_formula != ""]))
  right_formula <- unique(as.character(right_formula[!is.na(right_formula) & right_formula != ""]))
  shared_formula <- intersect(left_formula, right_formula)

  class <- if (length(left_formula) == 0 && length(right_formula) == 0) {
    "both_missing"
  } else if (length(left_formula) == 0 || length(right_formula) == 0) {
    "unilateral_missing"
  } else if (setequal(left_formula, right_formula)) {
    "exact_match"
  } else if (length(shared_formula) > 0) {
    "candidate_overlap"
  } else {
    "discordant"
  }

  list(class = class, shared = shared_formula)
}

best_formula_class <- function(...) {
  classes <- c(...)
  priority <- c(
    exact_match = 5L,
    candidate_overlap = 4L,
    discordant = 3L,
    unilateral_missing = 2L,
    both_missing = 1L
  )
  classes <- classes[classes %in% names(priority)]
  if (length(classes) == 0) return("both_missing")
  classes[which.max(priority[classes])]
}

cross_validate_annotations <- function(feature_csv,
                                       swisslipids_csv,
                                       metaspace_csv = NULL,
                                       ppm = 5,
                                       feature_id_col = "id",
                                       feature_mz_col = "mz",
                                       swiss_mz_col = "MZ_MminusH",
                                       swiss_formula_col = NULL,
                                       swiss_name_col = NULL,
                                       swiss_id_col = NULL,
                                       metaspace_mz_col = "mz",
                                       metaspace_formula_col = "formula",
                                       metaspace_name_col = "moleculeNames",
                                       metaspace_id_col = "moleculeIds",
                                       neutral_formula_correction = TRUE) {
  features <- utils::read.csv(feature_csv, check.names = FALSE, stringsAsFactors = FALSE)
  swiss <- utils::read.csv(swisslipids_csv, check.names = FALSE, stringsAsFactors = FALSE)
  metaspace <- if (!is.null(metaspace_csv)) {
    utils::read.csv(metaspace_csv, check.names = FALSE, stringsAsFactors = FALSE)
  } else {
    NULL
  }

  feature_id_col <- first_existing_column(features, feature_id_col, "Feature id column")
  feature_mz_col <- first_existing_column(features, feature_mz_col, "Feature m/z column")
  swiss_mz_col <- first_existing_column(swiss, swiss_mz_col, "SwissLipids m/z column")

  if (is.null(swiss_formula_col)) {
    swiss_formula_col <- first_existing_column(
      swiss,
      c("Formula", "formula", "MolecularFormula", "molecular_formula", "SumFormula", "sum_formula"),
      "SwissLipids formula column"
    )
  } else {
    swiss_formula_col <- first_existing_column(swiss, swiss_formula_col, "SwissLipids formula column")
  }
  if (is.null(swiss_name_col)) {
    swiss_name_col <- intersect(
      c("Name", "name", "LipidName", "lipid_name", "Species", "species", "Abbreviation", "abbreviation"),
      names(swiss)
    )
    swiss_name_col <- if (length(swiss_name_col) > 0) swiss_name_col[1] else NA_character_
  } else {
    swiss_name_col <- first_existing_column(swiss, swiss_name_col, "SwissLipids name column")
  }
  if (is.null(swiss_id_col)) {
    swiss_id_col <- intersect(
      c("LipidID", "lipid_id", "SwissLipids_ID", "SwissLipidsID", "LipidMaps_ID", "LipidMapsID", "ID", "id", "SLM_ID", "slm_id"),
      names(swiss)
    )
    swiss_id_col <- if (length(swiss_id_col) > 0) swiss_id_col[1] else NA_character_
  } else {
    swiss_id_col <- first_existing_column(swiss, swiss_id_col, "SwissLipids id column")
  }

  swiss_mz <- suppressWarnings(as.numeric(swiss[[swiss_mz_col]]))
  swiss_formula_all <- as.character(swiss[[swiss_formula_col]])
  swiss_neutral_formula_all <- neutralize_mminush_formulas(
    swiss_formula_all,
    enabled = neutral_formula_correction
  )
  feature_mz <- suppressWarnings(as.numeric(features[[feature_mz_col]]))

  if (!is.null(metaspace)) {
    metaspace_mz_col <- first_existing_column(metaspace, metaspace_mz_col, "METASPACE m/z column")
    metaspace_formula_col <- first_existing_column(metaspace, metaspace_formula_col, "METASPACE formula column")
    metaspace_name_col <- if (metaspace_name_col %in% names(metaspace)) metaspace_name_col else NA_character_
    metaspace_id_col <- if (metaspace_id_col %in% names(metaspace)) metaspace_id_col else NA_character_
    metaspace_mz <- suppressWarnings(as.numeric(metaspace[[metaspace_mz_col]]))
  } else {
    metaspace_mz <- numeric()
  }

  rows <- lapply(seq_len(nrow(features)), function(i) {
    mz <- feature_mz[i]
    swiss_idx <- match_rows_by_ppm(mz, swiss_mz, ppm)
    meta_idx <- if (!is.null(metaspace)) match_rows_by_ppm(mz, metaspace_mz, ppm) else integer()

    swiss_formula <- unique(as.character(swiss[[swiss_formula_col]][swiss_idx]))
    swiss_formula <- swiss_formula[!is.na(swiss_formula) & swiss_formula != ""]
    swiss_neutral_formula <- unique(as.character(swiss_neutral_formula_all[swiss_idx]))
    swiss_neutral_formula <- swiss_neutral_formula[!is.na(swiss_neutral_formula) & swiss_neutral_formula != ""]
    meta_formula <- if (!is.null(metaspace)) {
      unique(as.character(metaspace[[metaspace_formula_col]][meta_idx]))
    } else {
      character()
    }
    meta_formula <- meta_formula[!is.na(meta_formula) & meta_formula != ""]

    raw_overlap <- classify_formula_overlap(swiss_formula, meta_formula)
    hydrogen_overlap <- classify_formula_overlap(swiss_neutral_formula, meta_formula)
    corrected_class <- best_formula_class(raw_overlap$class, hydrogen_overlap$class)
    corrected_shared <- unique(c(raw_overlap$shared, hydrogen_overlap$shared))
    corrected_shared <- corrected_shared[!is.na(corrected_shared) & corrected_shared != ""]
    meta_swiss_idx <- which(swiss_formula_all %in% meta_formula | swiss_neutral_formula_all %in% meta_formula)

    data.frame(
      feature_id = features[[feature_id_col]][i],
      measured_mz = mz,
      ppm_tolerance_da = mz * ppm / 1e6,
      swiss_candidate_count = length(swiss_idx),
      swiss_formula_count = length(swiss_formula),
      swiss_formulas = collapse_unique_values(swiss_formula),
      swiss_neutral_formulas = collapse_unique_values(swiss_neutral_formula),
      swiss_names = if (!is.na(swiss_name_col)) collapse_unique_values(swiss[[swiss_name_col]][swiss_idx]) else NA_character_,
      swiss_ids = if (!is.na(swiss_id_col)) collapse_unique_values(swiss[[swiss_id_col]][swiss_idx]) else NA_character_,
      metaspace_candidate_count = length(meta_idx),
      metaspace_formula_count = length(meta_formula),
      metaspace_formulas = collapse_unique_values(meta_formula),
      metaspace_names = if (!is.null(metaspace) && !is.na(metaspace_name_col)) collapse_unique_values(metaspace[[metaspace_name_col]][meta_idx]) else NA_character_,
      metaspace_ids = if (!is.null(metaspace) && !is.na(metaspace_id_col)) collapse_unique_values(metaspace[[metaspace_id_col]][meta_idx]) else NA_character_,
      metaspace_swiss_names = if (!is.na(swiss_name_col)) collapse_unique_values(swiss[[swiss_name_col]][meta_swiss_idx]) else NA_character_,
      metaspace_swiss_ids = if (!is.na(swiss_id_col)) collapse_unique_values(swiss[[swiss_id_col]][meta_swiss_idx]) else NA_character_,
      raw_shared_formulas = collapse_unique_values(raw_overlap$shared),
      shared_formulas = collapse_unique_values(corrected_shared),
      raw_consistency_class = raw_overlap$class,
      consistency_class = corrected_class,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  })

  result <- do.call(rbind, rows)
  summary <- as.data.frame(table(result$consistency_class), stringsAsFactors = FALSE)
  names(summary) <- c("consistency_class", "feature_count")
  attr(result, "summary") <- summary
  result
}

run_spatial_metabolomics_pipeline <- function(msi_csv,
                                              x_col = "x",
                                              y_col = "y",
                                              bad_pixel_filter = TRUE,
                                              min_nonzero_count = 1,
                                              tic_normalize = TRUE,
                                              do_log = TRUE,
                                              do_scale = FALSE) {
  if (!file.exists(msi_csv)) {
    stop("MSI CSV file does not exist: ", msi_csv, call. = FALSE)
  }

  pixel_matrix <- import_peakpicked_msi_csv(msi_csv)
  required_columns(pixel_matrix, c(x_col, y_col), "MSI pixel matrix")
  fcols <- feature_columns(pixel_matrix)
  if (length(fcols) < 5) {
    stop("MSI data must contain at least 5 mz_ feature columns.", call. = FALSE)
  }

  coords <- pixel_matrix[, c(x_col, y_col), drop = FALSE]
  if (!all(is.finite(coords[[x_col]]) & is.finite(coords[[y_col]]))) {
    stop("MSI coordinates must be finite.", call. = FALSE)
  }

  values <- as.matrix(pixel_matrix[fcols])
  if (!all(is.finite(values))) {
    values[!is.finite(values)] <- 0
  }

  keep <- rep(TRUE, nrow(pixel_matrix))
  if (isTRUE(bad_pixel_filter)) {
    keep <- rowSums(values > 0, na.rm = TRUE) >= min_nonzero_count
  }
  filtered <- pixel_matrix[keep, , drop = FALSE]
  qc_summary <- list(
    raw_pixels = nrow(pixel_matrix),
    retained_pixels = nrow(filtered),
    removed_pixels = nrow(pixel_matrix) - nrow(filtered),
    input_feature_count = length(fcols),
    output_feature_count = length(fcols),
    feature_count = length(fcols)
  )

  if (nrow(filtered) == 0) {
    stop("No pixels remain after QC filtering.", call. = FALSE)
  }

  feature_meta <- data.frame(
    feature_id = seq_along(fcols),
    column_name = fcols,
    mz = suppressWarnings(as.numeric(sub("^mz_", "", fcols))),
    stringsAsFactors = FALSE
  )

  normalized <- as.matrix(filtered[fcols])
  if (isTRUE(tic_normalize)) {
    row_totals <- rowSums(normalized, na.rm = TRUE)
    normalized <- sweep(normalized, 1, pmax(row_totals, .Machine$double.eps), "/")
  }
  if (isTRUE(do_log)) {
    normalized <- log10(normalized + 1)
  }
  if (isTRUE(do_scale) && ncol(normalized) > 1) {
    normalized <- scale(normalized, center = TRUE, scale = TRUE)
    normalized[!is.finite(normalized)] <- 0
  }

  pixel_feature_matrix <- as.data.frame(normalized, stringsAsFactors = FALSE)
  pixel_feature_matrix$pixel_id <- filtered$pixel_id
  pixel_feature_matrix <- pixel_feature_matrix[, c("pixel_id", fcols), drop = FALSE]
  coordinates <- filtered[, c("pixel_id", x_col, y_col), drop = FALSE]
  if ("section_id" %in% names(filtered)) {
    coordinates$section_id <- filtered$section_id
  }
  coordinates <- coordinates[, intersect(c("pixel_id", x_col, y_col, "section_id"), names(coordinates)), drop = FALSE]

  list(
    pixel_feature_matrix = pixel_feature_matrix,
    coordinates = coordinates,
    feature_metadata = feature_meta,
    qc_summary = qc_summary,
    parameters = list(
      x_col = x_col,
      y_col = y_col,
      bad_pixel_filter = bad_pixel_filter,
      min_nonzero_count = min_nonzero_count,
      tic_normalize = tic_normalize,
      do_log = do_log,
      do_scale = do_scale
    )
  )
}

compute_spatially_variable_metabolites <- function(pixel_matrix,
                                                   coordinates = NULL,
                                                   x_col = "x",
                                                   y_col = "y",
                                                   fcols = feature_columns(pixel_matrix),
                                                   n_perm = 199,
                                                   alternative = c("greater", "two.sided"),
                                                   p_adjust_method = "BH",
                                                   seed = NULL,
                                                   neighbors = NULL,
                                                   neighbor_method = c("rook", "queen", "distance"),
                                                   distance_threshold = NULL,
                                                   weights = "binary",
                                                   symmetric = TRUE) {
  alternative <- match.arg(alternative)
  neighbor_method <- match.arg(neighbor_method)
  if (length(fcols) == 0) {
    stop("No mz_ features found for spatial analysis.", call. = FALSE)
  }

  if (!is.null(coordinates)) {
    required_columns(pixel_matrix, "pixel_id", "Pixel matrix")
    required_columns(coordinates, c("pixel_id", x_col, y_col), "Coordinates")
    coords <- merge(
      pixel_matrix[, "pixel_id", drop = FALSE],
      coordinates[, c("pixel_id", x_col, y_col), drop = FALSE],
      by = "pixel_id",
      sort = FALSE
    )
    if (nrow(coords) != nrow(pixel_matrix)) {
      stop("Coordinates must contain the same pixel_id values as pixel_matrix.", call. = FALSE)
    }
    coords <- coords[match(pixel_matrix$pixel_id, coords$pixel_id), , drop = FALSE]
  } else {
    required_columns(pixel_matrix, c("pixel_id", x_col, y_col), "Pixel matrix")
    coords <- pixel_matrix[, c("pixel_id", x_col, y_col), drop = FALSE]
  }

  if (is.null(neighbors)) {
    neighbors <- build_spatial_neighbors(
      coords[[x_col]], coords[[y_col]], method = neighbor_method,
      distance_threshold = distance_threshold, weights = weights, symmetric = symmetric
    )
  }

  if (!is.null(seed)) set.seed(seed)
  rows <- lapply(fcols, function(feature) {
    stat <- compute_morans_i_grid(
      values = pixel_matrix[[feature]],
      x = coords[[x_col]],
      y = coords[[y_col]],
      n_perm = n_perm,
      alternative = alternative,
      seed = NULL,
      neighbors = neighbors
    )
    data.frame(
      feature = feature,
      morans_i = stat$I,
      p_value = stat$p_value,
      n = stat$n,
      n_total = stat$n_total,
      n_effective = stat$n_effective,
      n_isolated = stat$n_isolated,
      n_edges = stat$n_edges,
      stringsAsFactors = FALSE
    )
  })
  result <- do.call(rbind, rows)
  result$adj_p_value <- stats::p.adjust(result$p_value, method = p_adjust_method)
  result[order(result$adj_p_value, -abs(result$morans_i)), , drop = FALSE]
}

compute_spatially_variable_metabolites_geary <- function(
    pixel_matrix,
    coordinates = NULL,
    x_col = "x",
    y_col = "y",
    fcols = feature_columns(pixel_matrix),
    n_perm = 199,
    alternative = c("less", "greater", "two.sided"),
    p_adjust_method = "BH",
    seed = NULL,
    neighbors = NULL,
    neighbor_method = c("rook", "queen", "distance"),
    distance_threshold = NULL,
    weights = "binary",
    symmetric = TRUE) {
  alternative <- match.arg(alternative)
  neighbor_method <- match.arg(neighbor_method)
  if (!length(fcols)) stop("No mz_ features found for spatial analysis.", call. = FALSE)
  missing_features <- setdiff(fcols, names(pixel_matrix))
  if (length(missing_features)) {
    stop("Requested spatial feature(s) are absent: ", paste(missing_features, collapse = ", "), call. = FALSE)
  }
  required_columns(pixel_matrix, "pixel_id", "Pixel matrix")
  if (anyNA(pixel_matrix$pixel_id) || anyDuplicated(pixel_matrix$pixel_id)) {
    stop("Pixel matrix pixel_id values must be non-missing and unique.", call. = FALSE)
  }
  if (!is.null(coordinates)) {
    required_columns(coordinates, c("pixel_id", x_col, y_col), "Coordinates")
    if (anyNA(coordinates$pixel_id) || anyDuplicated(coordinates$pixel_id)) {
      stop("Coordinate pixel_id values must be non-missing and unique.", call. = FALSE)
    }
    if (nrow(coordinates) != nrow(pixel_matrix) ||
        !setequal(coordinates$pixel_id, pixel_matrix$pixel_id)) {
      stop("Coordinates must contain exactly the pixel_id values in pixel_matrix.", call. = FALSE)
    }
    coordinate_index <- match(pixel_matrix$pixel_id, coordinates$pixel_id)
    coords <- coordinates[coordinate_index, c("pixel_id", x_col, y_col), drop = FALSE]
  } else {
    required_columns(pixel_matrix, c(x_col, y_col), "Pixel matrix")
    coords <- pixel_matrix[, c("pixel_id", x_col, y_col), drop = FALSE]
  }
  if (any(!is.finite(coords[[x_col]]) | !is.finite(coords[[y_col]]))) {
    stop("Spatial coordinates must be finite and non-missing.", call. = FALSE)
  }
  if (is.null(neighbors)) {
    neighbors <- build_spatial_neighbors(
      coords[[x_col]], coords[[y_col]], method = neighbor_method,
      distance_threshold = distance_threshold, weights = weights,
      symmetric = symmetric
    )
  } else {
    if (!inherits(neighbors, "spatial_neighbors")) {
      stop("neighbors must come from build_spatial_neighbors().", call. = FALSE)
    }
    if (!isTRUE(all.equal(as.numeric(coords[[x_col]]), neighbors$x)) ||
        !isTRUE(all.equal(as.numeric(coords[[y_col]]), neighbors$y))) {
      stop("neighbors coordinates and pixel_matrix row order do not match.", call. = FALSE)
    }
  }

  if (!is.null(seed)) set.seed(seed)
  rows <- lapply(fcols, function(feature) {
    stat <- compute_gearys_c_grid(
      values = pixel_matrix[[feature]], n_perm = n_perm,
      alternative = alternative, seed = NULL, neighbors = neighbors
    )
    data.frame(
      feature = feature, gearys_c = stat$C, expected_c = stat$expected_C,
      p_value = stat$p_value, alternative = stat$alternative,
      n = stat$n, n_total = stat$n_total, n_effective = stat$n_effective,
      n_isolated = stat$n_isolated, n_edges = stat$n_edges,
      stringsAsFactors = FALSE
    )
  })
  result <- do.call(rbind, rows)
  result$adj_p_value <- stats::p.adjust(result$p_value, method = p_adjust_method)
  rownames(result) <- NULL
  result[order(result$adj_p_value, -abs(result$gearys_c - 1), na.last = TRUE), , drop = FALSE]
}

compare_moran_geary <- function(moran_result,
                                geary_result,
                                alpha = 0.05) {
  required_columns(moran_result, c("feature", "morans_i", "adj_p_value"), "Moran result")
  required_columns(geary_result, c("feature", "gearys_c", "adj_p_value"), "Geary result")
  if (anyDuplicated(moran_result$feature) || anyDuplicated(geary_result$feature)) {
    stop("Each result must contain at most one row per feature.", call. = FALSE)
  }
  if (length(alpha) != 1L || !is.finite(alpha) || alpha <= 0 || alpha >= 1) {
    stop("alpha must be one finite value strictly between 0 and 1.", call. = FALSE)
  }
  merged <- merge(
    moran_result[, c("feature", "morans_i", "p_value", "adj_p_value"), drop = FALSE],
    geary_result[, c("feature", "gearys_c", "p_value", "adj_p_value"), drop = FALSE],
    by = "feature", suffixes = c("_moran", "_geary"), sort = FALSE
  )
  merged$moran_significant <- is.finite(merged$adj_p_value_moran) & merged$adj_p_value_moran < alpha
  merged$geary_significant <- is.finite(merged$adj_p_value_geary) & merged$adj_p_value_geary < alpha
  merged$agreement <- merged$moran_significant == merged$geary_significant
  merged$both_significant <- merged$moran_significant & merged$geary_significant
  merged$only_moran <- merged$moran_significant & !merged$geary_significant
  merged$only_geary <- merged$geary_significant & !merged$moran_significant
  merged$concordance_class <- ifelse(
    merged$both_significant, "both_significant",
    ifelse(merged$only_moran, "moran_only",
      ifelse(merged$only_geary, "geary_only", "neither_significant"))
  )
  merged$alpha <- alpha
  output <- merged[order(
    -merged$both_significant,
    pmin(merged$adj_p_value_moran, merged$adj_p_value_geary),
    na.last = TRUE
  ), , drop = FALSE]
  rownames(output) <- NULL
  attr(output, "interpretation") <- paste(
    "Moran and Geary results share the same measurements and spatial graph;",
    "agreement is a method-sensitivity result, not independent validation."
  )
  output
}

compute_binspect_feature <- function(values,
                                     neighbors,
                                     bin_method = c("kmeans", "rank"),
                                     percentage_rank = 30,
                                     inference = c("permutation", "fisher_approx"),
                                     n_perm = 199,
                                     kmeans_nstart = 10,
                                     kmeans_iter_max = 100,
                                     hub_min_neighbors = 1,
                                     seed = NULL) {
  bin_method <- match.arg(bin_method)
  inference <- match.arg(inference)
  if (!inherits(neighbors, "spatial_neighbors")) {
    stop("neighbors must come from build_spatial_neighbors().", call. = FALSE)
  }
  if (length(percentage_rank) != 1L || !is.finite(percentage_rank) ||
      percentage_rank <= 0 || percentage_rank >= 100) {
    stop("percentage_rank must be strictly between 0 and 100.", call. = FALSE)
  }
  if (length(n_perm) != 1L || !is.finite(n_perm) || n_perm < 0 || n_perm != as.integer(n_perm)) {
    stop("n_perm must be one non-negative integer.", call. = FALSE)
  }
  n_perm <- as.integer(n_perm)
  if (inference == "permutation" && n_perm < 1L) {
    stop("Permutation inference requires n_perm >= 1.", call. = FALSE)
  }
  if (length(hub_min_neighbors) != 1L || !is.finite(hub_min_neighbors) ||
      hub_min_neighbors < 1 || hub_min_neighbors != as.integer(hub_min_neighbors)) {
    stop("hub_min_neighbors must be one positive integer.", call. = FALSE)
  }
  values <- as.numeric(values)
  n_total <- length(values)
  if (n_total != length(neighbors$x)) {
    stop("values must contain one value per spatial-neighbor node.", call. = FALSE)
  }
  finite <- is.finite(values)
  if (sum(finite) < 3L) stop("At least three finite values are required for binSpect.", call. = FALSE)
  finite_values <- values[finite]
  if (length(unique(finite_values)) < 2L) {
    stop("Values are constant; binarization is undefined.", call. = FALSE)
  }
  if (!is.null(seed)) set.seed(seed)
  binary <- rep(NA_integer_, n_total)
  threshold <- NA_real_
  if (bin_method == "kmeans") {
    fit <- stats::kmeans(
      finite_values, centers = 2L, nstart = kmeans_nstart,
      iter.max = kmeans_iter_max, algorithm = "Lloyd"
    )
    cluster_means <- tapply(finite_values, fit$cluster, mean)
    high_cluster <- as.integer(names(cluster_means)[which.max(cluster_means)])
    binary[finite] <- as.integer(fit$cluster == high_cluster)
    threshold <- mean(range(cluster_means))
  } else {
    threshold <- unname(stats::quantile(
      finite_values, probs = 1 - percentage_rank / 100,
      na.rm = TRUE, names = FALSE, type = 7
    ))
    binary[finite] <- as.integer(finite_values >= threshold)
  }
  if (length(unique(binary[finite])) < 2L) {
    stop("Binarization did not produce both low and high states.", call. = FALSE)
  }

  edges <- neighbors$undirected_edges
  if (nrow(edges)) edges <- edges[finite[edges$from] & finite[edges$to], , drop = FALSE]
  if (!nrow(edges) || sum(edges$weight) <= 0) {
    stop("binSpect is undefined because the effective graph has no positive-weight edges.", call. = FALSE)
  }
  if (any(!is.finite(edges$weight)) || any(edges$weight < 0)) {
    stop("Spatial neighbor weights must be finite and non-negative.", call. = FALSE)
  }
  degree <- tabulate(c(edges$from, edges$to), nbins = n_total)
  effective <- finite & degree > 0L
  effective_nodes <- which(effective)
  effective_binary <- binary[effective_nodes]
  node_map <- integer(n_total)
  node_map[effective_nodes] <- seq_along(effective_nodes)
  edge_i <- node_map[edges$from]
  edge_j <- node_map[edges$to]
  edge_weight <- edges$weight

  symmetric_contingency <- function(state) {
    forward <- table(
      factor(state[edge_i], levels = c(0, 1)),
      factor(state[edge_j], levels = c(0, 1))
    )
    forward + t(forward)
  }
  contingency <- symmetric_contingency(effective_binary)
  observed_high_high <- sum(edge_weight * (effective_binary[edge_i] == 1L) *
    (effective_binary[edge_j] == 1L))
  odds_ratio <- if (contingency[1, 2] * contingency[2, 1] == 0) {
    if (contingency[1, 1] * contingency[2, 2] > 0) Inf else NA_real_
  } else {
    unname(contingency[1, 1] * contingency[2, 2] /
      (contingency[1, 2] * contingency[2, 1]))
  }
  if (inference == "permutation") {
    permuted <- replicate(n_perm, {
      perm <- sample(effective_binary, length(effective_binary), replace = FALSE)
      sum(edge_weight * (perm[edge_i] == 1L) * (perm[edge_j] == 1L))
    })
    p_value <- (sum(permuted >= observed_high_high) + 1) / (n_perm + 1)
  } else {
    # This reproduces the Fisher-table style used by binSpect, but edges share
    # nodes and both directions of each undirected edge are represented. The
    # resulting p-value is therefore an approximation, not an exact graph test.
    p_value <- stats::fisher.test(contingency, alternative = "greater")$p.value
  }

  high_nodes <- effective_nodes[effective_binary == 1L]
  high_neighbor_count <- vapply(high_nodes, function(node) {
    incident <- unique(c(
      edges$to[edges$from == node],
      edges$from[edges$to == node]
    ))
    sum(binary[incident] == 1L, na.rm = TRUE)
  }, integer(1))
  list(
    p_value = p_value, odds_ratio = odds_ratio,
    high_high_edge_weight = observed_high_high,
    contingency_table = contingency,
    n_high = sum(effective_binary == 1L),
    high_fraction = mean(effective_binary == 1L),
    n_effective = length(effective_nodes),
    n_isolated_or_nonfinite = n_total - length(effective_nodes),
    n_edges = nrow(edges),
    n_hub = sum(high_neighbor_count >= hub_min_neighbors),
    bin_method = bin_method, threshold = threshold,
    inference = inference, n_perm = if (inference == "permutation") n_perm else 0L,
    mean_expression = mean(finite_values),
    mean_high_expression = mean(values[high_nodes])
  )
}

compute_spatially_variable_metabolites_binspect <- function(
    pixel_matrix,
    coordinates = NULL,
    x_col = "x",
    y_col = "y",
    fcols = feature_columns(pixel_matrix),
    bin_method = c("kmeans", "rank"),
    percentage_rank = 30,
    inference = c("permutation", "fisher_approx"),
    n_perm = 199,
    p_adjust_method = "BH",
    seed = NULL,
    neighbors = NULL,
    neighbor_method = c("rook", "queen", "distance"),
    distance_threshold = NULL,
    weights = "binary",
    symmetric = TRUE) {
  bin_method <- match.arg(bin_method)
  inference <- match.arg(inference)
  neighbor_method <- match.arg(neighbor_method)
  if (!length(fcols)) stop("No mz_ features found for spatial analysis.", call. = FALSE)
  missing_features <- setdiff(fcols, names(pixel_matrix))
  if (length(missing_features)) {
    stop("Requested spatial feature(s) are absent: ", paste(missing_features, collapse = ", "), call. = FALSE)
  }
  required_columns(pixel_matrix, "pixel_id", "Pixel matrix")
  if (anyNA(pixel_matrix$pixel_id) || anyDuplicated(pixel_matrix$pixel_id)) {
    stop("Pixel matrix pixel_id values must be non-missing and unique.", call. = FALSE)
  }
  if (!is.null(coordinates)) {
    required_columns(coordinates, c("pixel_id", x_col, y_col), "Coordinates")
    if (anyNA(coordinates$pixel_id) || anyDuplicated(coordinates$pixel_id)) {
      stop("Coordinate pixel_id values must be non-missing and unique.", call. = FALSE)
    }
    if (nrow(coordinates) != nrow(pixel_matrix) ||
        !setequal(coordinates$pixel_id, pixel_matrix$pixel_id)) {
      stop("Coordinates must contain exactly the pixel_id values in pixel_matrix.", call. = FALSE)
    }
    coords <- coordinates[match(pixel_matrix$pixel_id, coordinates$pixel_id),
      c("pixel_id", x_col, y_col), drop = FALSE]
  } else {
    required_columns(pixel_matrix, c(x_col, y_col), "Pixel matrix")
    coords <- pixel_matrix[, c("pixel_id", x_col, y_col), drop = FALSE]
  }
  if (any(!is.finite(coords[[x_col]]) | !is.finite(coords[[y_col]]))) {
    stop("Spatial coordinates must be finite and non-missing.", call. = FALSE)
  }
  if (is.null(neighbors)) {
    neighbors <- build_spatial_neighbors(
      coords[[x_col]], coords[[y_col]], method = neighbor_method,
      distance_threshold = distance_threshold, weights = weights,
      symmetric = symmetric
    )
  } else {
    if (!inherits(neighbors, "spatial_neighbors")) {
      stop("neighbors must come from build_spatial_neighbors().", call. = FALSE)
    }
    if (!isTRUE(all.equal(as.numeric(coords[[x_col]]), neighbors$x)) ||
        !isTRUE(all.equal(as.numeric(coords[[y_col]]), neighbors$y))) {
      stop("neighbors coordinates and pixel_matrix row order do not match.", call. = FALSE)
    }
  }
  feature_seeds <- rep(NA_integer_, length(fcols))
  if (!is.null(seed)) {
    set.seed(seed)
    feature_seeds <- sample.int(.Machine$integer.max, length(fcols))
  }
  rows <- lapply(seq_along(fcols), function(feature_index) {
    feature <- fcols[feature_index]
    calculation <- tryCatch(
      compute_binspect_feature(
        values = pixel_matrix[[feature]], neighbors = neighbors,
        bin_method = bin_method, percentage_rank = percentage_rank,
        inference = inference, n_perm = n_perm,
        seed = if (is.na(feature_seeds[feature_index])) NULL else feature_seeds[feature_index]
      ),
      error = function(error) error
    )
    if (inherits(calculation, "error")) {
      return(data.frame(
        feature = feature, p_value = NA_real_, odds_ratio = NA_real_,
        high_high_edge_weight = NA_real_, n_high = NA_integer_,
        high_fraction = NA_real_, n_effective = NA_integer_, n_hub = NA_integer_,
        mean_expression = NA_real_, mean_high_expression = NA_real_,
        bin_method = bin_method, inference = inference,
        status = "not_tested", error_message = conditionMessage(calculation),
        stringsAsFactors = FALSE
      ))
    }
    data.frame(
      feature = feature, p_value = calculation$p_value,
      odds_ratio = calculation$odds_ratio,
      high_high_edge_weight = calculation$high_high_edge_weight,
      n_high = calculation$n_high, high_fraction = calculation$high_fraction,
      n_effective = calculation$n_effective, n_hub = calculation$n_hub,
      mean_expression = calculation$mean_expression,
      mean_high_expression = calculation$mean_high_expression,
      bin_method = calculation$bin_method, inference = calculation$inference,
      status = "tested", error_message = NA_character_, stringsAsFactors = FALSE
    )
  })
  result <- do.call(rbind, rows)
  result$adj_p_value <- stats::p.adjust(result$p_value, method = p_adjust_method)
  rownames(result) <- NULL
  result[order(result$adj_p_value, -result$odds_ratio, na.last = TRUE), , drop = FALSE]
}

compare_svg_methods <- function(moran_result,
                                geary_result,
                                binspect_result,
                                alpha = 0.05) {
  required_columns(moran_result, c("feature", "morans_i", "adj_p_value"), "Moran result")
  required_columns(geary_result, c("feature", "gearys_c", "adj_p_value"), "Geary result")
  required_columns(binspect_result, c("feature", "odds_ratio", "adj_p_value"), "binSpect result")
  if (anyDuplicated(moran_result$feature) || anyDuplicated(geary_result$feature) ||
      anyDuplicated(binspect_result$feature)) {
    stop("Each result must contain at most one row per feature.", call. = FALSE)
  }
  if (length(alpha) != 1L || !is.finite(alpha) || alpha <= 0 || alpha >= 1) {
    stop("alpha must be one finite value strictly between 0 and 1.", call. = FALSE)
  }
  merged <- Reduce(function(left, right) merge(left, right, by = "feature", sort = FALSE), list(
    data.frame(feature = moran_result$feature, morans_i = moran_result$morans_i,
      moran_adj_p = moran_result$adj_p_value),
    data.frame(feature = geary_result$feature, gearys_c = geary_result$gearys_c,
      geary_adj_p = geary_result$adj_p_value),
    data.frame(feature = binspect_result$feature, odds_ratio = binspect_result$odds_ratio,
      binspect_adj_p = binspect_result$adj_p_value)
  ))
  p_columns <- c("moran_adj_p", "geary_adj_p", "binspect_adj_p")
  tested <- is.finite(as.matrix(merged[p_columns]))
  significant <- tested & as.matrix(merged[p_columns]) < alpha
  merged$moran_significant <- ifelse(tested[, 1], significant[, 1], NA)
  merged$geary_significant <- ifelse(tested[, 2], significant[, 2], NA)
  merged$binspect_significant <- ifelse(tested[, 3], significant[, 3], NA)
  merged$n_methods_tested <- rowSums(tested)
  merged$n_methods_significant <- rowSums(significant)
  merged$all_tested_methods_significant <- merged$n_methods_tested > 0 &
    merged$n_methods_significant == merged$n_methods_tested
  merged$alpha <- alpha
  output <- merged[order(
    -merged$n_methods_significant, -merged$n_methods_tested,
    merged$moran_adj_p, na.last = TRUE
  ), , drop = FALSE]
  rownames(output) <- NULL
  attr(output, "interpretation") <- paste(
    "The three methods share measurements and a spatial graph. Concordance is",
    "a method-sensitivity result, not independent validation or proof of identity."
  )
  output
}

find_msi_lcms_candidates <- function(msi_features,
                                     lcms_features,
                                     msi_mz_col = "mz",
                                     lcms_mz_col = "mz",
                                     msi_mode_col = "ion_mode",
                                     lcms_mode_col = "ion_mode",
                                     ppm = 5) {
  if (!is.data.frame(msi_features) || !is.data.frame(lcms_features)) {
    stop("Both MSI features and LC-MS features must be data frames.", call. = FALSE)
  }
  if (!is.numeric(ppm) || length(ppm) != 1L || !is.finite(ppm) || ppm < 0) {
    stop("ppm must be one finite non-negative number.", call. = FALSE)
  }
  msi_mz_col <- first_existing_column(msi_features, msi_mz_col, "MSI m/z column")
  lcms_mz_col <- first_existing_column(lcms_features, lcms_mz_col, "LC-MS m/z column")
  msi_mz <- suppressWarnings(as.numeric(msi_features[[msi_mz_col]]))
  lcms_mz <- suppressWarnings(as.numeric(lcms_features[[lcms_mz_col]]))
  if (any(!is.finite(msi_mz) | msi_mz <= 0)) {
    stop("MSI m/z values must be numeric, finite, and positive.", call. = FALSE)
  }
  if (any(!is.finite(lcms_mz) | lcms_mz <= 0)) {
    stop("LC-MS m/z values must be numeric, finite, and positive.", call. = FALSE)
  }

  choose_optional <- function(data, candidates) {
    columns <- candidates[candidates %in% names(data)]
    if (length(columns)) columns[1L] else NULL
  }
  msi_mode_col <- choose_optional(msi_features, msi_mode_col)
  lcms_mode_col <- choose_optional(lcms_features, lcms_mode_col)
  mode_match <- !is.null(msi_mode_col) && !is.null(lcms_mode_col)
  tolerance <- ppm + .Machine$double.eps * 1e6
  rows <- lapply(seq_along(msi_mz), function(i) {
    candidate_index <- seq_along(lcms_mz)
    if (mode_match) {
      msi_mode <- as.character(msi_features[[msi_mode_col]][i])
      lcms_mode <- as.character(lcms_features[[lcms_mode_col]])
      candidate_index <- candidate_index[!is.na(msi_mode) & !is.na(lcms_mode) & lcms_mode == msi_mode]
    }
    if (!length(candidate_index)) return(NULL)
    error <- abs(lcms_mz[candidate_index] - msi_mz[i]) / msi_mz[i] * 1e6
    keep <- is.finite(error) & error <= tolerance
    if (!any(keep)) return(NULL)
    candidate_index <- candidate_index[keep]
    error <- error[keep]
    order_index <- order(error, candidate_index)
    data.frame(
      msi_row = i,
      lcms_row = candidate_index[order_index],
      msi_feature_id = if ("feature_id" %in% names(msi_features)) as.character(msi_features$feature_id[i]) else as.character(i),
      lcms_feature_id = as.character(candidate_index[order_index]),
      msi_mz = msi_mz[i],
      lcms_mz = lcms_mz[candidate_index[order_index]],
      ion_mode = if (mode_match) as.character(msi_features[[msi_mode_col]][i]) else NA_character_,
      ppm_error = error[order_index],
      candidate_rank = seq_along(order_index),
      stringsAsFactors = FALSE
    )
  })
  rows <- Filter(Negate(is.null), rows)
  if (!length(rows)) return(data.frame(
    msi_row = integer(), lcms_row = integer(), msi_feature_id = character(),
    lcms_feature_id = character(), msi_mz = numeric(), lcms_mz = numeric(),
    ion_mode = character(), ppm_error = numeric(), candidate_rank = integer(),
    candidate_count = integer(), stringsAsFactors = FALSE
  ))
  result <- do.call(rbind, rows)
  counts <- table(result$msi_row)
  result$candidate_count <- as.integer(counts[as.character(result$msi_row)])
  rownames(result) <- NULL
  result
}

assign_msi_lcms_candidates <- function(candidates,
                                       n_msi,
                                       n_lcms,
                                       ppm,
                                       method = c("greedy", "optimal")) {
  method <- match.arg(method)
  if (!is.data.frame(candidates) ||
      !all(c("msi_row", "lcms_row", "ppm_error") %in% names(candidates))) {
    stop("candidates must come from find_msi_lcms_candidates().", call. = FALSE)
  }
  if (!nrow(candidates)) return(integer())
  if (method == "greedy") {
    used_lcms <- rep(FALSE, n_lcms)
    selected <- integer()
    for (i in seq_len(n_msi)) {
      available <- which(candidates$msi_row == i & !used_lcms[candidates$lcms_row])
      if (!length(available)) next
      chosen <- available[which.min(candidates$ppm_error[available])]
      selected <- c(selected, chosen)
      used_lcms[candidates$lcms_row[chosen]] <- TRUE
    }
    return(selected)
  }
  # Every MSI row receives either a unique real LC-MS column or one of n_msi
  # rejection columns.  The rejection penalty is larger than the maximum
  # possible change in total feasible ppm cost, so the solution first
  # maximizes match cardinality and then minimizes total ppm error.
  rejection_penalty <- (n_msi + 1) * (ppm + 1)
  infeasible_penalty <- (n_msi + 1) * rejection_penalty
  cost <- matrix(infeasible_penalty, nrow = n_msi, ncol = n_lcms + n_msi)
  cost[, n_lcms + seq_len(n_msi)] <- rejection_penalty
  cost[cbind(candidates$msi_row, candidates$lcms_row)] <- candidates$ppm_error
  # Rectangular Hungarian algorithm. Keeping this small implementation in base
  # R avoids making scientifically important matching depend on an optional
  # package. The matrix always has at least as many columns as rows.
  solve_rectangular_lsap <- function(x) {
    nr <- nrow(x); nc <- ncol(x)
    if (nr > nc || any(!is.finite(x)) || any(x < 0)) {
      stop("Internal assignment costs must be finite, non-negative, and rectangular.", call. = FALSE)
    }
    u <- numeric(nr)
    v <- numeric(nc + 1L)
    p <- integer(nc + 1L)
    way <- integer(nc + 1L)
    for (i in seq_len(nr)) {
      p[1L] <- i
      j0 <- 1L
      minv <- rep(Inf, nc)
      used <- rep(FALSE, nc + 1L)
      repeat {
        used[j0] <- TRUE
        i0 <- p[j0]
        columns <- which(!used[-1L])
        current <- x[i0, columns] - u[i0] - v[columns + 1L]
        improve <- current < minv[columns]
        if (any(improve)) {
          improved_columns <- columns[improve]
          minv[improved_columns] <- current[improve]
          way[improved_columns + 1L] <- j0
        }
        delta <- min(minv[columns])
        j1_column <- columns[which.min(minv[columns])]
        used_indices <- which(used)
        real_rows <- p[used_indices]
        real <- real_rows > 0L
        u[real_rows[real]] <- u[real_rows[real]] + delta
        v[used_indices] <- v[used_indices] - delta
        minv[columns] <- minv[columns] - delta
        j0 <- j1_column + 1L
        if (p[j0] == 0L) break
      }
      repeat {
        j1 <- way[j0]
        p[j0] <- p[j1]
        j0 <- j1
        if (j0 == 1L) break
      }
    }
    assignment <- integer(nr)
    for (j in 2:(nc + 1L)) if (p[j] > 0L) assignment[p[j]] <- j - 1L
    assignment
  }
  assignment <- solve_rectangular_lsap(cost)
  real_rows <- which(assignment <= n_lcms)
  if (!length(real_rows)) return(integer())
  keys <- paste(candidates$msi_row, candidates$lcms_row, sep = ":")
  selected <- match(paste(real_rows, assignment[real_rows], sep = ":"), keys)
  selected[!is.na(selected)]
}

cross_validate_msi_lcms <- function(msi_features,
                                    lcms_features,
                                    msi_mz_col = "mz",
                                    lcms_mz_col = "mz",
                                    msi_mode_col = "ion_mode",
                                    lcms_mode_col = "ion_mode",
                                    msi_log2fc_col = "log2fc",
                                    lcms_log2fc_col = "log2fc",
                                    lcms_id_col = c("id", "compound_id", "feature_id", "name"),
                                    lcms_name_col = c("name", "compound", "compound_name"),
                                    lcms_rt_col = c("rt", "retention_time"),
                                    ppm = 5,
                                    assignment_method = c("greedy", "optimal")) {
  assignment_method <- match.arg(assignment_method)
  if (!is.data.frame(msi_features) || !is.data.frame(lcms_features)) {
    stop("Both MSI features and LC-MS features must be data frames.", call. = FALSE)
  }

  msi_mz_col <- first_existing_column(msi_features, msi_mz_col, "MSI m/z column")
  lcms_mz_col <- first_existing_column(lcms_features, lcms_mz_col, "LC-MS m/z column")
  required_columns(msi_features, c(msi_mz_col), "MSI feature table")
  required_columns(lcms_features, c(lcms_mz_col), "LC-MS feature table")

  choose_optional <- function(data, candidates) {
    cols <- candidates[candidates %in% names(data)]
    if (length(cols) > 0) cols[1] else NULL
  }

  msi_mode_col <- choose_optional(msi_features, msi_mode_col)
  lcms_mode_col <- choose_optional(lcms_features, lcms_mode_col)
  msi_log2fc_col <- choose_optional(msi_features, msi_log2fc_col)
  lcms_log2fc_col <- choose_optional(lcms_features, lcms_log2fc_col)
  lcms_id_col <- choose_optional(lcms_features, lcms_id_col)
  lcms_name_col <- choose_optional(lcms_features, lcms_name_col)
  lcms_rt_col <- choose_optional(lcms_features, lcms_rt_col)

  msi_mz <- suppressWarnings(as.numeric(msi_features[[msi_mz_col]]))
  lcms_mz <- suppressWarnings(as.numeric(lcms_features[[lcms_mz_col]]))
  if (any(!is.finite(msi_mz))) {
    stop("MSI m/z values must be numeric and finite.", call. = FALSE)
  }
  if (any(!is.finite(lcms_mz))) {
    stop("LC-MS m/z values must be numeric and finite.", call. = FALSE)
  }

  mode_match <- !is.null(msi_mode_col) && !is.null(lcms_mode_col)
  candidates <- find_msi_lcms_candidates(
    msi_features, lcms_features, msi_mz_col, lcms_mz_col,
    msi_mode_col, lcms_mode_col, ppm
  )
  selected <- assign_msi_lcms_candidates(
    candidates, nrow(msi_features), nrow(lcms_features), ppm,
    method = assignment_method
  )
  rows <- lapply(selected, function(candidate_row) {
    candidate <- candidates[candidate_row, , drop = FALSE]
    i <- candidate$msi_row
    chosen <- candidate$lcms_row

    msi_log2fc <- if (!is.null(msi_log2fc_col)) suppressWarnings(as.numeric(msi_features[[msi_log2fc_col]][i])) else NA_real_
    lcms_log2fc <- if (!is.null(lcms_log2fc_col)) suppressWarnings(as.numeric(lcms_features[[lcms_log2fc_col]][chosen])) else NA_real_
    direction_agreement <- NA
    if (is.finite(msi_log2fc) && is.finite(lcms_log2fc)) {
      direction_agreement <- sign(msi_log2fc) == sign(lcms_log2fc)
    }

    lcms_id <- if (!is.null(lcms_id_col)) as.character(lcms_features[[lcms_id_col]][chosen]) else NA_character_
    lcms_name <- if (!is.null(lcms_name_col)) as.character(lcms_features[[lcms_name_col]][chosen]) else NA_character_
    match_type <- "feature_level_orthogonal_support"
    lcms_rt <- if (!is.null(lcms_rt_col)) lcms_features[[lcms_rt_col]][chosen] else NA_real_

    data.frame(
      msi_feature_id = if ("feature_id" %in% names(msi_features)) as.character(msi_features$feature_id[i]) else NA_character_,
      msi_mz = msi_mz[i],
      lcms_feature_id = if (!is.null(lcms_id_col)) as.character(lcms_features[[lcms_id_col]][chosen]) else NA_character_,
      lcms_mz = lcms_mz[chosen],
      ion_mode = if (mode_match) as.character(msi_features[[msi_mode_col]][i]) else NA_character_,
      ppm_error = candidate$ppm_error,
      msi_log2fc = msi_log2fc,
      lcms_log2fc = lcms_log2fc,
      direction_agreement = direction_agreement,
      match_type = match_type,
      lcms_rt = if (!is.null(lcms_rt_col)) lcms_rt else NA_real_,
      lcms_name = lcms_name,
      candidate_count = candidate$candidate_count,
      assignment_method = assignment_method,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  })

  result <- if (length(rows) > 0) do.call(rbind, rows) else structure(data.frame(
    msi_feature_id = character(),
    msi_mz = numeric(),
    lcms_feature_id = character(),
    lcms_mz = numeric(),
    ion_mode = character(),
    ppm_error = numeric(),
    msi_log2fc = numeric(),
    lcms_log2fc = numeric(),
    direction_agreement = logical(),
    match_type = character(),
    lcms_rt = numeric(),
    lcms_name = character(),
    candidate_count = integer(),
    assignment_method = character(),
    stringsAsFactors = FALSE,
    check.names = FALSE
  ), class = "data.frame")

  summary <- data.frame(
    total_msi_features = nrow(msi_features),
    matched_features = nrow(result),
    ambiguous_msi_features = sum(result$candidate_count > 1L, na.rm = TRUE),
    assignment_method = assignment_method,
    feature_level_supported_matches = sum(result$match_type == "feature_level_orthogonal_support", na.rm = TRUE),
    agreement_rate = if (sum(!is.na(result$direction_agreement)) == 0) NA_real_ else mean(result$direction_agreement, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
  attr(result, "summary") <- summary
  result
}

#' Validate explicitly defined diagnostic-ion evidence in MS/MS scans
#'
#' This function evaluates user-supplied diagnostic ions without assigning or
#' inferring chemical identities. A core scan must satisfy both the precursor
#' tolerance and an isolation window that covers the target precursor.
#'
#' @param precursor_scans Data frame containing `scan_number`, `precursor_mz`,
#'   `isolation_target_mz`, `isolation_lower_offset`,
#'   `isolation_upper_offset`, and `base_peak_intensity`.
#' @param fragment_peaks Data frame containing `scan_number`, `fragment_mz`,
#'   and `fragment_intensity`. Optional provenance columns are retained.
#' @param diagnostic_ions Data frame with unique `label`, numeric `target_mz`,
#'   and `evidence_role` columns.
#' @param precursor_target Numeric target precursor m/z.
#' @param precursor_ppm_tolerance Non-negative precursor tolerance in ppm.
#' @param fragment_tolerance Non-negative fragment tolerance.
#' @param fragment_tolerance_unit Either `"Da"` or `"ppm"`.
#' @param min_relative_intensity Minimum fragment/base-peak intensity fraction.
#' @param min_core_scans Minimum number of qualifying core scans for repeated
#'   support.
#' @param required_diagnostic_labels Labels that must all achieve repeated
#'   diagnostic support for the overall rule to pass.
#' @return A list containing scan-level evidence, diagnostic summaries, and an
#'   explicit overall support result.
validate_msms_fragment_evidence <- function(
    precursor_scans,
    fragment_peaks,
    diagnostic_ions,
    precursor_target,
    precursor_ppm_tolerance = 10,
    fragment_tolerance = 0.02,
    fragment_tolerance_unit = c("Da", "ppm"),
    min_relative_intensity = 0.01,
    min_core_scans = 2L,
    required_diagnostic_labels) {
  fragment_tolerance_unit <- match.arg(fragment_tolerance_unit)
  if (!is.data.frame(precursor_scans) || !is.data.frame(fragment_peaks) ||
      !is.data.frame(diagnostic_ions)) {
    stop("Precursor scans, fragment peaks, and diagnostic ions must be data frames.", call. = FALSE)
  }
  scan_required <- c(
    "scan_number", "precursor_mz", "isolation_target_mz",
    "isolation_lower_offset", "isolation_upper_offset", "base_peak_intensity"
  )
  peak_required <- c("scan_number", "fragment_mz", "fragment_intensity")
  diagnostic_required <- c("label", "target_mz", "evidence_role")
  required_columns(precursor_scans, scan_required, "precursor scan metadata")
  required_columns(fragment_peaks, peak_required, "fragment peak table")
  required_columns(diagnostic_ions, diagnostic_required, "diagnostic ion definitions")

  scalar_nonnegative <- function(value, name, positive = FALSE) {
    value <- suppressWarnings(as.numeric(value))
    valid <- length(value) == 1L && is.finite(value) &&
      if (positive) value > 0 else value >= 0
    if (!valid) stop(name, if (positive) " must be a positive finite number." else " must be a non-negative finite number.", call. = FALSE)
    value
  }
  precursor_target <- scalar_nonnegative(precursor_target, "precursor_target", positive = TRUE)
  precursor_ppm_tolerance <- scalar_nonnegative(precursor_ppm_tolerance, "precursor_ppm_tolerance")
  fragment_tolerance <- scalar_nonnegative(fragment_tolerance, "fragment_tolerance")
  min_relative_intensity <- scalar_nonnegative(min_relative_intensity, "min_relative_intensity")
  min_core_scans <- suppressWarnings(as.integer(min_core_scans))
  if (length(min_core_scans) != 1L || is.na(min_core_scans) || min_core_scans < 2L) {
    stop("min_core_scans must be a single integer of at least 2.", call. = FALSE)
  }
  required_diagnostic_labels <- as.character(required_diagnostic_labels)
  if (length(required_diagnostic_labels) == 0L || anyNA(required_diagnostic_labels) ||
      any(!nzchar(required_diagnostic_labels))) {
    stop("required_diagnostic_labels must explicitly name at least one diagnostic label.", call. = FALSE)
  }

  diagnostic_ions$label <- as.character(diagnostic_ions$label)
  diagnostic_ions$evidence_role <- as.character(diagnostic_ions$evidence_role)
  diagnostic_ions$target_mz <- suppressWarnings(as.numeric(diagnostic_ions$target_mz))
  if (anyNA(diagnostic_ions$label) || any(!nzchar(diagnostic_ions$label)) ||
      anyDuplicated(diagnostic_ions$label)) {
    stop("Diagnostic ion labels must be non-missing, non-empty, and unique.", call. = FALSE)
  }
  if (anyNA(diagnostic_ions$evidence_role) || any(!nzchar(diagnostic_ions$evidence_role))) {
    stop("Each diagnostic ion must have an explicit evidence_role.", call. = FALSE)
  }
  if (any(!is.finite(diagnostic_ions$target_mz) | diagnostic_ions$target_mz <= 0)) {
    stop("Diagnostic target m/z values must be positive and finite.", call. = FALSE)
  }
  missing_required <- setdiff(required_diagnostic_labels, diagnostic_ions$label)
  if (length(missing_required) > 0L) {
    stop("Required diagnostic labels are absent from diagnostic_ions: ",
         paste(missing_required, collapse = ", "), call. = FALSE)
  }

  numeric_scan_columns <- setdiff(scan_required, "scan_number")
  for (column in numeric_scan_columns) {
    precursor_scans[[column]] <- suppressWarnings(as.numeric(precursor_scans[[column]]))
  }
  if (anyDuplicated(as.character(precursor_scans$scan_number))) {
    stop("precursor_scans must contain one row per scan_number.", call. = FALSE)
  }
  if (any(!is.finite(as.matrix(precursor_scans[numeric_scan_columns])))) {
    stop("Required precursor scan metadata must be numeric and finite.", call. = FALSE)
  }
  if (any(precursor_scans$isolation_lower_offset < 0 |
          precursor_scans$isolation_upper_offset < 0)) {
    stop("Isolation window offsets must be non-negative.", call. = FALSE)
  }
  if (any(precursor_scans$base_peak_intensity <= 0)) {
    stop("base_peak_intensity must be positive.", call. = FALSE)
  }

  fragment_peaks$fragment_mz <- suppressWarnings(as.numeric(fragment_peaks$fragment_mz))
  fragment_peaks$fragment_intensity <- suppressWarnings(as.numeric(fragment_peaks$fragment_intensity))
  if (any(!is.finite(fragment_peaks$fragment_mz) | fragment_peaks$fragment_mz <= 0) ||
      any(!is.finite(fragment_peaks$fragment_intensity) | fragment_peaks$fragment_intensity < 0)) {
    stop("Fragment m/z and intensity values must be finite; m/z must be positive and intensity non-negative.", call. = FALSE)
  }
  unknown_scans <- setdiff(as.character(fragment_peaks$scan_number), as.character(precursor_scans$scan_number))
  if (length(unknown_scans) > 0L) {
    stop("Fragment peaks reference scan_number values absent from precursor_scans.", call. = FALSE)
  }

  provenance_columns <- c("peak_source", "extraction_method", "extraction_parameters", "software_version")
  for (column in provenance_columns) {
    if (!column %in% names(fragment_peaks)) fragment_peaks[[column]] <- NA_character_
    fragment_peaks[[column]] <- as.character(fragment_peaks[[column]])
  }
  profile_rows <- !is.na(fragment_peaks$peak_source) &
    grepl("profile", fragment_peaks$peak_source, ignore.case = TRUE)
  if (any(profile_rows)) {
    incomplete <- vapply(provenance_columns[-1L], function(column) {
      any(is.na(fragment_peaks[[column]][profile_rows]) | !nzchar(fragment_peaks[[column]][profile_rows]))
    }, logical(1))
    if (any(incomplete)) {
      stop("Profile-derived peaks require extraction_method, extraction_parameters, and software_version.", call. = FALSE)
    }
  }

  precursor_scans$precursor_ppm_error <-
    (precursor_scans$precursor_mz - precursor_target) / precursor_target * 1e6
  precursor_scans$selected_precursor_within_tolerance <-
    abs(precursor_scans$precursor_ppm_error) <= precursor_ppm_tolerance + .Machine$double.eps * 1e6
  precursor_scans$isolation_window_covers_target <-
    precursor_target >= precursor_scans$isolation_target_mz - precursor_scans$isolation_lower_offset &
    precursor_target <= precursor_scans$isolation_target_mz + precursor_scans$isolation_upper_offset
  precursor_scans$core_scan <- precursor_scans$selected_precursor_within_tolerance &
    precursor_scans$isolation_window_covers_target

  scan_order <- order(as.character(precursor_scans$scan_number))
  diagnostic_order <- order(diagnostic_ions$label)
  precursor_scans <- precursor_scans[scan_order, , drop = FALSE]
  diagnostic_ions <- diagnostic_ions[diagnostic_order, , drop = FALSE]
  fragment_peaks <- fragment_peaks[order(as.character(fragment_peaks$scan_number), fragment_peaks$fragment_mz), , drop = FALSE]

  evidence_rows <- vector("list", nrow(precursor_scans) * nrow(diagnostic_ions))
  row_index <- 1L
  for (scan_index in seq_len(nrow(precursor_scans))) {
    scan <- precursor_scans[scan_index, , drop = FALSE]
    peaks <- fragment_peaks[as.character(fragment_peaks$scan_number) == as.character(scan$scan_number), , drop = FALSE]
    for (diagnostic_index in seq_len(nrow(diagnostic_ions))) {
      diagnostic <- diagnostic_ions[diagnostic_index, , drop = FALSE]
      if (fragment_tolerance_unit == "Da") {
        errors <- peaks$fragment_mz - diagnostic$target_mz
        within <- abs(errors) <= fragment_tolerance + .Machine$double.eps
      } else {
        errors <- (peaks$fragment_mz - diagnostic$target_mz) / diagnostic$target_mz * 1e6
        within <- abs(errors) <= fragment_tolerance + .Machine$double.eps * 1e6
      }
      candidates <- peaks[within, , drop = FALSE]
      if (nrow(candidates) > 0L) {
        chosen_index <- which.max(candidates$fragment_intensity)
        chosen <- candidates[chosen_index, , drop = FALSE]
        observed_mz <- chosen$fragment_mz
        error_da <- observed_mz - diagnostic$target_mz
        error_ppm <- error_da / diagnostic$target_mz * 1e6
        raw_intensity <- chosen$fragment_intensity
        relative_intensity <- raw_intensity / scan$base_peak_intensity
        provenance <- chosen[1L, provenance_columns, drop = FALSE]
        mass_match <- TRUE
      } else {
        observed_mz <- error_da <- error_ppm <- raw_intensity <- relative_intensity <- NA_real_
        provenance <- as.data.frame(stats::setNames(rep(list(NA_character_), length(provenance_columns)), provenance_columns), stringsAsFactors = FALSE)
        mass_match <- FALSE
      }
      qualifies <- mass_match && isTRUE(scan$core_scan) &&
        is.finite(relative_intensity) && relative_intensity >= min_relative_intensity
      evidence_rows[[row_index]] <- data.frame(
        scan_number = scan$scan_number,
        diagnostic_label = diagnostic$label,
        evidence_role = diagnostic$evidence_role,
        diagnostic_target_mz = diagnostic$target_mz,
        observed_mz = observed_mz,
        mass_error_da = error_da,
        mass_error_ppm = error_ppm,
        fragment_tolerance = fragment_tolerance,
        fragment_tolerance_unit = fragment_tolerance_unit,
        raw_intensity = raw_intensity,
        base_peak_intensity = scan$base_peak_intensity,
        relative_base_peak_intensity = relative_intensity,
        selected_precursor_mz = scan$precursor_mz,
        precursor_ppm_error = scan$precursor_ppm_error,
        selected_precursor_within_tolerance = scan$selected_precursor_within_tolerance,
        isolation_window_covers_target = scan$isolation_window_covers_target,
        core_scan = scan$core_scan,
        mass_match = mass_match,
        meets_relative_intensity = mass_match && is.finite(relative_intensity) && relative_intensity >= min_relative_intensity,
        qualifies_core_strength_rule = qualifies,
        peak_source = provenance$peak_source,
        extraction_method = provenance$extraction_method,
        extraction_parameters = provenance$extraction_parameters,
        software_version = provenance$software_version,
        stringsAsFactors = FALSE
      )
      row_index <- row_index + 1L
    }
  }
  scan_evidence <- do.call(rbind, evidence_rows)
  rownames(scan_evidence) <- NULL

  summary_rows <- lapply(seq_len(nrow(diagnostic_ions)), function(index) {
    diagnostic <- diagnostic_ions[index, , drop = FALSE]
    rows <- scan_evidence[scan_evidence$diagnostic_label == diagnostic$label, , drop = FALSE]
    n_qualifying <- sum(rows$qualifies_core_strength_rule)
    n_mass_matches <- sum(rows$mass_match)
    grade <- if (n_qualifying >= min_core_scans) {
      "repeated_diagnostic_support"
    } else if (n_qualifying == 1L) {
      "single_scan_support"
    } else if (n_mass_matches > 0L) {
      "trace_match"
    } else {
      "not_detected"
    }
    data.frame(
      diagnostic_label = diagnostic$label,
      evidence_role = diagnostic$evidence_role,
      diagnostic_target_mz = diagnostic$target_mz,
      n_mass_matched_scans = n_mass_matches,
      n_qualifying_core_scans = n_qualifying,
      evidence_grade = grade,
      stringsAsFactors = FALSE
    )
  })
  diagnostic_summary <- do.call(rbind, summary_rows)
  rownames(diagnostic_summary) <- NULL
  required_grades <- diagnostic_summary$evidence_grade[
    match(required_diagnostic_labels, diagnostic_summary$diagnostic_label)
  ]
  overall <- data.frame(
    support_label = "fragment_level_orthogonal_support",
    support = all(required_grades == "repeated_diagnostic_support"),
    required_diagnostic_labels = paste(required_diagnostic_labels, collapse = ";"),
    minimum_core_scan_repeats = min_core_scans,
    stringsAsFactors = FALSE
  )
  result <- list(
    scan_diagnostic_evidence = scan_evidence,
    diagnostic_summary = diagnostic_summary,
    overall = overall
  )
  result
}

backmap_sample_scores <- function(score_data, sample_mapping, score_column) {
  required_columns(score_data, c("Sample", score_column), "Score data")
  required_columns(sample_mapping, c("sample_id", "pixel_ids"), "Sample mapping")

  merged <- merge(
    sample_mapping,
    score_data[, c("Sample", score_column), drop = FALSE],
    by.x = "sample_id",
    by.y = "Sample",
    all.x = FALSE
  )
  if (nrow(merged) == 0) {
    stop(
      "No samples matched between score data and sample_mapping. Check sample name formatting.",
      call. = FALSE
    )
  }

  rows <- lapply(seq_len(nrow(merged)), function(i) {
    ids <- strsplit(merged$pixel_ids[i], ",", fixed = TRUE)[[1]]
    data.frame(
      pixel_id = ids,
      sample_id = merged$sample_id[i],
      score = merged[[score_column]][i],
      stringsAsFactors = FALSE
    )
  })

  do.call(rbind, rows)
}

resolve_feature_column <- function(pixel_matrix, feature) {
  fcols <- feature_columns(pixel_matrix)
  if (length(fcols) == 0) {
    stop("No mz feature columns are available.", call. = FALSE)
  }

  if (length(feature) != 1) {
    stop("feature must be a single mz value or feature column name.", call. = FALSE)
  }

  if (is.character(feature) && feature %in% fcols) {
    return(feature)
  }

  target_mz <- suppressWarnings(as.numeric(sub("^mz_", "", as.character(feature))))
  feature_mz <- suppressWarnings(as.numeric(sub("^mz_", "", fcols)))
  if (!is.finite(target_mz) || all(!is.finite(feature_mz))) {
    stop("Cannot match feature to an mz_* column.", call. = FALSE)
  }

  fcols[which.min(abs(feature_mz - target_mz))]
}

plot_ion_image <- function(pixel_matrix,
                           feature,
                           transform = c("identity", "log10"),
                           title = NULL) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for plotting.", call. = FALSE)
  }
  if (!requireNamespace("viridis", quietly = TRUE)) {
    stop("Package 'viridis' is required for plotting.", call. = FALSE)
  }

  transform <- match.arg(transform)
  required_columns(pixel_matrix, c("x", "y"), "Pixel matrix")
  column_name <- resolve_feature_column(pixel_matrix, feature)
  values <- pixel_matrix[[column_name]]
  if (transform == "log10") values <- log10(values + 1)

  plot_data <- data.frame(
    x = pixel_matrix$x,
    y = pixel_matrix$y,
    intensity = values
  )
  x_steps <- diff(sort(unique(plot_data$x)))
  y_steps <- diff(sort(unique(plot_data$y)))
  tile_width <- if (any(x_steps > 0)) min(x_steps[x_steps > 0]) else 1
  tile_height <- if (any(y_steps > 0)) min(y_steps[y_steps > 0]) else 1

  if (is.null(title)) title <- column_name

  ggplot2::ggplot(plot_data, ggplot2::aes(x = .data[["x"]], y = .data[["y"]], fill = .data[["intensity"]])) +
    ggplot2::geom_tile(width = tile_width, height = tile_height) +
    ggplot2::coord_fixed() +
    viridis::scale_fill_viridis(option = "viridis") +
    ggplot2::labs(title = title, x = "x", y = "y", fill = "Intensity") +
    ggplot2::theme_minimal()
}

plot_feature_thumbnails <- function(pixel_matrix,
                                    page = 1,
                                    per_page = 16,
                                    features = NULL,
                                    transform = c("identity", "log10")) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for plotting.", call. = FALSE)
  }
  if (!requireNamespace("viridis", quietly = TRUE)) {
    stop("Package 'viridis' is required for plotting.", call. = FALSE)
  }

  transform <- match.arg(transform)
  required_columns(pixel_matrix, c("x", "y"), "Pixel matrix")
  fcols <- if (is.null(features)) feature_columns(pixel_matrix) else vapply(
    features,
    function(feature) resolve_feature_column(pixel_matrix, feature),
    character(1)
  )
  fcols <- unique(fcols)
  if (length(fcols) == 0) {
    stop("No mz feature columns are available for thumbnails.", call. = FALSE)
  }

  page <- max(1, as.integer(page))
  per_page <- max(1, as.integer(per_page))
  start <- (page - 1) * per_page + 1
  end <- min(start + per_page - 1, length(fcols))
  if (start > length(fcols)) {
    stop("page is out of range. Maximum page is ", ceiling(length(fcols) / per_page), ".", call. = FALSE)
  }

  page_cols <- fcols[start:end]
  rows <- lapply(page_cols, function(column_name) {
    values <- pixel_matrix[[column_name]]
    if (transform == "log10") values <- log10(values + 1)
    data.frame(
      x = pixel_matrix$x,
      y = pixel_matrix$y,
      feature = column_name,
      intensity = values,
      stringsAsFactors = FALSE
    )
  })
  plot_data <- do.call(rbind, rows)

  ggplot2::ggplot(plot_data, ggplot2::aes(x = .data[["x"]], y = .data[["y"]], fill = .data[["intensity"]])) +
    ggplot2::geom_raster() +
    ggplot2::coord_fixed() +
    ggplot2::facet_wrap(~feature, ncol = 4) +
    viridis::scale_fill_viridis(option = "viridis") +
    ggplot2::labs(
      title = paste0("Feature thumbnails page ", page, " / ", ceiling(length(fcols) / per_page)),
      x = "x",
      y = "y",
      fill = "Intensity"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(strip.text = ggplot2::element_text(size = 8))
}

extract_metaboanalyst_pathway_hits <- function(metabo_object, pathway = NULL) {
  if (!is.list(metabo_object)) {
    stop("metabo_object must be a MetaboAnalystR mSet-like list.", call. = FALSE)
  }

  rows <- list()

  if (!is.null(metabo_object$analSet$ora.hits)) {
    ora_hits <- metabo_object$analSet$ora.hits
    pathway_ids <- names(ora_hits)
    pathway_names <- pathway_ids
    if (!is.null(metabo_object$analSet$pathwayMemberTable)) {
      member_table <- metabo_object$analSet$pathwayMemberTable
      if (all(c("Pathway", "Hit IDs") %in% names(member_table))) {
        hit_id_lookup <- stats::setNames(member_table$Pathway, member_table[["Hit IDs"]])
        pathway_names <- vapply(seq_along(ora_hits), function(i) {
          ids <- paste(unname(ora_hits[[i]]), collapse = "; ")
          if (ids %in% names(hit_id_lookup)) hit_id_lookup[[ids]] else pathway_ids[[i]]
        }, character(1))
      }
    }

    rows <- c(rows, lapply(seq_along(ora_hits), function(i) {
      hits <- ora_hits[[i]]
      if (length(hits) == 0) return(NULL)
      data.frame(
        source = "pathora",
        pathway_id = pathway_ids[[i]],
        pathway = pathway_names[[i]],
        compound_id = unname(as.character(hits)),
        compound_name = names(hits),
        mz = NA_real_,
        stringsAsFactors = FALSE
      )
    }))
  }

  if (!is.null(metabo_object$path.hits) && !is.null(metabo_object$path.nms)) {
    rows <- c(rows, lapply(seq_along(metabo_object$path.hits), function(i) {
      hits <- unique(as.character(unlist(metabo_object$path.hits[[i]])))
      hits <- hits[!is.na(hits) & hits != ""]
      if (length(hits) == 0) return(NULL)

      hit_mz <- lapply(hits, function(compound_id) {
        vals <- metabo_object$cpd2mz_dict[[compound_id]]
        vals <- suppressWarnings(as.numeric(vals))
        vals[is.finite(vals)]
      })
      data.frame(
        source = "mummichog",
        pathway_id = paste0("P", i),
        pathway = metabo_object$path.nms[[i]],
        compound_id = rep(hits, lengths(hit_mz)),
        compound_name = NA_character_,
        mz = unlist(hit_mz, use.names = FALSE),
        stringsAsFactors = FALSE
      )
    }))
  }

  rows <- rows[!vapply(rows, is.null, logical(1))]
  if (length(rows) == 0) {
    stop("No pathway hits were found in metabo_object.", call. = FALSE)
  }

  hits <- unique(do.call(rbind, rows))
  if (!is.null(pathway)) {
    keep <- tolower(hits$pathway) == tolower(pathway) |
      tolower(hits$pathway_id) == tolower(pathway)
    if (!any(keep)) {
      keep <- grepl(pathway, hits$pathway, ignore.case = TRUE, fixed = TRUE)
    }
    hits <- hits[keep, , drop = FALSE]
    if (nrow(hits) == 0) {
      stop("No pathway matched: ", pathway, call. = FALSE)
    }
  }

  rownames(hits) <- NULL
  hits
}

link_metaboanalyst_pathway_features <- function(metabo_object,
                                                pathway,
                                                pixel_matrix = NULL,
                                                feature_metadata = NULL,
                                                feature_col = "feature",
                                                mz_col = "mz",
                                                compound_id_col = "compound_id",
                                                compound_name_col = "compound_name") {
  hits <- extract_metaboanalyst_pathway_hits(metabo_object, pathway)

  if (!is.null(feature_metadata)) {
    if (!feature_col %in% names(feature_metadata)) {
      stop("feature_metadata is missing feature_col: ", feature_col, call. = FALSE)
    }

    link_data <- hits
    link_data$feature <- NA_character_
    if (compound_id_col %in% names(feature_metadata) && "compound_id" %in% names(link_data)) {
      id_map <- feature_metadata[, c(compound_id_col, feature_col), drop = FALSE]
      names(id_map) <- c("compound_id", "feature")
      link_data <- merge(link_data, id_map, by = "compound_id", all.x = TRUE, sort = FALSE)
      link_data$feature <- ifelse(is.na(link_data$feature.x), link_data$feature.y, link_data$feature.x)
      link_data$feature.x <- NULL
      link_data$feature.y <- NULL
    }
    if (compound_name_col %in% names(feature_metadata) && any(is.na(link_data$feature))) {
      name_map <- feature_metadata[, c(compound_name_col, feature_col), drop = FALSE]
      names(name_map) <- c("compound_name", "feature_by_name")
      link_data <- merge(link_data, name_map, by = "compound_name", all.x = TRUE, sort = FALSE)
      link_data$feature <- ifelse(is.na(link_data$feature), link_data$feature_by_name, link_data$feature)
      link_data$feature_by_name <- NULL
    }
    if (mz_col %in% names(feature_metadata) && any(is.na(link_data$feature))) {
      mz_map <- feature_metadata[, c(mz_col, feature_col), drop = FALSE]
      names(mz_map) <- c("mz", "feature_by_mz")
      link_data <- merge(link_data, mz_map, by = "mz", all.x = TRUE, sort = FALSE)
      link_data$feature <- ifelse(is.na(link_data$feature), link_data$feature_by_mz, link_data$feature)
      link_data$feature_by_mz <- NULL
    }
    hits <- link_data
  } else {
    hits$feature <- NA_character_
  }

  if (!is.null(pixel_matrix)) {
    hits$feature <- vapply(seq_len(nrow(hits)), function(i) {
      if (!is.na(hits$feature[i]) && hits$feature[i] %in% feature_columns(pixel_matrix)) {
        return(hits$feature[i])
      }
      if (is.finite(hits$mz[i])) {
        return(resolve_feature_column(pixel_matrix, hits$mz[i]))
      }
      NA_character_
    }, character(1))
  }

  hits <- hits[!is.na(hits$feature) & hits$feature != "", , drop = FALSE]
  hits <- unique(hits)
  if (nrow(hits) == 0) {
    stop(
      "No pathway hits could be linked to mz_* feature columns. ",
      "For ORA/pathway-ID results, pass feature_metadata with compound IDs or names.",
      call. = FALSE
    )
  }

  rownames(hits) <- NULL
  hits
}

pathway_feature_score_data <- function(pixel_matrix,
                                       features,
                                       transform = c("log10", "identity"),
                                       scale_features = FALSE,
                                       aggregate = c("mean", "sum")) {
  transform <- match.arg(transform)
  aggregate <- match.arg(aggregate)
  required_columns(pixel_matrix, c("x", "y"), "Pixel matrix")

  feature_cols <- unique(vapply(features, function(feature) {
    resolve_feature_column(pixel_matrix, feature)
  }, character(1)))
  if (length(feature_cols) == 0) {
    stop("No feature columns were supplied for pathway scoring.", call. = FALSE)
  }

  values <- as.matrix(pixel_matrix[, feature_cols, drop = FALSE])
  if (transform == "log10") values <- log10(values + 1)
  if (isTRUE(scale_features) && ncol(values) > 1) {
    values <- scale(values)
    values[!is.finite(values)] <- 0
  }
  score <- if (aggregate == "sum") rowSums(values, na.rm = TRUE) else rowMeans(values, na.rm = TRUE)

  data.frame(
    pixel_id = if ("pixel_id" %in% names(pixel_matrix)) pixel_matrix$pixel_id else seq_len(nrow(pixel_matrix)),
    x = pixel_matrix$x,
    y = pixel_matrix$y,
    pathway_score = score,
    n_features = length(feature_cols),
    features = paste(feature_cols, collapse = "; "),
    stringsAsFactors = FALSE
  )
}

pathway_sample_score_data <- function(sample_matrix,
                                      features,
                                      feature_weights = NULL,
                                      weight_col = "signed_log2fc",
                                      feature_col = "feature",
                                      sample_transform = c("identity", "log10"),
                                      score_method = c("weighted_mean", "weighted_sum", "mean", "sum"),
                                      score_column = "pathway_score") {
  sample_transform <- match.arg(sample_transform)
  score_method <- match.arg(score_method)
  required_columns(sample_matrix, "sample_id", "Sample matrix")

  feature_cols <- unique(vapply(features, function(feature) {
    resolve_feature_column(sample_matrix, feature)
  }, character(1)))
  if (length(feature_cols) == 0) {
    stop("No feature columns were supplied for pathway sample scoring.", call. = FALSE)
  }

  values <- as.matrix(sample_matrix[, feature_cols, drop = FALSE])
  if (sample_transform == "log10") values <- log10(values + 1)

  weights <- rep(1, length(feature_cols))
  names(weights) <- feature_cols
  if (!is.null(feature_weights)) {
    required_columns(feature_weights, c(feature_col, weight_col), "Feature weights")
    weight_features <- vapply(feature_weights[[feature_col]], function(feature) {
      resolve_feature_column(sample_matrix, feature)
    }, character(1))
    weight_values <- suppressWarnings(as.numeric(feature_weights[[weight_col]]))
    weight_lookup <- stats::setNames(weight_values, weight_features)
    matched <- intersect(feature_cols, names(weight_lookup))
    if (length(matched) > 0) {
      weights[matched] <- weight_lookup[matched]
    }
  }
  weights[!is.finite(weights)] <- 0

  if (score_method %in% c("weighted_mean", "weighted_sum")) {
    weighted <- sweep(values, 2, weights, "*")
    score <- if (score_method == "weighted_sum") {
      rowSums(weighted, na.rm = TRUE)
    } else {
      denom <- sum(abs(weights), na.rm = TRUE)
      if (!is.finite(denom) || denom == 0) denom <- length(weights)
      rowSums(weighted, na.rm = TRUE) / denom
    }
  } else {
    score <- if (score_method == "sum") rowSums(values, na.rm = TRUE) else rowMeans(values, na.rm = TRUE)
  }

  out <- data.frame(
    Sample = sample_matrix$sample_id,
    stringsAsFactors = FALSE
  )
  out[[score_column]] <- score
  attr(out, "features") <- feature_cols
  attr(out, "weights") <- weights
  out
}

plot_pathway_score_map <- function(pixel_matrix,
                                   features = NULL,
                                   pathway_name = NULL,
                                   transform = c("log10", "identity"),
                                   scale_features = FALSE,
                                   aggregate = c("mean", "sum"),
                                   sample_matrix = NULL,
                                   sample_mapping = NULL,
                                   feature_weights = NULL,
                                   weight_col = "signed_log2fc",
                                   feature_col = "feature",
                                   roi_column = "roi_id",
                                   sample_transform = c("identity", "log10"),
                                   score_method = c("weighted_mean", "weighted_sum", "mean", "sum"),
                                   map_level = c("sample", "roi")) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for plotting.", call. = FALSE)
  }

  transform <- match.arg(transform)
  aggregate <- match.arg(aggregate)
  sample_transform <- match.arg(sample_transform)
  score_method <- match.arg(score_method)
  map_level <- match.arg(map_level)

  if (is.null(features)) {
    if (!is.null(feature_weights) && feature_col %in% names(feature_weights)) {
      features <- feature_weights[[feature_col]]
    } else {
      stop("features must be supplied unless feature_weights contains feature_col.", call. = FALSE)
    }
  }

  if (!is.null(sample_matrix) || !is.null(sample_mapping) || !is.null(feature_weights)) {
    if (is.null(sample_matrix) || is.null(sample_mapping)) {
      stop(
        "sample_matrix and sample_mapping must both be supplied for sample-level pathway score mapping.",
        call. = FALSE
      )
    }
    score_data <- pathway_sample_score_data(
      sample_matrix = sample_matrix,
      features = features,
      feature_weights = feature_weights,
      weight_col = weight_col,
      feature_col = feature_col,
      sample_transform = sample_transform,
      score_method = score_method,
      score_column = "pathway_score"
    )
    mapped <- backmap_sample_scores(score_data, sample_mapping, "pathway_score")
    required_columns(pixel_matrix, c("pixel_id", "x", "y"), "Pixel matrix")
    if (map_level == "roi") {
      required_columns(sample_mapping, c("sample_id", roi_column), "Sample mapping")
      sample_scores <- merge(
        sample_mapping[, c("sample_id", roi_column), drop = FALSE],
        score_data,
        by.x = "sample_id",
        by.y = "Sample",
        all.x = FALSE,
        sort = FALSE
      )
      roi_scores <- stats::aggregate(
        sample_scores$pathway_score,
        by = list(roi_id = sample_scores[[roi_column]]),
        FUN = mean,
        na.rm = TRUE
      )
      names(roi_scores) <- c(roi_column, "pathway_score")
      if (roi_column %in% names(pixel_matrix)) {
        plot_data <- merge(
          pixel_matrix[, c("pixel_id", "x", "y", roi_column), drop = FALSE],
          roi_scores,
          by = roi_column,
          all.x = FALSE,
          sort = FALSE
        )
      } else {
        mapped_roi <- merge(
          mapped,
          sample_mapping[, c("sample_id", roi_column), drop = FALSE],
          by = "sample_id",
          all.x = TRUE,
          sort = FALSE
        )
        mapped_roi <- merge(
          mapped_roi[, c("pixel_id", roi_column), drop = FALSE],
          roi_scores,
          by = roi_column,
          all.x = FALSE,
          sort = FALSE
        )
        plot_data <- merge(
          pixel_matrix[, c("pixel_id", "x", "y"), drop = FALSE],
          mapped_roi[, c("pixel_id", roi_column, "pathway_score"), drop = FALSE],
          by = "pixel_id",
          all.x = FALSE,
          sort = FALSE
        )
      }
    } else {
    plot_data <- merge(
      pixel_matrix[, c("pixel_id", "x", "y"), drop = FALSE],
      mapped,
      by = "pixel_id",
      all.x = FALSE,
      sort = FALSE
    )
    if ("score" %in% names(plot_data) && !"pathway_score" %in% names(plot_data)) {
      names(plot_data)[names(plot_data) == "score"] <- "pathway_score"
    }
    }
    plot_data$n_features <- length(attr(score_data, "features"))
    plot_data$features <- paste(attr(score_data, "features"), collapse = "; ")
  } else {
    plot_data <- pathway_feature_score_data(
      pixel_matrix,
      features,
      transform = transform,
      scale_features = scale_features,
      aggregate = aggregate
    )
  }
  score_midpoint <- stats::median(plot_data$pathway_score, na.rm = TRUE)
  if (!is.finite(score_midpoint)) score_midpoint <- 0
  if (is.null(pathway_name)) pathway_name <- "Pathway feature score"

  ggplot2::ggplot(plot_data, ggplot2::aes(x = .data[["x"]], y = .data[["y"]], fill = .data[["pathway_score"]])) +
    ggplot2::geom_raster() +
    ggplot2::coord_fixed() +
    ggplot2::scale_fill_gradient2(
      low = "#2563eb",
      mid = "white",
      high = "#dc2626",
      midpoint = score_midpoint
    ) +
    ggplot2::labs(
      title = pathway_name,
      subtitle = paste0(length(unique(unlist(strsplit(plot_data$features[1], "; ", fixed = TRUE)))), " linked feature(s)"),
      x = "x",
      y = "y",
      fill = "Score"
    ) +
    ggplot2::theme_minimal()
}

plot_pathway_feature_images <- function(pixel_matrix,
                                        features,
                                        pathway_name = NULL,
                                        transform = c("log10", "identity"),
                                        ncol = 3) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for plotting.", call. = FALSE)
  }
  if (!requireNamespace("viridis", quietly = TRUE)) {
    stop("Package 'viridis' is required for plotting.", call. = FALSE)
  }

  transform <- match.arg(transform)
  required_columns(pixel_matrix, c("x", "y"), "Pixel matrix")
  feature_cols <- unique(vapply(features, function(feature) {
    resolve_feature_column(pixel_matrix, feature)
  }, character(1)))
  rows <- lapply(feature_cols, function(column_name) {
    values <- pixel_matrix[[column_name]]
    if (transform == "log10") values <- log10(values + 1)
    data.frame(
      x = pixel_matrix$x,
      y = pixel_matrix$y,
      feature = column_name,
      intensity = values,
      stringsAsFactors = FALSE
    )
  })
  plot_data <- do.call(rbind, rows)
  if (is.null(pathway_name)) pathway_name <- "Linked pathway feature ion images"

  ggplot2::ggplot(plot_data, ggplot2::aes(x = .data[["x"]], y = .data[["y"]], fill = .data[["intensity"]])) +
    ggplot2::geom_raster() +
    ggplot2::coord_fixed() +
    ggplot2::facet_wrap(~feature, ncol = ncol) +
    viridis::scale_fill_viridis(option = "viridis") +
    ggplot2::labs(title = pathway_name, x = "x", y = "y", fill = "Intensity") +
    ggplot2::theme_minimal() +
    ggplot2::theme(strip.text = ggplot2::element_text(size = 8))
}

plot_metabo_feature_rank <- function(result_data,
                                     top_n = 20,
                                     metric = NULL,
                                     decreasing = NULL,
                                     source_name = "") {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for plotting.", call. = FALSE)
  }

  result_data <- normalize_metaboanalyst_result(result_data, source_name = source_name)
  result_type <- detect_metaboanalyst_result(result_data)
  required_columns(result_data, "Feature", "MetaboAnalyst result")

  if (is.null(metric)) {
    metric <- if ("VIP" %in% names(result_data)) {
      "VIP"
    } else if ("FDR" %in% names(result_data)) {
      "FDR"
    } else if ("p.value" %in% names(result_data)) {
      "p.value"
    } else {
      numeric_cols <- names(result_data)[vapply(result_data, is.numeric, logical(1))]
      if (length(numeric_cols) == 0) {
        stop("No numeric metric column found in MetaboAnalyst result.", call. = FALSE)
      }
      numeric_cols[1]
    }
  }

  required_columns(result_data, metric, "MetaboAnalyst result")
  if (is.null(decreasing)) {
    decreasing <- !(metric %in% c("p.value", "FDR"))
  }

  keep <- is.finite(result_data[[metric]])
  ranked <- result_data[keep, , drop = FALSE]
  ranked <- ranked[order(ranked[[metric]], decreasing = decreasing), , drop = FALSE]
  ranked <- utils::head(ranked, max(1, as.integer(top_n)))
  ranked$Feature <- factor(ranked$Feature, levels = rev(ranked$Feature))

  ggplot2::ggplot(ranked, ggplot2::aes(x = .data[[metric]], y = .data[["Feature"]])) +
    ggplot2::geom_col(fill = "#3b82f6", width = 0.75) +
    ggplot2::labs(
      title = paste0("Top ", nrow(ranked), " MetaboAnalyst features"),
      subtitle = paste0("Detected result type: ", result_type),
      x = metric,
      y = NULL
    ) +
    ggplot2::theme_minimal()
}

plot_metabo_score_plot <- function(score_data,
                                   sample_matrix = NULL,
                                   x = NULL,
                                   y = NULL,
                                   source_name = "") {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for plotting.", call. = FALSE)
  }

  score_data <- normalize_metaboanalyst_result(score_data, source_name = source_name)
  required_columns(score_data, "Sample", "Score data")

  score_columns <- names(score_data)[grepl("^(PC|Comp)[0-9]+$", names(score_data))]
  if (length(score_columns) < 2) {
    stop("Score data must contain at least two score columns, such as PC1 and PC2.", call. = FALSE)
  }
  if (is.null(x)) x <- score_columns[1]
  if (is.null(y)) y <- score_columns[2]
  required_columns(score_data, c(x, y), "Score data")

  plot_data <- score_data
  if (!"Group" %in% names(plot_data)) {
    if (!is.null(sample_matrix) && all(c("sample_id", "cluster") %in% names(sample_matrix))) {
      group_data <- sample_matrix[, c("sample_id", "cluster"), drop = FALSE]
      plot_data <- merge(
        plot_data,
        group_data,
        by.x = "Sample",
        by.y = "sample_id",
        all.x = TRUE,
        sort = FALSE
      )
      plot_data$Group <- paste0("Region_", plot_data$cluster)
    } else {
      plot_data$Group <- "Sample"
    }
  }

  ggplot2::ggplot(plot_data, ggplot2::aes(x = .data[[x]], y = .data[[y]], color = .data[["Group"]], label = .data[["Sample"]])) +
    ggplot2::geom_point(size = 3, alpha = 0.9) +
    ggplot2::labs(
      title = "MetaboAnalyst score plot",
      x = x,
      y = y,
      color = "Group"
    ) +
    ggplot2::theme_minimal()
}

plot_sample_score_map <- function(score_data,
                                  sample_mapping,
                                  pixel_matrix,
                                  score_column = NULL,
                                  source_name = "") {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for plotting.", call. = FALSE)
  }

  score_data <- normalize_metaboanalyst_result(score_data, source_name = source_name)
  score_columns <- names(score_data)[grepl("^(PC|Comp)[0-9]+$", names(score_data))]
  if (length(score_columns) == 0) {
    stop("Score data must contain a score column, such as PC1 or Comp1.", call. = FALSE)
  }
  if (is.null(score_column)) score_column <- score_columns[1]

  mapped <- backmap_sample_scores(score_data, sample_mapping, score_column)
  required_columns(pixel_matrix, c("pixel_id", "x", "y"), "Pixel matrix")
  plot_data <- merge(
    pixel_matrix[, c("pixel_id", "x", "y"), drop = FALSE],
    mapped,
    by = "pixel_id",
    all.x = FALSE,
    sort = FALSE
  )
  score_midpoint <- stats::median(plot_data$score, na.rm = TRUE)
  if (!is.finite(score_midpoint)) score_midpoint <- 0

  x_steps <- diff(sort(unique(plot_data$x)))
  y_steps <- diff(sort(unique(plot_data$y)))
  tile_width <- stats::median(x_steps[x_steps > 0], na.rm = TRUE)
  tile_height <- stats::median(y_steps[y_steps > 0], na.rm = TRUE)
  if (!is.finite(tile_width)) tile_width <- 1
  if (!is.finite(tile_height)) tile_height <- 1

  ggplot2::ggplot(plot_data, ggplot2::aes(x = .data[["x"]], y = .data[["y"]], fill = .data[["score"]])) +
    ggplot2::geom_tile(width = tile_width, height = tile_height) +
    ggplot2::coord_fixed() +
    ggplot2::scale_fill_gradient2(
      low = "#2563eb",
      mid = "white",
      high = "#dc2626",
      midpoint = score_midpoint
    ) +
    ggplot2::labs(
      title = paste0("Spatial back-map: ", score_column),
      x = "x",
      y = "y",
      fill = score_column
    ) +
    ggplot2::theme_minimal()
}

plot_metabo_feature_view <- function(result_data,
                                     pixel_matrix,
                                     feature = NULL,
                                     rank = 1,
                                     metric = NULL,
                                     top_n = 20,
                                     source_name = "",
                                     transform = c("identity", "log10")) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for plotting.", call. = FALSE)
  }
  if (!requireNamespace("grid", quietly = TRUE)) {
    stop("Package 'grid' is required for arranging plots.", call. = FALSE)
  }

  transform <- match.arg(transform)
  result_data <- normalize_metaboanalyst_result(result_data, source_name = source_name)
  required_columns(result_data, "Feature", "MetaboAnalyst result")

  if (is.null(metric)) {
    metric <- if ("VIP" %in% names(result_data)) {
      "VIP"
    } else if ("FDR" %in% names(result_data)) {
      "FDR"
    } else if ("p.value" %in% names(result_data)) {
      "p.value"
    } else {
      numeric_cols <- names(result_data)[vapply(result_data, is.numeric, logical(1))]
      if (length(numeric_cols) == 0) {
        stop("No numeric metric column found in MetaboAnalyst result.", call. = FALSE)
      }
      numeric_cols[1]
    }
  }
  required_columns(result_data, metric, "MetaboAnalyst result")

  decreasing <- !(metric %in% c("p.value", "FDR"))
  ranked <- result_data[is.finite(result_data[[metric]]), , drop = FALSE]
  ranked <- ranked[order(ranked[[metric]], decreasing = decreasing), , drop = FALSE]
  if (nrow(ranked) == 0) {
    stop("No finite metric values found in MetaboAnalyst result.", call. = FALSE)
  }

  if (is.null(feature)) {
    rank <- max(1, min(as.integer(rank), nrow(ranked)))
    feature <- ranked$Feature[rank]
  }

  selected_row <- ranked[ranked$Feature == feature, , drop = FALSE]
  if (nrow(selected_row) == 0) {
    stop("Feature '", feature, "' was not found in the MetaboAnalyst result.", call. = FALSE)
  }
  selected_row <- selected_row[1, , drop = FALSE]

  ranked_top <- utils::head(ranked, max(1, as.integer(top_n)))
  ranked_top$selected <- ranked_top$Feature == feature
  ranked_top$Feature <- factor(ranked_top$Feature, levels = rev(ranked_top$Feature))

  rank_plot <- ggplot2::ggplot(ranked_top, ggplot2::aes(x = .data[[metric]], y = .data[["Feature"]], fill = .data[["selected"]])) +
    ggplot2::geom_col(width = 0.75) +
    ggplot2::scale_fill_manual(values = c("FALSE" = "#3b82f6", "TRUE" = "#dc2626"), guide = "none") +
    ggplot2::labs(
      title = paste0("Top ", nrow(ranked_top), " MetaboAnalyst features"),
      subtitle = paste0("Selected: ", feature),
      x = metric,
      y = NULL
    ) +
    ggplot2::theme_minimal()

  metric_label <- paste0(metric, "=", signif(selected_row[[metric]], 3))
  ion_plot <- plot_ion_image(
    pixel_matrix,
    feature = feature,
    transform = transform,
    title = paste(feature, metric_label)
  )

  grid::grid.newpage()
  grid::pushViewport(grid::viewport(layout = grid::grid.layout(1, 2, widths = grid::unit(c(0.42, 0.58), "npc"))))
  grid::pushViewport(grid::viewport(layout.pos.row = 1, layout.pos.col = 1))
  grid::grid.draw(ggplot2::ggplotGrob(rank_plot))
  grid::popViewport()
  grid::pushViewport(grid::viewport(layout.pos.row = 1, layout.pos.col = 2))
  grid::grid.draw(ggplot2::ggplotGrob(ion_plot))
  grid::popViewport(2)

  invisible(list(
    feature = feature,
    metric = metric,
    result_row = selected_row,
    rank_plot = rank_plot,
    ion_plot = ion_plot
  ))
}

plot_feature_region_boxplot <- function(sample_matrix,
                                        feature,
                                        region_column = "cluster",
                                        region_prefix = "Region_",
                                        title = NULL) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for plotting.", call. = FALSE)
  }

  required_columns(sample_matrix, region_column, "Sample matrix")
  column_name <- resolve_feature_column(sample_matrix, feature)
  region <- sample_matrix[[region_column]]
  if (identical(region_column, "cluster")) {
    region <- paste0(region_prefix, region)
  }

  plot_data <- data.frame(
    Region = factor(region),
    intensity = sample_matrix[[column_name]]
  )

  if (is.null(title)) {
    title <- paste0("Sample-level distribution: ", column_name)
  }

  ggplot2::ggplot(plot_data, ggplot2::aes(x = .data[["Region"]], y = .data[["intensity"]], color = .data[["Region"]])) +
    ggplot2::geom_boxplot(width = 0.55, outlier.shape = NA, fill = NA, linewidth = 0.8) +
    ggplot2::geom_jitter(width = 0.12, height = 0, size = 2.2, alpha = 0.9) +
    ggplot2::labs(
      title = title,
      x = "Region",
      y = column_name,
      color = "Region"
    ) +
    ggplot2::theme_minimal()
}

write_pipeline_outputs <- function(output_dir,
                                   pixel_matrix = NULL,
                                   feature_mapping = NULL,
                                   section_mapping = NULL,
                                   background_stats = NULL,
                                   selected_features = NULL,
                                   reduced_matrix = NULL,
                                   preprocessed_matrix = NULL,
                                   clustered_matrix = NULL,
                                   sample_matrix = NULL,
                                   sample_mapping = NULL,
                                   metaboanalyst_data = NULL,
                                   coordinates = NULL,
                                   feature_metadata = NULL,
                                   qc_summary = NULL,
                                   processing_parameters = NULL,
                                   provenance = NULL) {
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  }

  outputs <- list(
    pixel_feature_matrix.csv = pixel_matrix,
    coordinates.csv = coordinates,
    feature_metadata.csv = feature_metadata,
    feature_mapping.csv = feature_mapping,
    section_mapping.csv = section_mapping,
    background_stats.csv = background_stats,
    selected_features.csv = selected_features,
    reduced_pixel_matrix.csv = reduced_matrix,
    preprocessed_matrix.csv = preprocessed_matrix,
    clustered_matrix.csv = clustered_matrix,
    sample_matrix.csv = sample_matrix,
    sample_mapping.csv = sample_mapping,
    metaboanalyst_data.csv = metaboanalyst_data
  )

  flatten_values <- function(x, prefix = "") {
    if (is.null(x)) return(data.frame(key = character(), value = character(), stringsAsFactors = FALSE))
    if (is.data.frame(x)) return(x)
    if (!is.list(x)) {
      return(data.frame(key = prefix, value = paste(as.character(x), collapse = ";"), stringsAsFactors = FALSE))
    }
    rows <- lapply(names(x), function(name) {
      key <- if (nzchar(prefix)) paste(prefix, name, sep = ".") else name
      flatten_values(x[[name]], key)
    })
    rows <- Filter(function(item) nrow(item) > 0L, rows)
    if (!length(rows)) data.frame(key = character(), value = character(), stringsAsFactors = FALSE) else do.call(rbind, rows)
  }

  outputs[["qc_summary.csv"]] <- if (is.null(qc_summary) || is.data.frame(qc_summary)) qc_summary else flatten_values(qc_summary)
  outputs[["processing_parameters.csv"]] <- if (is.null(processing_parameters) || is.data.frame(processing_parameters)) processing_parameters else flatten_values(processing_parameters)
  outputs[["provenance_manifest.csv"]] <- if (is.null(provenance) || is.data.frame(provenance)) provenance else flatten_values(provenance)

  written <- character()
  for (file_name in names(outputs)) {
    data <- outputs[[file_name]]
    if (is.null(data)) next
    path <- file.path(output_dir, file_name)
    utils::write.csv(data, path, row.names = FALSE)
    written <- c(written, path)
  }

  written
}

make_pipeline_manifest <- function(input_files,
                                   input_type,
                                   ion_mode = NA_character_,
                                   ppm = NA_real_,
                                   parameters = list(),
                                   package_name = "SpatialOmicsMSI") {
  input_files <- normalizePath(input_files, mustWork = TRUE)
  info <- file.info(input_files)
  package_version <- tryCatch(
    as.character(utils::packageVersion(package_name)),
    error = function(e) NA_character_
  )
  manifest <- data.frame(
    record_type = "input_file",
    key = basename(input_files),
    value = input_files,
    size_bytes = as.numeric(info$size),
    md5 = unname(tools::md5sum(input_files)),
    stringsAsFactors = FALSE
  )
  metadata <- data.frame(
    record_type = "run_metadata",
    key = c("input_type", "ion_mode", "ppm", "R_version", "package_version"),
    value = c(input_type, ion_mode, as.character(ppm), R.version.string, package_version),
    size_bytes = NA_real_,
    md5 = NA_character_,
    stringsAsFactors = FALSE
  )
  parameter_rows <- if (length(parameters)) {
    data.frame(
      record_type = "processing_parameter",
      key = names(parameters),
      value = vapply(parameters, function(x) paste(as.character(x), collapse = ";"), character(1)),
      size_bytes = NA_real_,
      md5 = NA_character_,
      stringsAsFactors = FALSE
    )
  } else {
    manifest[FALSE, ]
  }
  rbind(manifest, metadata, parameter_rows)
}

write_validation_plots <- function(output_dir,
                                   pixel_matrix = NULL,
                                   sample_matrix = NULL,
                                   sample_mapping = NULL,
                                   metaboanalyst_data = NULL,
                                   feature = NULL,
                                   width = 8,
                                   height = 6,
                                   dpi = 150) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for validation plots.", call. = FALSE)
  }

  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  }
  plot_dir <- file.path(output_dir, "plots")
  if (!dir.exists(plot_dir)) {
    dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
  }

  written <- character()

  if (!is.null(pixel_matrix) && nrow(pixel_matrix) > 0 && all(c("x", "y") %in% names(pixel_matrix))) {
    fcols <- feature_columns(pixel_matrix)
    if (length(fcols) > 0) {
      pixel_feature <- if (is.null(feature)) fcols[1] else resolve_feature_column(pixel_matrix, feature)
      ion_plot <- plot_ion_image(
        pixel_matrix,
        feature = pixel_feature,
        transform = "log10",
        title = paste0("Ion image: ", pixel_feature)
      )
      path <- file.path(plot_dir, "ion_image.png")
      ggplot2::ggsave(path, ion_plot, width = width, height = height, dpi = dpi)
      written <- c(written, path)
    }
  }

  if (!is.null(sample_mapping) && nrow(sample_mapping) > 0 &&
      !is.null(pixel_matrix) && all(c("pixel_id", "x", "y") %in% names(pixel_matrix))) {
    sample_region <- if ("matched_region_label" %in% names(sample_mapping)) {
      sample_mapping[, c("sample_id", "matched_region_label", "pixel_ids"), drop = FALSE]
    } else if ("roi_id" %in% names(sample_mapping)) {
      tmp <- sample_mapping[, c("sample_id", "roi_id", "pixel_ids"), drop = FALSE]
      tmp$matched_region_label <- tmp$roi_id
      tmp[, c("sample_id", "matched_region_label", "pixel_ids"), drop = FALSE]
    } else if ("cluster" %in% names(sample_mapping)) {
      tmp <- sample_mapping[, c("sample_id", "cluster", "pixel_ids"), drop = FALSE]
      tmp$matched_region_label <- paste0("Region_", tmp$cluster)
      tmp[, c("sample_id", "matched_region_label", "pixel_ids"), drop = FALSE]
    } else {
      NULL
    }

    if (!is.null(sample_region)) {
      mapped <- do.call(rbind, lapply(seq_len(nrow(sample_region)), function(i) {
        ids <- strsplit(as.character(sample_region$pixel_ids[i]), ",", fixed = TRUE)[[1]]
        data.frame(
          pixel_id = as.integer(ids),
          sample_id = sample_region$sample_id[i],
          region = sample_region$matched_region_label[i],
          stringsAsFactors = FALSE
        )
      }))
      plot_data <- merge(
        pixel_matrix[, c("pixel_id", "x", "y"), drop = FALSE],
        mapped,
        by = "pixel_id",
        all.x = FALSE,
        sort = FALSE
      )
      region_plot <- ggplot2::ggplot(plot_data, ggplot2::aes(x = .data[["x"]], y = .data[["y"]], fill = .data[["region"]])) +
        ggplot2::geom_raster() +
        ggplot2::coord_fixed() +
        ggplot2::labs(title = "Constructed spatial samples", x = "x", y = "y", fill = "Region") +
        ggplot2::theme_minimal()
      path <- file.path(plot_dir, "sample_regions.png")
      ggplot2::ggsave(path, region_plot, width = width, height = height, dpi = dpi)
      written <- c(written, path)
    }
  }

  if (!is.null(sample_matrix) && nrow(sample_matrix) > 0) {
    fcols <- feature_columns(sample_matrix)
    if (length(fcols) > 0) {
      sample_feature <- if (is.null(feature)) fcols[1] else resolve_feature_column(sample_matrix, feature)
      group_values <- if ("matched_region_label" %in% names(sample_matrix)) {
        sample_matrix$matched_region_label
      } else if ("cluster" %in% names(sample_matrix)) {
        paste0("Region_", sample_matrix$cluster)
      } else {
        "Sample"
      }
      box_data <- data.frame(
        group = group_values,
        intensity = sample_matrix[[sample_feature]],
        stringsAsFactors = FALSE
      )
      box_plot <- ggplot2::ggplot(box_data, ggplot2::aes(x = .data[["group"]], y = .data[["intensity"]], color = .data[["group"]])) +
        ggplot2::geom_boxplot(outlier.shape = NA) +
        ggplot2::geom_jitter(width = 0.12, height = 0, size = 2) +
        ggplot2::labs(title = paste0("Sample-level feature distribution: ", sample_feature), x = "Group", y = sample_feature) +
        ggplot2::theme_minimal()
      path <- file.path(plot_dir, "sample_feature_boxplot.png")
      ggplot2::ggsave(path, box_plot, width = width, height = height, dpi = dpi)
      written <- c(written, path)
    }
  }

  if (!is.null(metaboanalyst_data) && nrow(metaboanalyst_data) > 0 && "Group" %in% names(metaboanalyst_data)) {
    group_counts <- as.data.frame(table(metaboanalyst_data$Group), stringsAsFactors = FALSE)
    names(group_counts) <- c("Group", "Samples")
    count_plot <- ggplot2::ggplot(group_counts, ggplot2::aes(x = .data[["Group"]], y = .data[["Samples"]], fill = .data[["Group"]])) +
      ggplot2::geom_col(width = 0.7) +
      ggplot2::labs(title = "MetaboAnalyst group sample counts", x = "Group", y = "Samples") +
      ggplot2::theme_minimal()
    path <- file.path(plot_dir, "metaboanalyst_group_counts.png")
    ggplot2::ggsave(path, count_plot, width = width, height = height, dpi = dpi)
    written <- c(written, path)
  }

  written
}
