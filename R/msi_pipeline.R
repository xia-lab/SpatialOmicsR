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

compute_morans_i_grid <- function(values,
                                  x,
                                  y,
                                  n_perm = 0,
                                  alternative = c("greater", "two.sided"),
                                  seed = NULL) {
  alternative <- match.arg(alternative)
  keep <- is.finite(values) & is.finite(x) & is.finite(y)
  values <- as.numeric(values[keep])
  x <- x[keep]
  y <- y[keep]
  n <- length(values)
  if (n < 3 || stats::var(values, na.rm = TRUE) == 0) {
    return(list(I = 0, p_value = NA_real_, n = n, n_edges = 0L))
  }

  coords <- data.frame(
    idx = seq_len(n),
    x = x,
    y = y,
    key = paste(x, y, sep = "\r"),
    right_key = paste(x + 1, y, sep = "\r"),
    down_key = paste(x, y + 1, sep = "\r"),
    stringsAsFactors = FALSE
  )
  key_to_idx <- stats::setNames(coords$idx, coords$key)
  right_idx <- unname(key_to_idx[coords$right_key])
  down_idx <- unname(key_to_idx[coords$down_key])

  edge_i <- c(coords$idx[!is.na(right_idx)], coords$idx[!is.na(down_idx)])
  edge_j <- c(right_idx[!is.na(right_idx)], down_idx[!is.na(down_idx)])
  n_edges <- length(edge_i)
  if (n_edges == 0) {
    return(list(I = 0, p_value = NA_real_, n = n, n_edges = 0L))
  }

  moran_stat <- function(v) {
    centered <- v - mean(v, na.rm = TRUE)
    denom <- sum(centered^2, na.rm = TRUE)
    if (!is.finite(denom) || denom == 0) return(0)
    numerator <- 2 * sum(centered[edge_i] * centered[edge_j], na.rm = TRUE)
    w <- 2 * n_edges
    (n / w) * numerator / denom
  }

  observed <- moran_stat(values)
  p_value <- NA_real_
  if (n_perm > 0) {
    if (!is.null(seed)) set.seed(seed)
    permuted <- replicate(n_perm, moran_stat(sample(values, length(values), replace = FALSE)))
    if (alternative == "greater") {
      p_value <- (sum(permuted >= observed, na.rm = TRUE) + 1) / (sum(is.finite(permuted)) + 1)
    } else {
      p_value <- (sum(abs(permuted) >= abs(observed), na.rm = TRUE) + 1) / (sum(is.finite(permuted)) + 1)
    }
  }

  list(I = observed, p_value = p_value, n = n, n_edges = n_edges)
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
                                                   seed = NULL) {
  alternative <- match.arg(alternative)
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

  if (!is.null(seed)) set.seed(seed)
  rows <- lapply(fcols, function(feature) {
    stat <- compute_morans_i_grid(
      values = pixel_matrix[[feature]],
      x = coords[[x_col]],
      y = coords[[y_col]],
      n_perm = n_perm,
      alternative = alternative,
      seed = NULL
    )
    data.frame(
      feature = feature,
      morans_i = stat$I,
      p_value = stat$p_value,
      n = stat$n,
      n_edges = stat$n_edges,
      stringsAsFactors = FALSE
    )
  })
  result <- do.call(rbind, rows)
  result$adj_p_value <- stats::p.adjust(result$p_value, method = p_adjust_method)
  result[order(result$adj_p_value, -abs(result$morans_i)), , drop = FALSE]
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
                                    ppm = 5) {
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
  used_lcms <- rep(FALSE, nrow(lcms_features))
  rows <- list()

  for (i in seq_len(nrow(msi_features))) {
    mz <- msi_mz[i]
    candidate_idx <- which(!used_lcms)
    if (mode_match) {
      candidate_idx <- candidate_idx[lcms_features[[lcms_mode_col]][candidate_idx] == msi_features[[msi_mode_col]][i]]
    }
    if (length(candidate_idx) == 0) next

    ppm_error <- abs(lcms_mz[candidate_idx] - mz) / mz * 1e6
    # Allow a tiny numerical tolerance to account for floating-point rounding
    # when comparing ppm errors to the threshold.
    within <- which(ppm_error <= (ppm + (.Machine$double.eps * 1e6)))
    if (length(within) == 0) next

    chosen <- candidate_idx[within[which.min(ppm_error[within])]]
    used_lcms[chosen] <- TRUE

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

    rows[[length(rows) + 1L]] <- data.frame(
      msi_feature_id = if ("feature_id" %in% names(msi_features)) as.character(msi_features$feature_id[i]) else NA_character_,
      msi_mz = mz,
      lcms_feature_id = if (!is.null(lcms_id_col)) as.character(lcms_features[[lcms_id_col]][chosen]) else NA_character_,
      lcms_mz = lcms_mz[chosen],
      ion_mode = if (mode_match) as.character(msi_features[[msi_mode_col]][i]) else NA_character_,
      ppm_error = ppm_error[within[which.min(ppm_error[within])]],
      msi_log2fc = msi_log2fc,
      lcms_log2fc = lcms_log2fc,
      direction_agreement = direction_agreement,
      match_type = match_type,
      lcms_rt = if (!is.null(lcms_rt_col)) lcms_rt else NA_real_,
      lcms_name = lcms_name,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  }

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
    stringsAsFactors = FALSE,
    check.names = FALSE
  ), class = "data.frame")

  summary <- data.frame(
    total_msi_features = nrow(msi_features),
    matched_features = nrow(result),
    feature_level_supported_matches = sum(result$match_type == "feature_level_orthogonal_support", na.rm = TRUE),
    agreement_rate = if (sum(!is.na(result$direction_agreement)) == 0) NA_real_ else mean(result$direction_agreement, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
  attr(result, "summary") <- summary
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

  if (is.null(title)) title <- column_name

  ggplot2::ggplot(plot_data, ggplot2::aes(x = .data[["x"]], y = .data[["y"]], fill = .data[["intensity"]])) +
    ggplot2::geom_raster() +
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
