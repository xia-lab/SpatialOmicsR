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

  for (i in seq_along(targets)) {
    tolerance <- targets[i] * ppm / 1e6
    hits <- which(abs(spectrum_mz - targets[i]) <= tolerance)
    if (length(hits) > 0) {
      closest <- hits[which.min(abs(spectrum_mz[hits] - targets[i]))]
      out[i] <- spectrum_intensity[closest]
    }
  }

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

extract_target_intensity_list <- function(spectrum_mz, intensity_data, targets, ppm, n_pixels) {
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
  }

  out
}

extract_target_intensity_matrix <- function(msi, targets, ppm = 10) {
  raw_mz <- Cardinal::mz(msi)
  intensity_data <- Cardinal::intensity(msi)
  n_pixels <- nrow(as.data.frame(Cardinal::coord(msi)))

  if (is_list_like_spectra(intensity_data) || is.null(dim(intensity_data))) {
    return(extract_target_intensity_list(raw_mz, intensity_data, targets, ppm, n_pixels))
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
    if (length(hits) == 0) next

    closest <- hits[which.min(abs(spectrum_mz[hits] - targets[i]))]
    values <- if (feature_by_pixel) {
      as_numeric_vector(intensity_data[closest, ])
    } else {
      as_numeric_vector(intensity_data[, closest])
    }
    out[, i] <- values
  }

  out
}

load_msi_target_features <- function(imzml_path,
                                     mzmine_features,
                                     ppm = 10,
                                     section_id = "Section01",
                                     section_order = 1) {
  if (!requireNamespace("Cardinal", quietly = TRUE)) {
    stop("Package 'Cardinal' is required for imzML loading.", call. = FALSE)
  }

  mapping <- make_feature_mapping(mzmine_features)
  msi <- Cardinal::readMSIData(imzml_path)
  coordinates <- as.data.frame(Cardinal::coord(msi))
  required_columns(coordinates, c("x", "y"), "Cardinal coordinates")

  n_pixels <- nrow(coordinates)
  intensity_matrix <- extract_target_intensity_matrix(msi, mapping$mz, ppm = ppm)
  colnames(intensity_matrix) <- mapping$column_name

  pixel_matrix <- data.frame(
    pixel_id = seq_len(n_pixels),
    section_id = section_id,
    section_order = section_order,
    x = coordinates$x,
    y = coordinates$y,
    intensity_matrix,
    check.names = FALSE
  )

  list(pixel_matrix = pixel_matrix, feature_mapping = mapping)
}

load_serial_msi_target_features <- function(imzml_paths,
                                            mzmine_features,
                                            ppm = 10,
                                            section_ids = NULL,
                                            section_order = NULL) {
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
      section_order = section_order[i]
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
                              background_method = c("log_tic_kmeans", "tic_percentile", "median_fraction"),
                              background_percent = 1,
                              background_percentile = 10,
                              min_nonzero_features = 1,
                              do_tic = TRUE,
                              do_log = TRUE,
                              do_scale = FALSE) {
  background_method <- match.arg(background_method)
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
    stringsAsFactors = FALSE
  )

  if (length(fcols) == 0) return(list(matrix = out, distributions = distributions, background_stats = background_stats))

  if (do_background) {
    row_tic <- rowSums(out[fcols], na.rm = TRUE)
    nonzero_features <- rowSums(out[fcols] > 0, na.rm = TRUE)
    positive_tic <- row_tic[row_tic > 0]
    threshold <- NA_real_

    if (background_method == "log_tic_kmeans" && length(unique(row_tic)) > 1) {
      set.seed(42)
      log_tic <- log10(row_tic + 1)
      fit <- stats::kmeans(log_tic, centers = 2, nstart = 25)
      cluster_means <- tapply(log_tic, fit$cluster, mean, na.rm = TRUE)
      tissue_cluster <- as.integer(names(cluster_means)[which.max(cluster_means)])
      keep_rows <- fit$cluster == tissue_cluster & nonzero_features >= min_nonzero_features
      background_stats$tissue_cluster <- tissue_cluster
      background_stats$cluster1_mean_log_tic <- unname(cluster_means[["1"]])
      background_stats$cluster2_mean_log_tic <- unname(cluster_means[["2"]])
    } else {
      threshold <- if (background_method == "tic_percentile" && length(positive_tic) > 0) {
        stats::quantile(positive_tic, probs = background_percentile / 100, na.rm = TRUE)
      } else {
        stats::median(row_tic, na.rm = TRUE) * background_percent / 100
      }
      keep_rows <- row_tic >= threshold & nonzero_features >= min_nonzero_features
    }

    out <- out[keep_rows, , drop = FALSE]
    background_stats$pixels_after <- nrow(out)
    background_stats$pixels_removed <- background_stats$pixels_before - background_stats$pixels_after
    background_stats$removed_fraction <- background_stats$pixels_removed / background_stats$pixels_before
    background_stats$threshold <- threshold
  }
  distributions$background <- unlist(out[fcols], use.names = FALSE)

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

preprocess_sections <- function(pixel_matrix,
                                section_column = "section_id",
                                do_background = TRUE,
                                background_method = c("log_tic_kmeans", "tic_percentile", "median_fraction"),
                                background_percent = 1,
                                background_percentile = 10,
                                min_nonzero_features = 1,
                                do_tic = TRUE,
                                do_log = TRUE,
                                do_scale = FALSE,
                                scale_scope = c("global", "section")) {
  background_method <- match.arg(background_method)
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
                                       background_method = c("log_tic_kmeans", "tic_percentile", "median_fraction"),
                                       background_percent = 1,
                                       background_percentile = 10,
                                       min_nonzero_features = 1,
                                       do_tic = TRUE,
                                       do_log = TRUE) {
  combine_mode <- match.arg(combine_mode)
  background_method <- match.arg(background_method)
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

cluster_pixels <- function(pixel_matrix, k = 3, nstart = 25, pca_components = NULL) {
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
  set.seed(1)
  fit <- stats::kmeans(
    matrix_data,
    centers = k,
    nstart = nstart,
    iter.max = 100,
    algorithm = "Lloyd"
  )
  out <- pixel_matrix
  out$cluster <- fit$cluster
  list(matrix = out, fit = fit, pca = pca, pca_components = if (is.null(pca)) 0 else ncol(matrix_data))
}

cluster_diagnostics <- function(pixel_matrix, max_k = 10, pca_components = NULL) {
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
  if (nrow(matrix_data_full) > 2000) {
    set.seed(1)
    sample_index <- sample(seq_len(nrow(matrix_data_full)), 2000)
    subsampled <- TRUE
  }
  matrix_data_used <- matrix_data_full[sample_index, , drop = FALSE]

  rows <- lapply(2:max_k, function(k) {
    set.seed(1)
    matrix_data <- matrix_data_used
    fit <- stats::kmeans(matrix_data, centers = k, nstart = 10, iter.max = 50, algorithm = "Lloyd")
    sil_mean <- NA_real_
    if (nrow(matrix_data) <= 2000) {
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

sample_subregions <- function(clustered_matrix, grid_size = 5, min_pixels = 30) {
  required_columns(clustered_matrix, c("pixel_id", "x", "y", "cluster"), "Clustered matrix")
  fcols <- feature_columns(clustered_matrix)
  grid_size <- max(2L, as.integer(grid_size))
  min_pixels <- max(1L, as.integer(min_pixels))
  work <- clustered_matrix
  work$grid_x <- as.integer(cut(work$x, breaks = grid_size, labels = FALSE, include.lowest = TRUE))
  work$grid_y <- as.integer(cut(work$y, breaks = grid_size, labels = FALSE, include.lowest = TRUE))
  work$grid_cell <- paste(work$grid_x, work$grid_y, sep = "_")

  split_groups <- split(work, interaction(work$cluster, work$grid_cell, drop = TRUE), drop = TRUE)
  sample_rows <- list()
  mapping_rows <- list()

  for (group_name in names(split_groups)) {
    group_data <- split_groups[[group_name]]
    n_pixels <- nrow(group_data)
    if (n_pixels < min_pixels) next

    cluster_value <- group_data$cluster[1]
    grid_cell <- group_data$grid_cell[1]
    sample_id <- paste0("region", cluster_value, "_", grid_cell)
    means <- colMeans(group_data[fcols], na.rm = TRUE)

    sample_rows[[sample_id]] <- data.frame(
      sample_id = sample_id,
      cluster = cluster_value,
      grid_cell = grid_cell,
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
      cluster = cluster_value,
      grid_cell = grid_cell,
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

make_matched_sample_id <- function(section_id, matched_region_label) {
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
  if (is.null(group_column)) {
    group_column <- if ("cluster" %in% names(sample_matrix)) "cluster" else if ("matched_region_label" %in% names(sample_matrix)) "matched_region_label" else NULL
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

  ggplot2::ggplot(plot_data, ggplot2::aes(x = .data[["x"]], y = .data[["y"]], fill = .data[["score"]])) +
    ggplot2::geom_raster() +
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
                                   metaboanalyst_data = NULL) {
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  }

  outputs <- list(
    pixel_feature_matrix.csv = pixel_matrix,
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
