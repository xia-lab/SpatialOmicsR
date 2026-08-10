# Cell-type and morphology-label inputs for niche analysis -----------------

prepare_cell_type_labels <- function(labels,
                                     proportions = NULL,
                                     required_fields = c("object_id", "x", "y", "cell_type"),
                                     entity_type = NULL,
                                     label_mode = NULL,
                                     coordinate_system = NULL) {
  if (!is.data.frame(labels)) stop("labels must be a data frame.", call. = FALSE)
  required_columns(labels, required_fields, "Cell type label table")
  if (!nrow(labels)) {
    labels$confidence <- numeric()
    labels$source_method <- character()
  }
  if (!is.numeric(labels$x) || !is.numeric(labels$y) ||
      any(!is.finite(labels$x)) || any(!is.finite(labels$y))) {
    stop("Label x and y coordinates must be finite numeric values.", call. = FALSE)
  }
  if ("object_id" %in% names(labels) && (anyNA(labels$object_id) || anyDuplicated(labels$object_id))) {
    stop("object_id must be non-missing and unique.", call. = FALSE)
  }
  if (anyNA(labels$cell_type) || any(!nzchar(as.character(labels$cell_type)))) {
    stop("cell_type must not contain missing or empty labels; use 'unresolved' explicitly.", call. = FALSE)
  }
  if (!"confidence" %in% names(labels)) labels$confidence <- NA_real_
  if (!is.numeric(labels$confidence) ||
      any(is.finite(labels$confidence) & (labels$confidence < 0 | labels$confidence > 1))) {
    stop("confidence must contain NA or numeric values between zero and one.", call. = FALSE)
  }
  if (!"source_method" %in% names(labels)) labels$source_method <- "unspecified"
  if (!"region_id" %in% names(labels)) labels$region_id <- "region_1"
  if (anyNA(labels$region_id) || any(!nzchar(as.character(labels$region_id)))) {
    stop("region_id must not contain missing or empty values.", call. = FALSE)
  }
  if (!is.null(entity_type)) labels$entity_type <- as.character(entity_type)[1]
  if (!"entity_type" %in% names(labels)) labels$entity_type <- "unspecified"
  if (!is.null(label_mode)) labels$label_mode <- as.character(label_mode)[1]
  if (!"label_mode" %in% names(labels)) labels$label_mode <- if (is.null(proportions)) "hard_label" else "soft_proportion"
  if (!is.null(coordinate_system)) labels$coordinate_system <- as.character(coordinate_system)[1]
  if (!"coordinate_system" %in% names(labels)) labels$coordinate_system <- "unspecified"

  if (!is.null(proportions)) {
    proportions <- as.matrix(proportions)
    storage.mode(proportions) <- "double"
    if (nrow(proportions) != nrow(labels) || !ncol(proportions) || is.null(colnames(proportions)) ||
        any(!nzchar(colnames(proportions))) || anyDuplicated(colnames(proportions))) {
      stop("proportions must have one row per label and unique, non-empty cell-type columns.", call. = FALSE)
    }
    if (any(!is.finite(proportions)) || any(proportions < 0)) {
      stop("proportions must contain finite, non-negative values.", call. = FALSE)
    }
    totals <- rowSums(proportions)
    if (any(abs(totals - 1) > sqrt(.Machine$double.eps))) {
      stop("Every proportions row must sum to one.", call. = FALSE)
    }
  }
  structure(
    list(labels = labels, proportions = proportions),
    class = "cell_type_labels"
  )
}

deconvolve_spatial_transcriptomics <- function(count_matrix,
                                               spot_coordinates,
                                               reference_signature,
                                               method = c("nnls", "cosine_label_transfer"),
                                               min_confidence = 0.1,
                                               min_shared_genes = 10,
                                               region_column = NULL,
                                               coordinate_system = "spatial_transcriptomics") {
  method <- match.arg(method)
  required_columns(spot_coordinates, c("spot_id", "x", "y", region_column), "Spot coordinates")
  count_matrix <- as.matrix(count_matrix)
  reference_signature <- as.matrix(reference_signature)
  storage.mode(count_matrix) <- "double"
  storage.mode(reference_signature) <- "double"
  if (!nrow(count_matrix) || !ncol(count_matrix) || !nrow(reference_signature) || !ncol(reference_signature)) {
    stop("count_matrix and reference_signature must be non-empty matrices.", call. = FALSE)
  }
  if (is.null(colnames(count_matrix)) || is.null(rownames(reference_signature)) ||
      is.null(colnames(reference_signature)) || anyDuplicated(colnames(reference_signature))) {
    stop("Gene names and unique reference cell-type column names are required.", call. = FALSE)
  }
  if (any(!is.finite(count_matrix)) || any(count_matrix < 0) ||
      any(!is.finite(reference_signature)) || any(reference_signature < 0)) {
    stop("Expression and reference matrices must contain finite, non-negative linear-scale values.", call. = FALSE)
  }
  if (nrow(count_matrix) != nrow(spot_coordinates)) {
    stop("count_matrix must contain one row per spot coordinate.", call. = FALSE)
  }
  if (!is.null(rownames(count_matrix))) {
    if (anyDuplicated(rownames(count_matrix)) || !setequal(rownames(count_matrix), as.character(spot_coordinates$spot_id))) {
      stop("count_matrix row names and spot_id must match one-to-one.", call. = FALSE)
    }
    count_matrix <- count_matrix[match(as.character(spot_coordinates$spot_id), rownames(count_matrix)), , drop = FALSE]
  }
  shared_genes <- intersect(colnames(count_matrix), rownames(reference_signature))
  if (length(shared_genes) < min_shared_genes) {
    stop("Too few shared genes: ", length(shared_genes), "; require at least ", min_shared_genes, ".", call. = FALSE)
  }
  observed <- count_matrix[, shared_genes, drop = FALSE]
  reference <- reference_signature[shared_genes, , drop = FALSE]
  reference_totals <- colSums(reference)
  if (any(reference_totals <= 0)) {
    stop("Every reference cell type must have positive signal among shared genes.", call. = FALSE)
  }
  reference <- sweep(reference, 2, reference_totals, "/")
  observed_totals <- rowSums(observed)
  if (any(observed_totals <= 0)) {
    stop("Every spatial spot must have positive signal among shared genes.", call. = FALSE)
  }
  observed <- observed / observed_totals
  cell_types <- colnames(reference)
  scores <- matrix(0, nrow(observed), length(cell_types), dimnames = list(NULL, cell_types))
  residual_norm <- numeric(nrow(observed))
  for (i in seq_len(nrow(observed))) {
    target <- as.numeric(observed[i, ])
    if (method == "nnls") {
      if (!requireNamespace("nnls", quietly = TRUE)) {
        stop("Package 'nnls' is required for NNLS deconvolution.", call. = FALSE)
      }
      coefficients <- nnls::nnls(reference, target)$x
    } else {
      target_norm <- sqrt(sum(target^2))
      coefficients <- apply(reference, 2, function(profile) {
        sum(target * profile) / (target_norm * sqrt(sum(profile^2)))
      })
    }
    coefficients[!is.finite(coefficients) | coefficients < 0] <- 0
    total <- sum(coefficients)
    if (total <= 0) {
      stop("No positive cell-type score was obtained for spot ", spot_coordinates$spot_id[i], ".", call. = FALSE)
    }
    scores[i, ] <- coefficients / total
    fitted <- as.numeric(reference %*% scores[i, ])
    residual_norm[i] <- sqrt(sum((target - fitted)^2))
  }
  ordered_scores <- t(apply(scores, 1, sort, decreasing = TRUE))
  confidence <- ordered_scores[, 1]
  margin <- if (ncol(scores) > 1L) ordered_scores[, 1] - ordered_scores[, 2] else rep(1, nrow(scores))
  entropy <- -rowSums(ifelse(scores > 0, scores * log(scores), 0))
  dominant <- cell_types[max.col(scores, ties.method = "first")]
  dominant[confidence < min_confidence] <- "unresolved"
  region_id <- if (is.null(region_column)) rep("region_1", nrow(spot_coordinates)) else as.character(spot_coordinates[[region_column]])
  labels <- data.frame(
    object_id = as.character(spot_coordinates$spot_id),
    x = spot_coordinates$x, y = spot_coordinates$y, region_id = region_id,
    cell_type = dominant, confidence = confidence,
    confidence_margin = margin, composition_entropy = entropy,
    reconstruction_residual = residual_norm,
    source_method = paste0("spatial_transcriptomics_", method),
    stringsAsFactors = FALSE
  )
  result <- prepare_cell_type_labels(
    labels, proportions = scores, entity_type = "spot",
    label_mode = if (method == "nnls") "soft_proportion" else "similarity_score",
    coordinate_system = coordinate_system
  )
  result$shared_genes <- shared_genes
  result$settings <- list(method = method, min_confidence = min_confidence, min_shared_genes = min_shared_genes)
  if (method == "cosine_label_transfer") {
    result$interpretation <- "Normalized cosine similarities are label-transfer scores, not estimated cell fractions."
  }
  result
}

segment_cells_from_histology <- function(image_path,
                                         tissue_mask = NULL,
                                         min_nucleus_size = 20,
                                         max_nucleus_size = 500,
                                         n_morphology_classes = 3,
                                         opening_size = 3,
                                         watershed_tolerance = 1,
                                         watershed_ext = 1,
                                         seed = 1,
                                         region_id = "region_1",
                                         coordinate_system = "histology_pixel") {
  if (!requireNamespace("EBImage", quietly = TRUE)) {
    stop("Package 'EBImage' (Bioconductor) is required for histology segmentation.", call. = FALSE)
  }
  if (!file.exists(image_path)) stop("Histology image does not exist: ", image_path, call. = FALSE)
  if (min_nucleus_size < 1 || max_nucleus_size < min_nucleus_size) {
    stop("Nucleus size bounds must be positive and ordered.", call. = FALSE)
  }
  if (n_morphology_classes < 2 || n_morphology_classes != as.integer(n_morphology_classes)) {
    stop("n_morphology_classes must be an integer of at least two.", call. = FALSE)
  }
  image <- EBImage::readImage(image_path)
  gray <- if (EBImage::colorMode(image) == EBImage::Color) EBImage::channel(image, "luminance") else image
  if (length(dim(gray)) != 2L) stop("image_path must contain one two-dimensional image frame.", call. = FALSE)
  threshold <- EBImage::otsu(gray)
  binary_mask <- gray < threshold
  if (!is.null(tissue_mask)) {
    if (!identical(dim(binary_mask), dim(tissue_mask))) stop("tissue_mask dimensions must match the image.", call. = FALSE)
    if (anyNA(tissue_mask)) stop("tissue_mask must not contain missing values.", call. = FALSE)
    binary_mask <- binary_mask & as.logical(tissue_mask)
  }
  brush <- EBImage::makeBrush(as.integer(opening_size), shape = "disc")
  cleaned_mask <- EBImage::fillHull(EBImage::opening(binary_mask, brush))
  labeled_nuclei <- EBImage::watershed(
    EBImage::distmap(cleaned_mask), tolerance = watershed_tolerance, ext = watershed_ext
  )
  shape <- EBImage::computeFeatures.shape(labeled_nuclei)
  if (is.null(shape) || !nrow(shape)) {
    empty <- data.frame(object_id = character(), x = numeric(), y = numeric(),
      region_id = character(), cell_type = character(), confidence = numeric(),
      source_method = character(), stringsAsFactors = FALSE)
    result <- prepare_cell_type_labels(empty, entity_type = "nucleus_morphology",
      label_mode = "morphology_cluster", coordinate_system = coordinate_system)
    result$qc <- data.frame(n_detected_total = 0L, n_passed_filter = 0L, threshold = threshold)
    result$segmentation <- list(binary_mask = binary_mask, cleaned_mask = cleaned_mask, labeled_nuclei = labeled_nuclei)
    return(result)
  }
  moment <- EBImage::computeFeatures.moment(labeled_nuclei)
  intensity <- EBImage::computeFeatures.basic(labeled_nuclei, gray)
  size_ok <- shape[, "s.area"] >= min_nucleus_size & shape[, "s.area"] <= max_nucleus_size
  if (sum(size_ok) < n_morphology_classes) {
    stop("Too few nuclei pass size filtering to form ", n_morphology_classes, " classes.", call. = FALSE)
  }
  perimeter <- shape[size_ok, "s.perimeter"]
  circularity <- 4 * pi * shape[size_ok, "s.area"] / perimeter^2
  circularity[!is.finite(circularity)] <- 0
  features <- cbind(
    area = shape[size_ok, "s.area"], circularity = circularity,
    mean_intensity = intensity[size_ok, "b.mean"], intensity_sd = intensity[size_ok, "b.sd"]
  )
  scaled <- scale(features)
  scaled[!is.finite(scaled)] <- 0
  if (nrow(unique(as.data.frame(scaled))) < n_morphology_classes) {
    stop("Too few distinct morphology profiles for the requested class count.", call. = FALSE)
  }
  if (!is.null(seed)) set.seed(seed)
  fit <- stats::kmeans(scaled, centers = as.integer(n_morphology_classes), nstart = 25)
  class_order <- order(tapply(features[, "area"], fit$cluster, mean))
  label_map <- stats::setNames(paste0("morphology_class_", seq_len(n_morphology_classes)), class_order)
  morphology_label <- unname(label_map[as.character(fit$cluster)])
  distances <- vapply(seq_len(nrow(scaled)), function(i) {
    sqrt(rowSums((fit$centers - matrix(scaled[i, ], nrow(fit$centers), ncol(scaled), byrow = TRUE))^2))
  }, numeric(n_morphology_classes))
  distances <- t(distances)
  ordered_distance <- t(apply(distances, 1, sort))
  confidence <- ordered_distance[, 2] / (ordered_distance[, 1] + ordered_distance[, 2] + .Machine$double.eps)
  labels <- data.frame(
    object_id = paste0("nucleus_", seq_len(nrow(shape))[size_ok]),
    x = moment[size_ok, "m.cx"], y = moment[size_ok, "m.cy"], region_id = as.character(region_id)[1],
    cell_type = morphology_label, confidence = confidence,
    nucleus_area = features[, "area"], circularity = features[, "circularity"],
    mean_intensity = features[, "mean_intensity"], intensity_sd = features[, "intensity_sd"],
    source_method = "watershed_morphology_kmeans", stringsAsFactors = FALSE
  )
  result <- prepare_cell_type_labels(labels, entity_type = "nucleus_morphology",
    label_mode = "morphology_cluster", coordinate_system = coordinate_system)
  result$features <- features
  result$fit <- fit
  result$qc <- data.frame(
    n_detected_total = nrow(shape), n_passed_filter = sum(size_ok),
    n_removed_by_size = sum(!size_ok), otsu_threshold = threshold,
    min_nucleus_size = min_nucleus_size, max_nucleus_size = max_nucleus_size
  )
  result$segmentation <- list(binary_mask = binary_mask, cleaned_mask = cleaned_mask, labeled_nuclei = labeled_nuclei)
  result$interpretation <- "Classes are exploratory nuclear-morphology clusters, not biological cell types."
  result
}

transfer_cell_type_to_msi_pixels <- function(pixel_matrix,
                                             cell_type_labels,
                                             method = c("shared_coordinates", "nearest_neighbor"),
                                             max_distance = NULL,
                                             pixel_region_column = "section_id",
                                             registration = NULL) {
  method <- match.arg(method)
  required_columns(pixel_matrix, c("pixel_id", "x", "y"), "Pixel matrix")
  label_object <- if (inherits(cell_type_labels, "cell_type_labels")) cell_type_labels else prepare_cell_type_labels(cell_type_labels)
  labels <- label_object$labels
  if (!is.null(registration)) {
    labels <- transform_histology_coordinates(labels, registration,
      x_column = "x", y_column = "y", output_columns = c("x", "y"),
      section_column = if (identical(pixel_region_column, "section_id")) "region_id" else pixel_region_column)
  }
  pixel_region <- if (pixel_region_column %in% names(pixel_matrix)) as.character(pixel_matrix[[pixel_region_column]]) else rep("region_1", nrow(pixel_matrix))
  label_region <- as.character(labels$region_id)
  if (method == "nearest_neighbor" &&
      (is.null(max_distance) || length(max_distance) != 1L || !is.finite(max_distance) || max_distance <= 0)) {
    stop("nearest_neighbor transfer requires one positive finite max_distance in MSI coordinate units.", call. = FALSE)
  }
  assigned <- rep(NA_integer_, nrow(pixel_matrix))
  match_distance <- rep(NA_real_, nrow(pixel_matrix))
  if (method == "shared_coordinates") {
    label_key <- paste(label_region, labels$x, labels$y, sep = "\r")
    if (anyDuplicated(label_key)) stop("Label coordinates must be unique within region.", call. = FALSE)
    pixel_key <- paste(pixel_region, pixel_matrix$x, pixel_matrix$y, sep = "\r")
    assigned <- match(pixel_key, label_key)
    match_distance[!is.na(assigned)] <- 0
  } else {
    for (region in unique(pixel_region)) {
      pixel_idx <- which(pixel_region == region)
      label_idx <- which(label_region == region)
      if (!length(label_idx)) next
      for (i in pixel_idx) {
        distance <- sqrt((labels$x[label_idx] - pixel_matrix$x[i])^2 + (labels$y[label_idx] - pixel_matrix$y[i])^2)
        local <- which.min(distance)
        if (distance[local] <= max_distance) {
          assigned[i] <- label_idx[local]
          match_distance[i] <- distance[local]
        }
      }
    }
  }
  out <- pixel_matrix
  matched <- !is.na(assigned)
  out$cell_type <- NA_character_
  out$cell_type_confidence <- NA_real_
  out$cell_type_source_method <- NA_character_
  out$cell_type_match_distance <- match_distance
  out$cell_type[matched] <- as.character(labels$cell_type[assigned[matched]])
  out$cell_type_confidence[matched] <- labels$confidence[assigned[matched]]
  out$cell_type_source_method[matched] <- as.character(labels$source_method[assigned[matched]])
  transferred_proportions <- NULL
  if (!is.null(label_object$proportions)) {
    transferred_proportions <- matrix(NA_real_, nrow(pixel_matrix), ncol(label_object$proportions),
      dimnames = list(NULL, colnames(label_object$proportions)))
    transferred_proportions[matched, ] <- label_object$proportions[assigned[matched], , drop = FALSE]
  }
  list(
    matrix = out, proportions = transferred_proportions,
    match_table = data.frame(pixel_id = pixel_matrix$pixel_id,
      label_object_id = ifelse(matched, as.character(labels$object_id[assigned]), NA_character_),
      distance = match_distance, matched = matched, stringsAsFactors = FALSE),
    qc = data.frame(method = method, n_pixels = nrow(pixel_matrix), n_matched = sum(matched),
      matched_fraction = mean(matched), max_distance = if (is.null(max_distance)) NA_real_ else max_distance),
    settings = list(pixel_region_column = pixel_region_column, registration_used = !is.null(registration))
  )
}

compute_neighborhood_composition_soft <- function(proportions,
                                                  x,
                                                  y,
                                                  k = 10,
                                                  pixel_id = seq_len(nrow(proportions)),
                                                  region_id = rep("region_1", nrow(proportions)),
                                                  include_self = TRUE) {
  proportions <- as.matrix(proportions)
  storage.mode(proportions) <- "double"
  if (!nrow(proportions) || !ncol(proportions) || is.null(colnames(proportions)) ||
      any(!is.finite(proportions)) || any(proportions < 0) ||
      any(abs(rowSums(proportions) - 1) > sqrt(.Machine$double.eps))) {
    stop("proportions must be a finite non-negative matrix with named columns and rows summing to one.", call. = FALSE)
  }
  graph <- compute_neighborhood_composition(
    rep(".__soft__", nrow(proportions)), x, y, k = k,
    pixel_id = pixel_id, region_id = region_id, include_self = include_self
  )
  neighborhood <- t(vapply(seq_len(nrow(proportions)), function(i) {
    colMeans(proportions[graph$neighbor_indices[i, ], , drop = FALSE])
  }, numeric(ncol(proportions))))
  safe_names <- paste0("composition__", make.unique(make.names(colnames(proportions))))
  colnames(neighborhood) <- safe_names
  dominant <- colnames(proportions)[max.col(proportions, ties.method = "first")]
  matrix_output <- data.frame(pixel_id = pixel_id, region_id = as.character(region_id),
    focal_cell_type = dominant, window_size = k, neighborhood,
    check.names = FALSE, stringsAsFactors = FALSE)
  structure(list(
    matrix = matrix_output, counts = NULL, neighbor_indices = graph$neighbor_indices,
    cell_type_mapping = data.frame(cell_type = colnames(proportions), composition_column = safe_names),
    settings = list(k = k, include_self = include_self, input_mode = "soft_proportion")
  ), class = "neighborhood_composition")
}
