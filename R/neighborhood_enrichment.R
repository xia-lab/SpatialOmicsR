# Cell-type neighborhood enrichment ----------------------------------------

compute_neighborhood_enrichment <- function(cell_type_labels,
                                            neighbors,
                                            n_perm = 1000,
                                            seed = NULL,
                                            region_id = rep("region_1", length(cell_type_labels)),
                                            p_adjust_method = "BH",
                                            return_permutations = FALSE) {
  if (!inherits(neighbors, "spatial_neighbors")) {
    stop("neighbors must come from build_spatial_neighbors().", call. = FALSE)
  }
  if (!isTRUE(neighbors$symmetric)) {
    stop("Neighborhood enrichment requires a symmetric neighbor graph.", call. = FALSE)
  }
  labels <- as.character(cell_type_labels)
  n_total <- length(labels)
  if (n_total != length(neighbors$x) || length(region_id) != n_total) {
    stop("cell_type_labels and region_id must have one value per graph node.", call. = FALSE)
  }
  if (length(n_perm) != 1L || !is.finite(n_perm) || n_perm < 2 || n_perm != as.integer(n_perm)) {
    stop("n_perm must be one integer of at least two.", call. = FALSE)
  }
  if (length(return_permutations) != 1L || is.na(return_permutations)) {
    stop("return_permutations must be TRUE or FALSE.", call. = FALSE)
  }
  n_perm <- as.integer(n_perm)
  valid <- !is.na(labels) & nzchar(labels)
  region <- as.character(region_id)
  valid <- valid & !is.na(region) & nzchar(region)
  cell_types <- sort(unique(labels[valid]))
  if (length(cell_types) < 2L) {
    stop("At least two non-missing cell types are required.", call. = FALSE)
  }
  edges <- neighbors$edges
  edges <- edges[valid[edges$from] & valid[edges$to], , drop = FALSE]
  if (!nrow(edges)) {
    stop("No valid directed adjacency entries remain after label filtering.", call. = FALSE)
  }
  if (any(region[edges$from] != region[edges$to])) {
    stop("The neighbor graph contains cross-region edges; build graphs separately by region.", call. = FALSE)
  }
  type_factor <- function(x) factor(x, levels = cell_types)
  count_type_pairs <- function(current_labels) {
    unclass(table(
      type_factor(current_labels[edges$from]),
      type_factor(current_labels[edges$to])
    ))
  }
  observed <- count_type_pairs(labels)
  dimnames(observed) <- list(cell_types, cell_types)
  n_types <- length(cell_types)
  sum_perm <- matrix(0, n_types, n_types)
  sumsq_perm <- matrix(0, n_types, n_types)
  upper_extreme <- matrix(0L, n_types, n_types)
  lower_extreme <- matrix(0L, n_types, n_types)
  permutations <- if (isTRUE(return_permutations)) {
    array(0, dim = c(n_perm, n_types, n_types),
      dimnames = list(NULL, cell_types, cell_types))
  } else NULL
  shuffle_groups <- split(which(valid), region[valid])
  if (!is.null(seed)) set.seed(seed)
  for (permutation_index in seq_len(n_perm)) {
    shuffled <- labels
    for (indices in shuffle_groups) {
      shuffled[indices] <- sample(labels[indices], length(indices), replace = FALSE)
    }
    permuted <- count_type_pairs(shuffled)
    sum_perm <- sum_perm + permuted
    sumsq_perm <- sumsq_perm + permuted^2
    upper_extreme <- upper_extreme + (permuted >= observed)
    lower_extreme <- lower_extreme + (permuted <= observed)
    if (isTRUE(return_permutations)) permutations[permutation_index, , ] <- permuted
  }
  permutation_mean <- sum_perm / n_perm
  permutation_variance <- sumsq_perm / n_perm - permutation_mean^2
  permutation_variance[permutation_variance < 0] <- 0
  permutation_sd <- sqrt(permutation_variance)
  zero_sd <- permutation_sd <= sqrt(.Machine$double.eps)
  z_score <- (observed - permutation_mean) / permutation_sd
  z_score[zero_sd] <- NA_real_
  p_enriched <- (upper_extreme + 1) / (n_perm + 1)
  p_depleted <- (lower_extreme + 1) / (n_perm + 1)
  p_two_sided <- 2 * pmin(p_enriched, p_depleted)
  p_two_sided[p_two_sided > 1] <- 1
  dimnames(permutation_mean) <- dimnames(permutation_sd) <- dimnames(z_score) <-
    dimnames(zero_sd) <- dimnames(p_enriched) <- dimnames(p_depleted) <-
    dimnames(p_two_sided) <- list(cell_types, cell_types)
  pair_index <- expand.grid(
    from_cell_type = cell_types,
    to_cell_type = cell_types,
    stringsAsFactors = FALSE
  )
  matrix_index <- cbind(
    match(pair_index$from_cell_type, cell_types),
    match(pair_index$to_cell_type, cell_types)
  )
  pair_index$observed_count <- observed[matrix_index]
  pair_index$permutation_mean <- permutation_mean[matrix_index]
  pair_index$permutation_sd <- permutation_sd[matrix_index]
  pair_index$z_score <- z_score[matrix_index]
  pair_index$p_enriched <- p_enriched[matrix_index]
  pair_index$p_depleted <- p_depleted[matrix_index]
  pair_index$p_two_sided <- p_two_sided[matrix_index]
  pair_index$zero_permutation_sd <- zero_sd[matrix_index]
  # A symmetric graph contains duplicate orientations, so adjust only the
  # unique upper triangle to avoid pretending these are independent tests.
  unique_pair <- match(pair_index$from_cell_type, cell_types) <=
    match(pair_index$to_cell_type, cell_types)
  pair_index$adj_p_two_sided <- NA_real_
  pair_index$adj_p_two_sided[unique_pair] <- stats::p.adjust(
    pair_index$p_two_sided[unique_pair], method = p_adjust_method
  )
  reverse_key <- paste(pair_index$to_cell_type, pair_index$from_cell_type, sep = "\r")
  forward_key <- paste(pair_index$from_cell_type, pair_index$to_cell_type, sep = "\r")
  pair_index$adj_p_two_sided[!unique_pair] <- pair_index$adj_p_two_sided[
    match(reverse_key[!unique_pair], forward_key)
  ]
  list(
    z_score = z_score,
    observed_count = observed,
    permutation_mean = permutation_mean,
    permutation_sd = permutation_sd,
    p_enriched = p_enriched,
    p_depleted = p_depleted,
    p_two_sided = p_two_sided,
    zero_permutation_sd = zero_sd,
    pair_table = pair_index,
    cell_types = cell_types,
    permutations = permutations,
    settings = list(
      n_perm = n_perm, seed = seed, p_adjust_method = p_adjust_method,
      shuffle_scope = "within_region", directed_adjacency_count = TRUE
    ),
    interpretation = paste(
      "This is a global cell-type adjacency enrichment summary, not a per-position niche label.",
      "Positive z-scores indicate more adjacency entries than under within-region label permutations;",
      "zero-variance null entries remain NA."
    )
  )
}

group_cell_types_from_enrichment <- function(enrichment_result,
                                             k_groups,
                                             method = "ward.D2",
                                             zero_sd_action = c("error", "zero"),
                                             exclude_diagonal = FALSE) {
  zero_sd_action <- match.arg(zero_sd_action)
  if (!is.list(enrichment_result) || is.null(enrichment_result$z_score) ||
      is.null(enrichment_result$cell_types)) {
    stop("enrichment_result must come from compute_neighborhood_enrichment().", call. = FALSE)
  }
  profiles <- as.matrix(enrichment_result$z_score)
  cell_types <- as.character(enrichment_result$cell_types)
  if (!identical(dim(profiles), c(length(cell_types), length(cell_types)))) {
    stop("The enrichment z-score matrix dimensions are inconsistent with cell_types.", call. = FALSE)
  }
  if (isTRUE(exclude_diagonal)) diag(profiles) <- 0
  if (any(!is.finite(profiles))) {
    if (zero_sd_action == "error") {
      stop(
        "Non-finite z-scores cannot be clustered. Inspect zero_permutation_sd or set zero_sd_action='zero' explicitly.",
        call. = FALSE
      )
    }
    profiles[!is.finite(profiles)] <- 0
  }
  if (length(k_groups) != 1L || !is.finite(k_groups) || k_groups < 2 ||
      k_groups > length(cell_types) || k_groups != as.integer(k_groups)) {
    stop("k_groups must be an integer from two through the number of cell types.", call. = FALSE)
  }
  distance <- stats::dist(profiles, method = "euclidean")
  if (any(!is.finite(distance))) stop("Interaction-profile distances are not finite.", call. = FALSE)
  fit <- stats::hclust(distance, method = method)
  group <- stats::cutree(fit, k = as.integer(k_groups))
  list(
    groups = data.frame(
      cell_type = cell_types,
      interaction_profile_group = unname(group),
      stringsAsFactors = FALSE
    ),
    fit = fit,
    distance = distance,
    profiles = profiles,
    settings = list(
      k_groups = as.integer(k_groups), method = method,
      zero_sd_action = zero_sd_action, exclude_diagonal = isTRUE(exclude_diagonal)
    ),
    interpretation = paste(
      "Groups summarize similarity of whole enrichment profiles for heatmap ordering.",
      "They are not spatial niches and do not assign positions to a niche."
    )
  )
}
