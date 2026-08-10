# Cell-type proximity and ligand-receptor sensitivity analyses ------------

.lr_local_seed <- function(seed, code) {
  if (is.null(seed)) return(force(code))
  if (length(seed) != 1L || !is.finite(seed)) {
    stop("seed must be NULL or one finite value.", call. = FALSE)
  }
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  on.exit({
    if (had_seed) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  set.seed(as.integer(seed))
  force(code)
}

.lr_validate_common <- function(expression_matrix, cell_type_labels, permutation_strata = NULL) {
  expression_matrix <- as.matrix(expression_matrix)
  storage.mode(expression_matrix) <- "double"
  if (!nrow(expression_matrix) || !ncol(expression_matrix)) {
    stop("expression_matrix must contain at least one row and one gene.", call. = FALSE)
  }
  if (is.null(colnames(expression_matrix)) || anyNA(colnames(expression_matrix)) ||
      any(!nzchar(colnames(expression_matrix))) || anyDuplicated(colnames(expression_matrix))) {
    stop("expression_matrix must have unique, non-empty gene column names.", call. = FALSE)
  }
  if (any(!is.finite(expression_matrix))) {
    stop("expression_matrix must contain only finite values.", call. = FALSE)
  }
  cell_type_labels <- as.character(cell_type_labels)
  if (length(cell_type_labels) != nrow(expression_matrix)) {
    stop("cell_type_labels must contain one value per expression-matrix row.", call. = FALSE)
  }
  valid <- !is.na(cell_type_labels) & nzchar(cell_type_labels)
  if (sum(valid) < 2L || length(unique(cell_type_labels[valid])) < 1L) {
    stop("Too few valid cell-type labels.", call. = FALSE)
  }
  if (is.null(permutation_strata)) {
    permutation_strata <- rep("all", nrow(expression_matrix))
    stratified <- FALSE
    warning(
      "No permutation_strata were supplied. Randomization inference is limited to the ",
      "observed tissue field and is not population-level biological replication.",
      call. = FALSE
    )
  } else {
    if (is.data.frame(permutation_strata) || is.matrix(permutation_strata)) {
      if (nrow(permutation_strata) != nrow(expression_matrix)) {
        stop("permutation_strata must have one row per expression-matrix row.", call. = FALSE)
      }
      strata_frame <- as.data.frame(permutation_strata, stringsAsFactors = FALSE)
      permutation_strata <- do.call(
        interaction,
        c(as.list(strata_frame), list(drop = TRUE, lex.order = TRUE, sep = "::"))
      )
    }
    if (length(permutation_strata) != nrow(expression_matrix)) {
      stop("permutation_strata must contain one value per expression-matrix row.", call. = FALSE)
    }
    permutation_strata <- as.character(permutation_strata)
    if (any(is.na(permutation_strata[valid]) | !nzchar(permutation_strata[valid]))) {
      stop("permutation_strata must be non-missing for every labeled node.", call. = FALSE)
    }
    stratified <- TRUE
  }
  list(
    expression = expression_matrix,
    labels = cell_type_labels,
    strata = permutation_strata,
    valid = valid,
    stratified = stratified
  )
}

.lr_finite_mean <- function(x) {
  x <- x[is.finite(x)]
  if (length(x)) mean(x) else NA_real_
}

.lr_validate_iterations <- function(n_perm) {
  if (length(n_perm) != 1L || !is.finite(n_perm) || n_perm < 1L ||
      n_perm != as.integer(n_perm)) {
    stop("n_perm must be one positive integer.", call. = FALSE)
  }
  as.integer(n_perm)
}

.lr_prepare_database <- function(lr_database, available_genes,
                                 ligand_col, receptor_col, complex_delimiter) {
  required_columns(lr_database, c(ligand_col, receptor_col), "LR database")
  if (length(complex_delimiter) != 1L || is.na(complex_delimiter) || !nzchar(complex_delimiter)) {
    stop("complex_delimiter must be one non-empty fixed string.", call. = FALSE)
  }
  split_components <- function(value) {
    value <- trimws(as.character(value))
    components <- trimws(strsplit(value, complex_delimiter, fixed = TRUE)[[1]])
    unique(components[nzchar(components)])
  }
  ligand_names <- as.character(lr_database[[ligand_col]])
  receptor_names <- as.character(lr_database[[receptor_col]])
  keep_name <- !is.na(ligand_names) & nzchar(trimws(ligand_names)) &
    !is.na(receptor_names) & nzchar(trimws(receptor_names))
  pairs <- lr_database[keep_name, , drop = FALSE]
  ligand_names <- trimws(ligand_names[keep_name])
  receptor_names <- trimws(receptor_names[keep_name])
  ligand_components <- lapply(ligand_names, split_components)
  receptor_components <- lapply(receptor_names, split_components)
  available <- vapply(seq_along(ligand_components), function(i) {
    length(ligand_components[[i]]) > 0L && length(receptor_components[[i]]) > 0L &&
      all(ligand_components[[i]] %in% available_genes) &&
      all(receptor_components[[i]] %in% available_genes)
  }, logical(1))
  pairs <- pairs[available, , drop = FALSE]
  ligand_names <- ligand_names[available]
  receptor_names <- receptor_names[available]
  ligand_components <- ligand_components[available]
  receptor_components <- receptor_components[available]
  if (!nrow(pairs)) {
    stop("No complete ligand-receptor pairs match expression_matrix genes.", call. = FALSE)
  }
  key <- paste(ligand_names, receptor_names, sep = "\r")
  unique_pair <- !duplicated(key)
  data.frame_result <- data.frame(
    ligand = ligand_names[unique_pair], receptor = receptor_names[unique_pair],
    stringsAsFactors = FALSE
  )
  list(
    table = data.frame_result,
    ligand_components = ligand_components[unique_pair],
    receptor_components = receptor_components[unique_pair],
    n_input = nrow(lr_database),
    n_complete = sum(unique_pair),
    n_unavailable = nrow(lr_database) - sum(unique_pair)
  )
}

.lr_complex_expression <- function(expression_matrix, components) {
  values <- expression_matrix[, components, drop = FALSE]
  if (ncol(values) == 1L) as.numeric(values[, 1]) else apply(values, 1L, min)
}

.lr_score <- function(ligand_expression, receptor_expression, ligand_nodes,
                      receptor_nodes, expression_threshold, min_expression_fraction) {
  ligand_values <- ligand_expression[ligand_nodes]
  receptor_values <- receptor_expression[receptor_nodes]
  ligand_values <- ligand_values[is.finite(ligand_values)]
  receptor_values <- receptor_values[is.finite(receptor_values)]
  if (!length(ligand_values) || !length(receptor_values)) {
    return(c(score = NA_real_, ligand_mean = NA_real_, receptor_mean = NA_real_,
             ligand_fraction = NA_real_, receptor_fraction = NA_real_))
  }
  ligand_fraction <- mean(ligand_values > expression_threshold)
  receptor_fraction <- mean(receptor_values > expression_threshold)
  ligand_mean <- mean(ligand_values)
  receptor_mean <- mean(receptor_values)
  score <- if (ligand_fraction >= min_expression_fraction &&
               receptor_fraction >= min_expression_fraction) {
    ligand_mean + receptor_mean
  } else {
    0
  }
  c(
    score = score, ligand_mean = ligand_mean, receptor_mean = receptor_mean,
    ligand_fraction = ligand_fraction, receptor_fraction = receptor_fraction
  )
}

.lr_permutation_p <- function(observed, permuted) {
  permuted <- permuted[is.finite(permuted)]
  if (!is.finite(observed) || !length(permuted)) {
    return(c(p_enrichment = NA_real_, p_depletion = NA_real_, p_two_sided = NA_real_))
  }
  upper <- (sum(permuted >= observed) + 1) / (length(permuted) + 1)
  lower <- (sum(permuted <= observed) + 1) / (length(permuted) + 1)
  c(p_enrichment = upper, p_depletion = lower, p_two_sided = min(1, 2 * min(upper, lower)))
}

cell_proximity_enrichment <- function(
    cell_type_labels,
    neighbors,
    permutation_strata = NULL,
    n_perm = 1000L,
    p_adjust_method = "BH",
    cross_stratum_action = c("error", "drop"),
    seed = 1234) {
  cross_stratum_action <- match.arg(cross_stratum_action)
  n_perm <- .lr_validate_iterations(n_perm)
  if (!inherits(neighbors, "spatial_neighbors")) {
    stop("neighbors must come from build_spatial_neighbors().", call. = FALSE)
  }
  labels <- as.character(cell_type_labels)
  n_nodes <- length(labels)
  if (n_nodes != length(neighbors$x)) {
    stop("cell_type_labels must contain one value per spatial-neighbor node.", call. = FALSE)
  }
  valid <- !is.na(labels) & nzchar(labels)
  if (is.null(permutation_strata)) {
    strata <- rep("all", n_nodes)
    stratified <- FALSE
    warning(
      "No permutation_strata were supplied. Proximity inference is limited to the ",
      "observed tissue field and is not population-level biological replication.",
      call. = FALSE
    )
  } else {
    if (is.data.frame(permutation_strata) || is.matrix(permutation_strata)) {
      if (nrow(permutation_strata) != n_nodes) {
        stop("permutation_strata must have one row per node.", call. = FALSE)
      }
      strata_frame <- as.data.frame(permutation_strata, stringsAsFactors = FALSE)
      permutation_strata <- do.call(
        interaction,
        c(as.list(strata_frame), list(drop = TRUE, lex.order = TRUE, sep = "::"))
      )
    }
    if (length(permutation_strata) != n_nodes) {
      stop("permutation_strata must contain one value per node.", call. = FALSE)
    }
    strata <- as.character(permutation_strata)
    if (any(is.na(strata[valid]) | !nzchar(strata[valid]))) {
      stop("permutation_strata must be non-missing for every labeled node.", call. = FALSE)
    }
    stratified <- TRUE
  }
  cell_types <- sort(unique(labels[valid]))
  if (!length(cell_types)) stop("No valid cell-type labels remain.", call. = FALSE)
  edges <- neighbors$undirected_edges
  edges <- edges[valid[edges$from] & valid[edges$to], , drop = FALSE]
  if (!nrow(edges)) stop("No valid spatial edges remain.", call. = FALSE)
  cross_stratum <- strata[edges$from] != strata[edges$to]
  if (any(cross_stratum)) {
    if (cross_stratum_action == "error") {
      stop("Spatial graph contains edges crossing permutation strata.", call. = FALSE)
    }
    edges <- edges[!cross_stratum, , drop = FALSE]
  }
  if (!nrow(edges)) stop("No within-stratum spatial edges remain.", call. = FALSE)

  count_pairs <- function(current_labels) {
    from <- factor(current_labels[edges$from], levels = cell_types)
    to <- factor(current_labels[edges$to], levels = cell_types)
    directional <- as.matrix(table(from, to))
    directional + t(directional) - diag(diag(directional))
  }
  observed <- count_pairs(labels)
  strata_indices <- split(which(valid), strata[valid])
  permuted <- .lr_local_seed(seed, {
    replicate(n_perm, {
      shuffled <- labels
      for (index in strata_indices) shuffled[index] <- sample(labels[index], replace = FALSE)
      count_pairs(shuffled)
    }, simplify = "array")
  })
  if (length(cell_types) == 1L) {
    dim(permuted) <- c(1L, 1L, n_perm)
  }
  expected <- apply(permuted, c(1L, 2L), mean)
  ratio <- observed / expected
  ratio[expected == 0] <- NA_real_
  rows <- list()
  row_index <- 1L
  for (i in seq_along(cell_types)) {
    for (j in i:length(cell_types)) {
      p <- .lr_permutation_p(observed[i, j], permuted[i, j, ])
      rows[[row_index]] <- data.frame(
        type_a = cell_types[i], type_b = cell_types[j],
        observed = observed[i, j], expected = expected[i, j],
        enrichment_ratio = ratio[i, j], log2_enrichment = log2(ratio[i, j]),
        p_enrichment = p["p_enrichment"], p_depletion = p["p_depletion"],
        p_two_sided = p["p_two_sided"], stringsAsFactors = FALSE
      )
      row_index <- row_index + 1L
    }
  }
  result <- do.call(rbind, rows)
  result$adj_p_enrichment <- stats::p.adjust(result$p_enrichment, method = p_adjust_method)
  result$adj_p_depletion <- stats::p.adjust(result$p_depletion, method = p_adjust_method)
  result$adj_p_two_sided <- stats::p.adjust(result$p_two_sided, method = p_adjust_method)
  result <- result[order(result$adj_p_two_sided, -abs(result$log2_enrichment), na.last = TRUE), , drop = FALSE]
  rownames(result) <- NULL
  list(
    enrichment_table = result,
    enrichment_ratio_matrix = ratio,
    observed_matrix = observed,
    expected_matrix = expected,
    settings = list(
      n_perm = n_perm, stratified_permutation = stratified,
      n_strata = length(strata_indices), n_edges = nrow(edges),
      n_cross_stratum_edges_dropped = sum(cross_stratum),
      p_value_resolution = 1 / (n_perm + 1),
      null_model = "cell-type labels shuffled within permutation strata on a fixed graph"
    )
  )
}

expr_cell_cell_communication <- function(
    expression_matrix,
    cell_type_labels,
    lr_database,
    ligand_col = "ligand",
    receptor_col = "receptor",
    complex_delimiter = "|",
    permutation_strata = NULL,
    expression_threshold = 0,
    min_expression_fraction = 0.1,
    min_cells_per_type = 5L,
    n_perm = 1000L,
    p_adjust_method = "BH",
    seed = 1234) {
  common <- .lr_validate_common(expression_matrix, cell_type_labels, permutation_strata)
  n_perm <- .lr_validate_iterations(n_perm)
  if (length(min_cells_per_type) != 1L || !is.finite(min_cells_per_type) ||
      min_cells_per_type < 2L || min_cells_per_type != as.integer(min_cells_per_type)) {
    stop("min_cells_per_type must be one integer of at least two.", call. = FALSE)
  }
  if (length(min_expression_fraction) != 1L || !is.finite(min_expression_fraction) ||
      min_expression_fraction < 0 || min_expression_fraction > 1) {
    stop("min_expression_fraction must lie between zero and one.", call. = FALSE)
  }
  database <- .lr_prepare_database(
    lr_database, colnames(common$expression), ligand_col, receptor_col, complex_delimiter
  )
  cell_types <- sort(unique(common$labels[common$valid]))
  type_indices <- lapply(cell_types, function(type) which(common$valid & common$labels == type))
  names(type_indices) <- cell_types
  strata_indices <- split(which(common$valid), common$strata[common$valid])
  shuffled_labels <- .lr_local_seed(seed, {
    replicate(n_perm, {
      shuffled <- common$labels
      for (index in strata_indices) shuffled[index] <- sample(common$labels[index], replace = FALSE)
      shuffled
    }, simplify = "matrix")
  })
  if (n_perm == 1L) shuffled_labels <- matrix(shuffled_labels, ncol = 1L)

  rows <- list()
  row_index <- 1L
  for (pair_index in seq_len(nrow(database$table))) {
    ligand_expression <- .lr_complex_expression(common$expression, database$ligand_components[[pair_index]])
    receptor_expression <- .lr_complex_expression(common$expression, database$receptor_components[[pair_index]])
    for (ligand_type in cell_types) {
      ligand_nodes <- type_indices[[ligand_type]]
      for (receptor_type in cell_types) {
        receptor_nodes <- type_indices[[receptor_type]]
        eligible <- length(ligand_nodes) >= min_cells_per_type && length(receptor_nodes) >= min_cells_per_type
        observed <- if (eligible) {
          .lr_score(
            ligand_expression, receptor_expression, ligand_nodes, receptor_nodes,
            expression_threshold, min_expression_fraction
          )
        } else {
          c(score = NA_real_, ligand_mean = NA_real_, receptor_mean = NA_real_,
            ligand_fraction = NA_real_, receptor_fraction = NA_real_)
        }
        permuted <- if (eligible) vapply(seq_len(n_perm), function(p) {
          labels_p <- shuffled_labels[, p]
          nodes_l <- which(common$valid & labels_p == ligand_type)
          nodes_r <- which(common$valid & labels_p == receptor_type)
          .lr_score(
            ligand_expression, receptor_expression, nodes_l, nodes_r,
            expression_threshold, min_expression_fraction
          )["score"]
        }, numeric(1)) else rep(NA_real_, n_perm)
        p <- .lr_permutation_p(observed["score"], permuted)
        rows[[row_index]] <- data.frame(
          ligand = database$table$ligand[pair_index], receptor = database$table$receptor[pair_index],
          lig_cell_type = ligand_type, rec_cell_type = receptor_type,
          lig_expr = observed["ligand_mean"], rec_expr = observed["receptor_mean"],
          lig_expression_fraction = observed["ligand_fraction"],
          rec_expression_fraction = observed["receptor_fraction"],
          LR_expr = observed["score"], rand_expr = .lr_finite_mean(permuted),
          av_diff = observed["score"] - .lr_finite_mean(permuted),
          p_enrichment = p["p_enrichment"], p_depletion = p["p_depletion"],
          p_two_sided = p["p_two_sided"],
          lig_nr = length(ligand_nodes), rec_nr = length(receptor_nodes),
          status = if (eligible) "tested" else "too_few_cells",
          stringsAsFactors = FALSE
        )
        row_index <- row_index + 1L
      }
    }
  }
  result <- do.call(rbind, rows)
  result$adj_p_enrichment <- stats::p.adjust(result$p_enrichment, method = p_adjust_method)
  result$adj_p_depletion <- stats::p.adjust(result$p_depletion, method = p_adjust_method)
  result$adj_p_two_sided <- stats::p.adjust(result$p_two_sided, method = p_adjust_method)
  result <- result[order(result$adj_p_enrichment, -result$av_diff, na.last = TRUE), , drop = FALSE]
  rownames(result) <- NULL
  attr(result, "settings") <- list(
    analysis = "expression-only cell-type-specific LR sensitivity analysis",
    n_perm = n_perm, stratified_permutation = common$stratified,
    n_strata = length(strata_indices), p_value_resolution = 1 / (n_perm + 1),
    complex_summary = "minimum expression across required subunits",
    combined_score = "ligand mean plus receptor mean after expression-fraction filtering",
    fdr_scope = "all tested LR-pair by directed-cell-type combinations",
    inference_scope = "conditional randomization within observed tissue fields; not population-level replication",
    n_lr_input = database$n_input, n_lr_complete = database$n_complete
  )
  result
}

spat_cell_cell_communication <- function(
    expression_matrix,
    cell_type_labels,
    neighbors,
    lr_database,
    ligand_col = "ligand",
    receptor_col = "receptor",
    complex_delimiter = "|",
    permutation_strata = NULL,
    expression_threshold = 0,
    min_expression_fraction = 0.1,
    min_interacting_cells = 5L,
    n_perm = 1000L,
    p_adjust_method = "BH",
    cross_stratum_action = c("error", "drop"),
    seed = 1234) {
  cross_stratum_action <- match.arg(cross_stratum_action)
  common <- .lr_validate_common(expression_matrix, cell_type_labels, permutation_strata)
  n_perm <- .lr_validate_iterations(n_perm)
  if (!inherits(neighbors, "spatial_neighbors")) {
    stop("neighbors must come from build_spatial_neighbors().", call. = FALSE)
  }
  if (length(neighbors$x) != nrow(common$expression)) {
    stop("neighbors must contain one node per expression-matrix row.", call. = FALSE)
  }
  if (length(min_interacting_cells) != 1L || !is.finite(min_interacting_cells) ||
      min_interacting_cells < 2L || min_interacting_cells != as.integer(min_interacting_cells)) {
    stop("min_interacting_cells must be one integer of at least two.", call. = FALSE)
  }
  if (length(min_expression_fraction) != 1L || !is.finite(min_expression_fraction) ||
      min_expression_fraction < 0 || min_expression_fraction > 1) {
    stop("min_expression_fraction must lie between zero and one.", call. = FALSE)
  }
  database <- .lr_prepare_database(
    lr_database, colnames(common$expression), ligand_col, receptor_col, complex_delimiter
  )
  edges <- neighbors$undirected_edges
  edges <- edges[common$valid[edges$from] & common$valid[edges$to], , drop = FALSE]
  if (!nrow(edges)) stop("No valid spatial edges remain.", call. = FALSE)
  cross_stratum <- common$strata[edges$from] != common$strata[edges$to]
  if (any(cross_stratum)) {
    if (cross_stratum_action == "error") {
      stop("Spatial graph contains edges crossing permutation strata.", call. = FALSE)
    }
    edges <- edges[!cross_stratum, , drop = FALSE]
  }
  if (!nrow(edges)) stop("No within-stratum spatial edges remain.", call. = FALSE)
  cell_types <- sort(unique(common$labels[common$valid]))
  all_type_stratum <- split(
    which(common$valid),
    interaction(common$labels[common$valid], common$strata[common$valid], drop = TRUE, sep = "\r")
  )
  get_pool <- function(type, stratum) {
    key <- paste(type, stratum, sep = "\r")
    pool <- all_type_stratum[[key]]
    if (is.null(pool)) integer() else pool
  }

  rows <- list()
  row_index <- 1L
  .lr_local_seed(seed, {
    for (pair_index in seq_len(nrow(database$table))) {
      ligand_expression <- .lr_complex_expression(common$expression, database$ligand_components[[pair_index]])
      receptor_expression <- .lr_complex_expression(common$expression, database$receptor_components[[pair_index]])
      for (ligand_type in cell_types) {
        for (receptor_type in cell_types) {
          from_ligand <- common$labels[edges$from] == ligand_type & common$labels[edges$to] == receptor_type
          to_ligand <- common$labels[edges$to] == ligand_type & common$labels[edges$from] == receptor_type
          relevant <- from_ligand | to_ligand
          interaction_edges <- edges[relevant, , drop = FALSE]
          if (nrow(interaction_edges)) {
            ligand_nodes_by_edge <- ifelse(from_ligand[relevant], interaction_edges$from, interaction_edges$to)
            receptor_nodes_by_edge <- ifelse(from_ligand[relevant], interaction_edges$to, interaction_edges$from)
            if (ligand_type == receptor_type) {
              participating <- sort(unique(c(interaction_edges$from, interaction_edges$to)))
              ligand_nodes <- receptor_nodes <- participating
            } else {
              ligand_nodes <- sort(unique(ligand_nodes_by_edge))
              receptor_nodes <- sort(unique(receptor_nodes_by_edge))
            }
          } else {
            ligand_nodes <- receptor_nodes <- integer()
          }
          eligible <- length(ligand_nodes) >= min_interacting_cells &&
            length(receptor_nodes) >= min_interacting_cells
          observed <- if (eligible) {
            .lr_score(
              ligand_expression, receptor_expression, ligand_nodes, receptor_nodes,
              expression_threshold, min_expression_fraction
            )
          } else {
            c(score = NA_real_, ligand_mean = NA_real_, receptor_mean = NA_real_,
              ligand_fraction = NA_real_, receptor_fraction = NA_real_)
          }
          permuted <- rep(NA_real_, n_perm)
          if (eligible) {
            ligand_counts <- table(common$strata[ligand_nodes])
            receptor_counts <- table(common$strata[receptor_nodes])
            permuted <- vapply(seq_len(n_perm), function(p) {
              sampled_ligand <- unlist(lapply(names(ligand_counts), function(stratum) {
                pool <- get_pool(ligand_type, stratum)
                sample(pool, as.integer(ligand_counts[[stratum]]), replace = FALSE)
              }), use.names = FALSE)
              sampled_receptor <- if (ligand_type == receptor_type &&
                                      identical(ligand_counts, receptor_counts)) {
                sampled_ligand
              } else {
                unlist(lapply(names(receptor_counts), function(stratum) {
                  pool <- get_pool(receptor_type, stratum)
                  sample(pool, as.integer(receptor_counts[[stratum]]), replace = FALSE)
                }), use.names = FALSE)
              }
              .lr_score(
                ligand_expression, receptor_expression, sampled_ligand, sampled_receptor,
                expression_threshold, min_expression_fraction
              )["score"]
            }, numeric(1))
          }
          p <- .lr_permutation_p(observed["score"], permuted)
          rows[[row_index]] <- data.frame(
            ligand = database$table$ligand[pair_index], receptor = database$table$receptor[pair_index],
            lig_cell_type = ligand_type, rec_cell_type = receptor_type,
            lig_expr = observed["ligand_mean"], rec_expr = observed["receptor_mean"],
            lig_expression_fraction = observed["ligand_fraction"],
            rec_expression_fraction = observed["receptor_fraction"],
            LR_expr = observed["score"], rand_expr = .lr_finite_mean(permuted),
            av_diff = observed["score"] - .lr_finite_mean(permuted),
            log2fc = log2((observed["score"] + 0.1) / (.lr_finite_mean(permuted) + 0.1)),
            p_enrichment = p["p_enrichment"], p_depletion = p["p_depletion"],
            p_two_sided = p["p_two_sided"],
            lig_nr = length(ligand_nodes), rec_nr = length(receptor_nodes),
            n_edges = nrow(interaction_edges),
            status = if (eligible) "tested" else "too_few_interacting_cells",
            stringsAsFactors = FALSE
          )
          row_index <- row_index + 1L
        }
      }
    }
  })
  result <- do.call(rbind, rows)
  result$adj_p_enrichment <- stats::p.adjust(result$p_enrichment, method = p_adjust_method)
  result$adj_p_depletion <- stats::p.adjust(result$p_depletion, method = p_adjust_method)
  result$adj_p_two_sided <- stats::p.adjust(result$p_two_sided, method = p_adjust_method)
  result$PI <- result$log2fc * -log10(pmax(result$adj_p_enrichment, .Machine$double.xmin))
  result <- result[order(result$adj_p_enrichment, -result$av_diff, na.last = TRUE), , drop = FALSE]
  rownames(result) <- NULL
  attr(result, "settings") <- list(
    analysis = "spatially interacting-cell LR sensitivity analysis",
    n_perm = n_perm, stratified_permutation = common$stratified,
    n_strata = length(unique(common$strata[common$valid])),
    n_edges = nrow(edges), n_cross_stratum_edges_dropped = sum(cross_stratum),
    p_value_resolution = 1 / (n_perm + 1),
    null_model = paste(
      "within each stratum and cell type, sample the same number of cells as the",
      "unique cells participating in observed cross-type spatial edges"
    ),
    complex_summary = "minimum expression across required subunits",
    fdr_scope = "all tested LR-pair by directed-cell-type combinations",
    inference_scope = "conditional randomization within observed tissue fields; not population-level replication",
    n_lr_input = database$n_input, n_lr_complete = database$n_complete
  )
  result
}
