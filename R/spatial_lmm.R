# Gaussian spatial mixed-effects differential analysis --------------------

differential_region_analysis_spatial_lmm <- function(
    sample_matrix,
    group_column = "roi_id",
    subject_column,
    section_column,
    x_col = NULL,
    y_col = NULL,
    x_bounds = c("x_min", "x_max"),
    y_bounds = c("y_min", "y_max"),
    x_resolution,
    y_resolution,
    distance_unit,
    correlation_structure = c("exponential", "gaussian", "spherical"),
    reference_group = NULL,
    p_adjust_method = "BH",
    min_subjects_per_group = 5,
    min_unique_coordinates_per_field = 4,
    field_action = c("error", "drop"),
    confidence_level = 0.95,
    keep_fits = FALSE) {
  correlation_structure <- match.arg(correlation_structure)
  field_action <- match.arg(field_action)
  if (!requireNamespace("nlme", quietly = TRUE)) {
    stop("Package 'nlme' is required for Gaussian spatial LMM analysis.", call. = FALSE)
  }
  required_columns(sample_matrix, c(group_column, subject_column, section_column), "Sample matrix")
  if (is.null(x_col) != is.null(y_col)) {
    stop("x_col and y_col must either both be supplied or both be NULL.", call. = FALSE)
  }
  if (is.null(x_col)) {
    if (length(x_bounds) != 2L || length(y_bounds) != 2L) {
      stop("x_bounds and y_bounds must each contain minimum and maximum column names.", call. = FALSE)
    }
    required_columns(sample_matrix, c(x_bounds, y_bounds), "Sample matrix")
    coordinate_columns <- c(x_bounds, y_bounds)
    if (any(!vapply(sample_matrix[coordinate_columns], is.numeric, logical(1)))) {
      stop("Coordinate-bound columns must be numeric.", call. = FALSE)
    }
    raw_x <- (sample_matrix[[x_bounds[1]]] + sample_matrix[[x_bounds[2]]]) / 2
    raw_y <- (sample_matrix[[y_bounds[1]]] + sample_matrix[[y_bounds[2]]]) / 2
    coordinate_source <- "bounds_center"
  } else {
    required_columns(sample_matrix, c(x_col, y_col), "Sample matrix")
    if (!is.numeric(sample_matrix[[x_col]]) || !is.numeric(sample_matrix[[y_col]])) {
      stop("x_col and y_col must be numeric.", call. = FALSE)
    }
    raw_x <- sample_matrix[[x_col]]
    raw_y <- sample_matrix[[y_col]]
    coordinate_source <- "explicit_columns"
  }
  positive_scalar <- function(value, name) {
    if (length(value) != 1L || !is.finite(value) || value <= 0) {
      stop(name, " must be supplied as one positive finite value.", call. = FALSE)
    }
    as.numeric(value)
  }
  x_resolution <- positive_scalar(x_resolution, "x_resolution")
  y_resolution <- positive_scalar(y_resolution, "y_resolution")
  if (length(distance_unit) != 1L || is.na(distance_unit) || !nzchar(as.character(distance_unit))) {
    stop("distance_unit must be supplied explicitly.", call. = FALSE)
  }
  if (length(min_subjects_per_group) != 1L || !is.finite(min_subjects_per_group) ||
      min_subjects_per_group < 3 || min_subjects_per_group != as.integer(min_subjects_per_group)) {
    stop("min_subjects_per_group must be an integer of at least three.", call. = FALSE)
  }
  if (length(min_unique_coordinates_per_field) != 1L ||
      !is.finite(min_unique_coordinates_per_field) ||
      min_unique_coordinates_per_field < 3 ||
      min_unique_coordinates_per_field != as.integer(min_unique_coordinates_per_field)) {
    stop("min_unique_coordinates_per_field must be an integer of at least three.", call. = FALSE)
  }
  if (length(confidence_level) != 1L || !is.finite(confidence_level) ||
      confidence_level <= 0 || confidence_level >= 1) {
    stop("confidence_level must lie strictly between zero and one.", call. = FALSE)
  }
  fcols <- feature_columns(sample_matrix)
  if (!length(fcols)) stop("No mz_ feature columns found in sample_matrix.", call. = FALSE)
  if (any(!vapply(sample_matrix[fcols], is.numeric, logical(1)))) {
    stop("All mz_ feature columns must be numeric.", call. = FALSE)
  }
  base <- data.frame(
    row_index = seq_len(nrow(sample_matrix)),
    group = as.character(sample_matrix[[group_column]]),
    subject = as.character(sample_matrix[[subject_column]]),
    section = as.character(sample_matrix[[section_column]]),
    x = raw_x * x_resolution,
    y = raw_y * y_resolution,
    stringsAsFactors = FALSE
  )
  valid_identifier <- function(x) !is.na(x) & nzchar(x)
  if (any(!valid_identifier(base$group)) || any(!valid_identifier(base$subject)) ||
      any(!valid_identifier(base$section))) {
    stop("Group, subject, and section identifiers must be non-missing and non-empty.", call. = FALSE)
  }
  if (any(!is.finite(base$x)) || any(!is.finite(base$y))) {
    stop("Spatial center coordinates must be finite.", call. = FALSE)
  }
  groups <- sort(unique(base$group))
  if (length(groups) < 2L) stop("At least two groups are required.", call. = FALSE)
  if (!is.null(reference_group)) {
    reference_group <- as.character(reference_group)[1]
    if (!reference_group %in% groups) stop("reference_group is absent from group_column.", call. = FALSE)
    groups <- c(reference_group, setdiff(groups, reference_group))
  }
  base$group <- factor(base$group, levels = groups)
  base$subject <- factor(base$subject)
  base$section <- factor(base$section)
  base$spatial_field <- interaction(base$subject, base$section, drop = TRUE, sep = "::")
  subject_counts <- vapply(groups, function(group) {
    length(unique(as.character(base$subject[base$group == group])))
  }, integer(1))
  if (any(subject_counts < min_subjects_per_group)) {
    stop(
      "Every group must contain at least ", min_subjects_per_group,
      " biological subjects. Counts: ",
      paste(names(subject_counts), subject_counts, sep = "=", collapse = ", "),
      call. = FALSE
    )
  }
  field_indices <- split(seq_len(nrow(base)), base$spatial_field)
  field_qc <- do.call(rbind, lapply(names(field_indices), function(field) {
    index <- field_indices[[field]]
    coordinate_key <- paste(base$x[index], base$y[index], sep = "\r")
    data.frame(
      spatial_field = field,
      subject = as.character(base$subject[index[1]]),
      section = as.character(base$section[index[1]]),
      n_observations = length(index),
      n_unique_coordinates = length(unique(coordinate_key)),
      n_duplicate_coordinates = sum(duplicated(coordinate_key)),
      n_groups = length(unique(base$group[index])),
      stringsAsFactors = FALSE
    )
  }))
  field_qc$eligible <- field_qc$n_unique_coordinates >= min_unique_coordinates_per_field &
    field_qc$n_duplicate_coordinates == 0L
  field_qc$exclusion_reason <- ifelse(
    field_qc$n_duplicate_coordinates > 0L,
    "duplicate_coordinates_within_spatial_field",
    ifelse(
      field_qc$n_unique_coordinates < min_unique_coordinates_per_field,
      "too_few_unique_coordinates",
      NA_character_
    )
  )
  if (any(!field_qc$eligible) && field_action == "error") {
    failed <- field_qc$spatial_field[!field_qc$eligible]
    stop(
      "Ineligible spatial field(s): ", paste(failed, collapse = ", "),
      ". Aggregate true technical duplicates or choose field_action='drop' explicitly.",
      call. = FALSE
    )
  }
  retained_fields <- field_qc$spatial_field[field_qc$eligible]
  base <- droplevels(base[as.character(base$spatial_field) %in% retained_fields, , drop = FALSE])
  if (nlevels(base$subject) < min_subjects_per_group) {
    stop("Too few subjects remain after spatial-field filtering.", call. = FALSE)
  }
  retained_subject_counts <- vapply(groups, function(group) {
    length(unique(as.character(base$subject[base$group == group])))
  }, integer(1))
  if (any(retained_subject_counts < min_subjects_per_group)) {
    stop("Spatial-field filtering left too few subjects in at least one group.", call. = FALSE)
  }

  corr_fun <- switch(
    correlation_structure,
    exponential = nlme::corExp,
    gaussian = nlme::corGaus,
    spherical = nlme::corSpher
  )
  correlation_formula <- ~ x + y | subject/section
  random_formula <- ~ 1 | subject/section
  fit_control <- nlme::lmeControl(maxIter = 100, msMaxIter = 100, returnObject = FALSE)
  feature_rows <- vector("list", length(fcols))
  contrast_rows <- list()
  fits <- if (isTRUE(keep_fits)) vector("list", length(fcols)) else NULL
  contrast_index <- 1L

  for (feature_index in seq_along(fcols)) {
    feature <- fcols[feature_index]
    work <- base
    work$value <- as.numeric(sample_matrix[[feature]][work$row_index])
    work <- droplevels(work[is.finite(work$value), , drop = FALSE])
    skip_reason <- NA_character_
    if (nrow(work) < 6L || length(unique(work$value)) < 2L || nlevels(work$group) < 2L) {
      skip_reason <- "insufficient_or_constant_response"
    }
    feature_subject_counts <- vapply(groups, function(group) {
      length(unique(as.character(work$subject[work$group == group])))
    }, integer(1))
    if (is.na(skip_reason) && any(feature_subject_counts < min_subjects_per_group)) {
      skip_reason <- "feature_missingness_reduces_subject_replication"
    }
    feature_field_qc <- if (is.na(skip_reason)) {
      do.call(rbind, lapply(split(seq_len(nrow(work)), work$spatial_field), function(index) {
        data.frame(
          n = length(index),
          n_coordinates = nrow(unique(work[index, c("x", "y"), drop = FALSE]))
        )
      }))
    } else data.frame(n = integer(), n_coordinates = integer())
    if (is.na(skip_reason) && any(feature_field_qc$n_coordinates < min_unique_coordinates_per_field)) {
      skip_reason <- "feature_missingness_reduces_spatial_field"
    }
    warning_messages <- character()
    fit_independence <- fit_spatial <- NULL
    if (is.na(skip_reason)) {
      capture_fit <- function(expression) {
        withCallingHandlers(
          tryCatch(expression, error = function(error) error),
          warning = function(warning) {
            warning_messages <<- unique(c(warning_messages, conditionMessage(warning)))
            invokeRestart("muffleWarning")
          }
        )
      }
      fit_independence <- capture_fit(nlme::lme(
        fixed = value ~ group,
        random = random_formula,
        data = work,
        method = "REML",
        control = fit_control,
        na.action = stats::na.fail
      ))
      fit_spatial <- capture_fit(nlme::lme(
        fixed = value ~ group,
        random = random_formula,
        correlation = corr_fun(form = correlation_formula, nugget = TRUE),
        data = work,
        method = "REML",
        control = fit_control,
        na.action = stats::na.fail
      ))
      if (inherits(fit_independence, "error") || inherits(fit_spatial, "error")) {
        errors <- c(
          if (inherits(fit_independence, "error")) paste0("independence: ", conditionMessage(fit_independence)),
          if (inherits(fit_spatial, "error")) paste0("spatial: ", conditionMessage(fit_spatial))
        )
        skip_reason <- paste(errors, collapse = " | ")
      }
    }
    if (!is.na(skip_reason)) {
      feature_rows[[feature_index]] <- data.frame(
        feature = feature, converged = FALSE, status = "failed",
        error_message = skip_reason, group_p_value_spatial = NA_real_,
        group_p_value_independence = NA_real_, range_parameter = NA_real_,
        nugget = NA_real_, AIC_independence = NA_real_, AIC_spatial = NA_real_,
        BIC_independence = NA_real_, BIC_spatial = NA_real_, delta_AIC = NA_real_,
        n = nrow(work), n_subjects = nlevels(work$subject),
        n_sections = nlevels(work$spatial_field), warning_message = paste(warning_messages, collapse = " | "),
        stringsAsFactors = FALSE
      )
      next
    }
    spatial_anova <- nlme::anova.lme(fit_spatial, type = "marginal")
    independence_anova <- nlme::anova.lme(fit_independence, type = "marginal")
    group_p_spatial <- spatial_anova["group", "p-value"]
    group_p_independence <- independence_anova["group", "p-value"]
    correlation_parameters <- stats::coef(fit_spatial$modelStruct$corStruct, unconstrained = FALSE)
    range_parameter <- unname(correlation_parameters["range"])
    nugget <- unname(correlation_parameters["nugget"])
    aic_independence <- stats::AIC(fit_independence)
    aic_spatial <- stats::AIC(fit_spatial)
    bic_independence <- stats::BIC(fit_independence)
    bic_spatial <- stats::BIC(fit_spatial)
    feature_rows[[feature_index]] <- data.frame(
      feature = feature, converged = TRUE, status = "fitted", error_message = NA_character_,
      group_p_value_spatial = group_p_spatial,
      group_p_value_independence = group_p_independence,
      range_parameter = range_parameter, nugget = nugget,
      AIC_independence = aic_independence, AIC_spatial = aic_spatial,
      BIC_independence = bic_independence, BIC_spatial = bic_spatial,
      delta_AIC = aic_spatial - aic_independence,
      n = nrow(work), n_subjects = nlevels(work$subject),
      n_sections = nlevels(work$spatial_field),
      warning_message = paste(warning_messages, collapse = " | "),
      stringsAsFactors = FALSE
    )
    if (isTRUE(keep_fits)) fits[[feature_index]] <- fit_spatial

    fixed_effects <- nlme::fixef(fit_spatial)
    covariance <- stats::vcov(fit_spatial)
    design <- stats::model.matrix(
      ~ group,
      data = data.frame(group = factor(groups, levels = groups))
    )
    colnames(design) <- names(fixed_effects)
    comparisons <- utils::combn(seq_along(groups), 2L)
    fixed_df <- fit_spatial$fixDF$X
    for (comparison in seq_len(ncol(comparisons))) {
      left_index <- comparisons[1, comparison]
      right_index <- comparisons[2, comparison]
      contrast <- design[right_index, ] - design[left_index, ]
      estimate <- unname(sum(contrast * fixed_effects))
      standard_error <- sqrt(unname(t(contrast) %*% covariance %*% contrast))
      active_coefficients <- which(abs(contrast) > sqrt(.Machine$double.eps))
      degrees_freedom <- if (length(active_coefficients)) min(fixed_df[active_coefficients]) else min(fixed_df)
      t_value <- estimate / standard_error
      p_value <- if (is.finite(t_value) && degrees_freedom > 0) {
        2 * stats::pt(abs(t_value), df = degrees_freedom, lower.tail = FALSE)
      } else NA_real_
      critical <- stats::qt((1 + confidence_level) / 2, df = degrees_freedom)
      contrast_rows[[contrast_index]] <- data.frame(
        feature = feature,
        group_a = groups[left_index], group_b = groups[right_index],
        estimate = estimate, standard_error = standard_error,
        degrees_freedom = degrees_freedom, t_value = t_value,
        confidence_lower = estimate - critical * standard_error,
        confidence_upper = estimate + critical * standard_error,
        p_value = p_value, stringsAsFactors = FALSE
      )
      contrast_index <- contrast_index + 1L
    }
  }
  feature_table <- do.call(rbind, feature_rows)
  feature_table$adj_group_p_value_spatial <- stats::p.adjust(
    feature_table$group_p_value_spatial, method = p_adjust_method
  )
  feature_table$adj_group_p_value_independence <- stats::p.adjust(
    feature_table$group_p_value_independence, method = p_adjust_method
  )
  feature_table$distance_unit <- as.character(distance_unit)
  feature_table$correlation_structure <- correlation_structure
  feature_table <- feature_table[order(feature_table$adj_group_p_value_spatial, na.last = TRUE), , drop = FALSE]
  rownames(feature_table) <- NULL
  contrast_table <- if (length(contrast_rows)) do.call(rbind, contrast_rows) else data.frame()
  if (nrow(contrast_table)) {
    contrast_table$fdr <- stats::ave(
      contrast_table$p_value,
      interaction(contrast_table$group_a, contrast_table$group_b, drop = TRUE),
      FUN = function(p) stats::p.adjust(p, method = p_adjust_method)
    )
    contrast_table <- contrast_table[order(contrast_table$fdr, na.last = TRUE), , drop = FALSE]
    rownames(contrast_table) <- NULL
  }
  if (isTRUE(keep_fits)) names(fits) <- fcols
  list(
    features = feature_table,
    contrasts = contrast_table,
    field_qc = field_qc,
    subject_counts = data.frame(group = names(subject_counts), n_subjects = unname(subject_counts)),
    fits = fits,
    settings = list(
      model_class = "Gaussian spatial linear mixed-effects model",
      group_column = group_column, subject_column = subject_column,
      section_column = section_column, coordinate_source = coordinate_source,
      x_resolution = x_resolution, y_resolution = y_resolution,
      distance_unit = as.character(distance_unit),
      correlation_structure = correlation_structure,
      random_effects = "nested random intercepts: subject/section",
      spatial_correlation_scope = "within subject/section",
      covariance_comparison = "AIC/BIC sensitivity summary; no regular chi-square LRT claimed"
    ),
    interpretation = paste(
      "This function fits Gaussian spatial LMMs, not non-Gaussian GLMMs.",
      "The spatial group test is an omnibus marginal F-test; contrasts provide direction.",
      "Range and nugget are covariance parameters in the declared distance unit, and model diagnostics remain required."
    )
  )
}

differential_region_analysis_glmm <- function(...) {
  warning(
    "differential_region_analysis_glmm() is a compatibility alias. The implemented model is a Gaussian spatial LMM, not a non-Gaussian GLMM; use differential_region_analysis_spatial_lmm().",
    call. = FALSE
  )
  differential_region_analysis_spatial_lmm(...)
}
