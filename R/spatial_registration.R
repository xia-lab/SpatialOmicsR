# Histology-to-MSI coordinate registration --------------------------------

fit_histology_msi_registration <- function(control_points,
                                           histology_columns = c("histology_x", "histology_y"),
                                           msi_columns = c("msi_x", "msi_y"),
                                           section_column = NULL) {
  required_columns(
    control_points,
    c(histology_columns, msi_columns, section_column),
    "Histology/MSI control-point table"
  )
  if (length(histology_columns) != 2L || length(msi_columns) != 2L) {
    stop("histology_columns and msi_columns must each contain x and y column names.", call. = FALSE)
  }

  section_values <- if (is.null(section_column)) {
    rep(".__single__", nrow(control_points))
  } else {
    as.character(control_points[[section_column]])
  }
  section_values[is.na(section_values) | section_values == ""] <- ".__missing__"
  section_groups <- split(seq_len(nrow(control_points)), section_values)

  fits <- lapply(names(section_groups), function(section_id) {
    idx <- section_groups[[section_id]]
    points <- control_points[idx, c(histology_columns, msi_columns), drop = FALSE]
    points[] <- lapply(points, function(value) suppressWarnings(as.numeric(value)))
    finite <- stats::complete.cases(points) & apply(points, 1, function(value) all(is.finite(value)))
    points <- points[finite, , drop = FALSE]
    if (nrow(points) < 3L) {
      stop("At least three finite control points are required for section '", section_id, "'.", call. = FALSE)
    }

    design <- cbind(intercept = 1, x = points[[histology_columns[1]]], y = points[[histology_columns[2]]])
    if (qr(design)$rank < 3L) {
      stop("Histology control points must not be collinear in section '", section_id, "'.", call. = FALSE)
    }
    target <- as.matrix(points[, msi_columns, drop = FALSE])
    coefficients <- qr.solve(design, target)
    fitted <- design %*% coefficients
    residuals <- target - fitted
    distances <- sqrt(rowSums(residuals^2))

    list(
      section_id = section_id,
      coefficients = coefficients,
      n_control_points = nrow(points),
      rmse = sqrt(mean(rowSums(residuals^2))),
      max_error = max(distances),
      residuals = residuals
    )
  })
  names(fits) <- names(section_groups)

  structure(
    list(
      fits = fits,
      histology_columns = histology_columns,
      msi_columns = msi_columns,
      section_column = section_column
    ),
    class = "histology_msi_registration"
  )
}

transform_histology_coordinates <- function(coordinates,
                                            registration,
                                            x_column = "x",
                                            y_column = "y",
                                            output_columns = c("x", "y"),
                                            section_column = registration$section_column) {
  if (!inherits(registration, "histology_msi_registration")) {
    stop("registration must come from fit_histology_msi_registration().", call. = FALSE)
  }
  required_columns(coordinates, c(x_column, y_column, section_column), "Histology coordinate table")
  if (length(output_columns) != 2L) {
    stop("output_columns must contain the MSI x and y column names.", call. = FALSE)
  }

  out <- coordinates
  source_x <- suppressWarnings(as.numeric(out[[x_column]]))
  source_y <- suppressWarnings(as.numeric(out[[y_column]]))
  section_values <- if (is.null(section_column)) {
    rep(".__single__", nrow(out))
  } else {
    as.character(out[[section_column]])
  }
  section_values[is.na(section_values) | section_values == ""] <- ".__missing__"
  unknown <- setdiff(unique(section_values), names(registration$fits))
  if (length(unknown)) {
    stop("No registration is available for section(s): ", paste(unknown, collapse = ", "), call. = FALSE)
  }

  transformed_x <- rep(NA_real_, nrow(out))
  transformed_y <- rep(NA_real_, nrow(out))
  for (section_id in unique(section_values)) {
    idx <- which(section_values == section_id & is.finite(source_x) & is.finite(source_y))
    if (!length(idx)) next
    design <- cbind(intercept = 1, x = source_x[idx], y = source_y[idx])
    transformed <- design %*% registration$fits[[section_id]]$coefficients
    transformed_x[idx] <- transformed[, 1]
    transformed_y[idx] <- transformed[, 2]
  }
  out[[output_columns[1]]] <- transformed_x
  out[[output_columns[2]]] <- transformed_y
  out
}

registration_diagnostics <- function(registration) {
  if (!inherits(registration, "histology_msi_registration")) {
    stop("registration must come from fit_histology_msi_registration().", call. = FALSE)
  }
  do.call(rbind, lapply(registration$fits, function(fit) {
    data.frame(
      section_id = if (identical(fit$section_id, ".__single__")) NA_character_ else fit$section_id,
      n_control_points = fit$n_control_points,
      rmse = fit$rmse,
      max_error = fit$max_error,
      stringsAsFactors = FALSE
    )
  }))
}
