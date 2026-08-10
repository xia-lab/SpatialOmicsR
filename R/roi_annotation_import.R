# External ROI annotation import and LMD geometry QC ----------------------

import_qupath_geojson <- function(path, section_id = NULL) {
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    stop("Package 'jsonlite' is required to import QuPath GeoJSON.", call. = FALSE)
  }
  if (length(path) != 1L || !file.exists(path)) {
    stop("QuPath GeoJSON file does not exist: ", path, call. = FALSE)
  }
  document <- jsonlite::fromJSON(path, simplifyVector = FALSE)
  document_type <- if (is.list(document) && !is.null(document$type)) {
    as.character(document$type)
  } else NA_character_
  features <- if (identical(document_type, "FeatureCollection")) {
    document$features
  } else if (identical(document_type, "Feature")) {
    list(document)
  } else if (is.list(document) && length(document) &&
             all(vapply(document, function(item) identical(as.character(item$type), "Feature"), logical(1)))) {
    document
  } else {
    stop("Expected a GeoJSON FeatureCollection, Feature, or array of Features.", call. = FALSE)
  }
  empty_result <- function() data.frame(
    roi_id = character(), section_id = character(), polygon_part_id = integer(),
    ring_id = integer(), ring_role = character(), vertex_order = integer(),
    x = numeric(), y = numeric(), classification = character(),
    annotation_name = character(), source_feature_index = integer(),
    source_format = character(), coordinate_unit = character(),
    coordinate_origin = character(), y_axis_direction = character(),
    stringsAsFactors = FALSE
  )
  if (!length(features)) return(empty_result())

  property_value <- function(properties, name) {
    if (is.null(properties) || is.null(properties[[name]])) return(NA_character_)
    value <- properties[[name]]
    if (is.list(value)) {
      if (!is.null(value$name)) value <- value$name
      else if (!is.null(value$names)) value <- paste(unlist(value$names), collapse = ":")
      else return(NA_character_)
    }
    value <- as.character(unlist(value))[1L]
    if (is.na(value) || !nzchar(value)) NA_character_ else value
  }
  feature_metadata <- lapply(seq_along(features), function(i) {
    feature <- features[[i]]
    if (!identical(as.character(feature$type), "Feature") || is.null(feature$geometry)) {
      stop("Every item must be a GeoJSON Feature with geometry.", call. = FALSE)
    }
    geometry_type <- as.character(feature$geometry$type)
    if (!geometry_type %in% c("Polygon", "MultiPolygon")) {
      stop("Unsupported QuPath geometry type '", geometry_type,
           "'; only Polygon and MultiPolygon are accepted.", call. = FALSE)
    }
    classification <- property_value(feature$properties, "classification")
    annotation_name <- property_value(feature$properties, "name")
    feature_id <- if (!is.null(feature$id)) as.character(feature$id)[1L] else NA_character_
    base_id <- if (!is.na(annotation_name)) annotation_name else if (!is.na(feature_id) && nzchar(feature_id)) {
      feature_id
    } else if (!is.na(classification)) classification else "roi"
    list(feature = feature, geometry_type = geometry_type,
         classification = classification, annotation_name = annotation_name,
         base_id = base_id, feature_index = i)
  })
  roi_ids <- make.unique(vapply(feature_metadata, `[[`, character(1), "base_id"), sep = "_")
  rows <- list()
  for (i in seq_along(feature_metadata)) {
    metadata <- feature_metadata[[i]]
    geometry <- metadata$feature$geometry
    polygons <- if (metadata$geometry_type == "Polygon") {
      list(geometry$coordinates)
    } else geometry$coordinates
    if (!length(polygons)) stop("Feature ", i, " has no polygon coordinates.", call. = FALSE)
    for (part_index in seq_along(polygons)) {
      rings <- polygons[[part_index]]
      if (!length(rings)) stop("Feature ", i, " contains an empty polygon part.", call. = FALSE)
      for (ring_index in seq_along(rings)) {
        points <- rings[[ring_index]]
        coordinates <- do.call(rbind, lapply(points, function(point) {
          value <- suppressWarnings(as.numeric(unlist(point)))
          if (length(value) < 2L) c(NA_real_, NA_real_) else value[1:2]
        }))
        if (nrow(coordinates) < 4L || any(!is.finite(coordinates)) ||
            !isTRUE(all.equal(coordinates[1L, ], coordinates[nrow(coordinates), ], tolerance = 0))) {
          stop("Every GeoJSON polygon ring must be finite, closed, and contain at least four coordinates.", call. = FALSE)
        }
        if (nrow(unique(data.frame(x = coordinates[, 1L], y = coordinates[, 2L]))) < 3L) {
          stop("A GeoJSON polygon ring has fewer than three distinct vertices.", call. = FALSE)
        }
        rows[[length(rows) + 1L]] <- data.frame(
          roi_id = roi_ids[i],
          section_id = if (is.null(section_id)) NA_character_ else as.character(section_id),
          polygon_part_id = part_index,
          ring_id = ring_index,
          ring_role = if (ring_index == 1L) "outer" else "hole",
          vertex_order = seq_len(nrow(coordinates)),
          x = coordinates[, 1L], y = coordinates[, 2L],
          classification = metadata$classification,
          annotation_name = metadata$annotation_name,
          source_feature_index = metadata$feature_index,
          source_format = "qupath_geojson",
          coordinate_unit = "full_resolution_pixel",
          coordinate_origin = "top_left",
          y_axis_direction = "down",
          stringsAsFactors = FALSE
        )
      }
    }
  }
  result <- if (length(rows)) do.call(rbind, rows) else empty_result()
  rownames(result) <- NULL
  attr(result, "coordinate_warning") <- paste(
    "QuPath GeoJSON coordinates use full-resolution image pixels with a",
    "top-left origin. Register and verify axes before transfer to MSI or LMD."
  )
  result
}

validate_lmd_shapes <- function(
    polygon_table,
    max_area,
    min_spacing,
    min_collectable_area = 0,
    roi_column = "roi_id",
    section_column = NULL,
    part_column = "polygon_part_id",
    ring_column = "ring_id",
    role_column = "ring_role",
    x_column = "x",
    y_column = "y",
    coordinate_unit) {
  required_columns(
    polygon_table,
    c(roi_column, section_column, part_column, ring_column, role_column, x_column, y_column),
    "LMD polygon table"
  )
  limits <- c(max_area, min_spacing, min_collectable_area)
  if (any(!is.finite(limits)) || max_area <= 0 || min_spacing < 0 || min_collectable_area < 0) {
    stop("max_area must be positive; min_spacing and min_collectable_area must be non-negative.", call. = FALSE)
  }
  if (length(coordinate_unit) != 1L || is.na(coordinate_unit) || !nzchar(coordinate_unit)) {
    stop("coordinate_unit must be supplied explicitly, for example 'um'.", call. = FALSE)
  }
  role <- as.character(polygon_table[[role_column]])
  if (any(is.na(role) | !role %in% c("outer", "hole"))) {
    stop("ring roles must contain only 'outer' or 'hole'.", call. = FALSE)
  }
  x <- suppressWarnings(as.numeric(polygon_table[[x_column]]))
  y <- suppressWarnings(as.numeric(polygon_table[[y_column]]))
  if (any(!is.finite(x)) || any(!is.finite(y))) stop("Polygon coordinates must be finite.", call. = FALSE)
  work <- polygon_table
  work$.__x__ <- x; work$.__y__ <- y
  work$.__section__ <- if (is.null(section_column)) ".__single__" else as.character(work[[section_column]])
  shape_key <- interaction(
    as.character(work[[roi_column]]), work$.__section__, as.character(work[[part_column]]),
    drop = TRUE, lex.order = TRUE
  )
  shape_groups <- split(seq_len(nrow(work)), shape_key)

  ring_coordinates <- function(indices) {
    ring_groups <- split(indices, as.character(work[[ring_column]][indices]))
    lapply(ring_groups, function(idx) {
      points <- data.frame(x = work$.__x__[idx], y = work$.__y__[idx])
      if (nrow(points) > 1L &&
          !(points$x[1L] == points$x[nrow(points)] && points$y[1L] == points$y[nrow(points)])) {
        points <- rbind(points, points[1L, , drop = FALSE])
      }
      list(points = points, role = unique(role[idx]))
    })
  }
  signed_area <- function(points) {
    if (nrow(points) < 4L) return(NA_real_)
    sum(points$x[-nrow(points)] * points$y[-1L] -
          points$x[-1L] * points$y[-nrow(points)]) / 2
  }
  orientation <- function(ax, ay, bx, by, cx, cy) (bx - ax) * (cy - ay) - (by - ay) * (cx - ax)
  segments_intersect <- function(a, b, c, d) {
    o1 <- orientation(a[1], a[2], b[1], b[2], c[1], c[2])
    o2 <- orientation(a[1], a[2], b[1], b[2], d[1], d[2])
    o3 <- orientation(c[1], c[2], d[1], d[2], a[1], a[2])
    o4 <- orientation(c[1], c[2], d[1], d[2], b[1], b[2])
    (o1 * o2 < 0) && (o3 * o4 < 0)
  }
  self_intersects <- function(points) {
    n_segments <- nrow(points) - 1L
    if (n_segments < 3L) return(TRUE)
    for (i in seq_len(n_segments)) for (j in seq_len(n_segments)) {
      if (j <= i + 1L || (i == 1L && j == n_segments)) next
      if (segments_intersect(as.numeric(points[i, ]), as.numeric(points[i + 1L, ]),
                             as.numeric(points[j, ]), as.numeric(points[j + 1L, ]))) return(TRUE)
    }
    FALSE
  }
  point_segment_distance <- function(point, a, b) {
    delta <- b - a
    denom <- sum(delta^2)
    if (denom == 0) return(sqrt(sum((point - a)^2)))
    position <- max(0, min(1, sum((point - a) * delta) / denom))
    sqrt(sum((point - (a + position * delta))^2))
  }
  ring_distance <- function(first, second) {
    best <- Inf
    for (i in seq_len(nrow(first) - 1L)) for (j in seq_len(nrow(second) - 1L)) {
      a <- as.numeric(first[i, ]); b <- as.numeric(first[i + 1L, ])
      c <- as.numeric(second[j, ]); d <- as.numeric(second[j + 1L, ])
      if (segments_intersect(a, b, c, d)) return(0)
      best <- min(best, point_segment_distance(a, c, d), point_segment_distance(b, c, d),
                  point_segment_distance(c, a, b), point_segment_distance(d, a, b))
    }
    if (point_in_polygon(first$x[1L], first$y[1L], second$x, second$y) ||
        point_in_polygon(second$x[1L], second$y[1L], first$x, first$y)) 0 else best
  }

  shapes <- lapply(shape_groups, function(indices) {
    rings <- ring_coordinates(indices)
    role_unique <- vapply(rings, function(item) length(unique(item$role)) == 1L, logical(1))
    roles <- vapply(rings, function(item) item$role[1L], character(1))
    role_valid <- all(role_unique) && sum(roles == "outer") == 1L && all(roles %in% c("outer", "hole"))
    ring_valid <- vapply(rings, function(item) {
      nrow(unique(item$points)) >= 3L && !self_intersects(item$points)
    }, logical(1))
    areas <- vapply(rings, function(item) abs(signed_area(item$points)), numeric(1))
    outer <- if (sum(roles == "outer") == 1L) rings[[which(roles == "outer")]]$points else NULL
    holes_inside <- role_valid && all(vapply(rings[roles == "hole"], function(item) {
      all(point_in_polygon(item$points$x, item$points$y, outer$x, outer$y, boundary_inside = TRUE))
    }, logical(1)))
    net_area <- if (role_valid && holes_inside && all(is.finite(areas))) {
      sum(areas[roles == "outer"]) - sum(areas[roles == "hole"])
    } else NA_real_
    list(
      roi_id = as.character(work[[roi_column]][indices[1L]]),
      section_id = work$.__section__[indices[1L]],
      part_id = as.character(work[[part_column]][indices[1L]]),
      outer = outer,
      row = data.frame(
        roi_id = as.character(work[[roi_column]][indices[1L]]),
        section_id = if (is.null(section_column)) NA_character_ else work$.__section__[indices[1L]],
        polygon_part_id = as.character(work[[part_column]][indices[1L]]),
        area = net_area, n_rings = length(rings),
        n_vertices = sum(vapply(rings, function(item) nrow(item$points) - 1L, integer(1))),
        valid_ring_roles = role_valid, holes_inside_outer = holes_inside,
        self_intersection = any(!ring_valid),
        stringsAsFactors = FALSE
      )
    )
  })
  result <- do.call(rbind, lapply(shapes, `[[`, "row"))
  minimum_spacing <- rep(Inf, length(shapes))
  if (length(shapes) > 1L) for (i in seq_len(length(shapes) - 1L)) {
    for (j in seq.int(i + 1L, length(shapes))) {
      if (shapes[[i]]$section_id != shapes[[j]]$section_id ||
          is.null(shapes[[i]]$outer) || is.null(shapes[[j]]$outer)) next
      distance <- ring_distance(shapes[[i]]$outer, shapes[[j]]$outer)
      minimum_spacing[c(i, j)] <- pmin(minimum_spacing[c(i, j)], distance)
    }
  }
  minimum_spacing[!is.finite(minimum_spacing)] <- NA_real_
  result$minimum_spacing_to_other_shape <- minimum_spacing
  result$area_above_maximum <- is.finite(result$area) & result$area > max_area
  result$area_below_collectable_minimum <- is.finite(result$area) & result$area < min_collectable_area
  result$spacing_below_minimum <- is.finite(minimum_spacing) & minimum_spacing < min_spacing
  result$qc_pass <- result$valid_ring_roles & result$holes_inside_outer & !result$self_intersection &
    is.finite(result$area) & result$area > 0 & !result$area_above_maximum &
    !result$area_below_collectable_minimum & !result$spacing_below_minimum
  result$qc_reason <- vapply(seq_len(nrow(result)), function(i) {
    reason <- character()
    if (!result$valid_ring_roles[i]) reason <- c(reason, "invalid_ring_roles")
    if (!result$holes_inside_outer[i]) reason <- c(reason, "hole_outside_outer")
    if (result$self_intersection[i]) reason <- c(reason, "self_intersection")
    if (!is.finite(result$area[i]) || result$area[i] <= 0) reason <- c(reason, "invalid_area")
    if (result$area_above_maximum[i]) reason <- c(reason, "area_above_maximum")
    if (result$area_below_collectable_minimum[i]) reason <- c(reason, "area_below_collectable_minimum")
    if (result$spacing_below_minimum[i]) reason <- c(reason, "spacing_below_minimum")
    if (length(reason)) paste(reason, collapse = ";") else "pass"
  }, character(1))
  result$coordinate_unit <- coordinate_unit
  result$max_area <- max_area
  result$min_spacing <- min_spacing
  result$min_collectable_area <- min_collectable_area
  attr(result, "scope") <- paste(
    "Read-only geometry QC. No polygon is split, buffered, simplified, or",
    "otherwise altered for instrument cutting."
  )
  rownames(result) <- NULL
  result
}
