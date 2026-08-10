# Real-data adapters -------------------------------------------------------

validate_imzml_ibd_pair <- function(imzml_path, ibd_path = NULL) {
  imzml_path <- normalizePath(imzml_path, mustWork = TRUE)
  if (!identical(tolower(tools::file_ext(imzml_path)), "imzml")) {
    stop("imzml_path must identify an imzML file.", call. = FALSE)
  }
  if (is.null(ibd_path)) {
    stem <- substr(imzml_path, 1L, nchar(imzml_path) - nchar(tools::file_ext(imzml_path)))
    candidates <- c(paste0(stem, "ibd"), paste0(stem, "IBD"))
    ibd_path <- candidates[file.exists(candidates)][1]
  }
  if (length(ibd_path) != 1L || is.na(ibd_path) || !file.exists(ibd_path)) {
    stop("A readable companion ibd_path is required.", call. = FALSE)
  }
  ibd_path <- normalizePath(ibd_path, mustWork = TRUE)
  xml_head <- paste(readLines(imzml_path, n = 80L, warn = FALSE), collapse = " ")
  uuid <- sub('.*<fileContent[^>]*>.*', '', xml_head)
  uuid_hit <- regmatches(xml_head, regexpr("[0-9A-Fa-f]{8}-[0-9A-Fa-f]{4}-[0-9A-Fa-f]{4}-[0-9A-Fa-f]{4}-[0-9A-Fa-f]{12}", xml_head))
  ibd_uuid <- readBin(ibd_path, what = "raw", n = 16L)
  if (length(ibd_uuid) != 16L) stop("The ibd file is too short to contain its UUID.", call. = FALSE)
  uuid_match <- NA
  if (length(uuid_hit) == 1L && nzchar(uuid_hit)) {
    h <- gsub("-", "", tolower(uuid_hit))
    expected <- c(substr(h, 7, 8), substr(h, 5, 6), substr(h, 3, 4), substr(h, 1, 2),
                  substr(h, 11, 12), substr(h, 9, 10), substr(h, 15, 16), substr(h, 13, 14),
                  substring(h, seq(17, 31, 2), seq(18, 32, 2)))
    uuid_match <- identical(paste(sprintf("%02x", as.integer(ibd_uuid)), collapse = ""), paste(expected, collapse = ""))
    if (!uuid_match) stop("The imzML UUID does not match the ibd UUID.", call. = FALSE)
  }
  list(imzml_path = imzml_path, ibd_path = ibd_path, imzml_uuid = uuid_hit,
       ibd_uuid_checked = !is.na(uuid_match), uuid_match = uuid_match)
}

cardinal_aligned_matrix <- function(object, n_pixels) {
  values <- as.matrix(Cardinal::intensity(object))
  mz_values <- as.numeric(Cardinal::mz(object))
  if (nrow(values) == length(mz_values) && ncol(values) == n_pixels) {
    feature_by_pixel <- values
    orientation <- "features_x_pixels"
  } else if (ncol(values) == length(mz_values) && nrow(values) == n_pixels) {
    feature_by_pixel <- t(values)
    orientation <- "pixels_x_features"
  } else {
    stop("Cardinal aligned intensity dimensions do not match pixels and features.", call. = FALSE)
  }
  storage.mode(feature_by_pixel) <- "double"
  list(values = feature_by_pixel, mz = mz_values, source_orientation = orientation)
}

load_variable_mz_msi_features <- function(imzml_path, ibd_path = NULL,
                                           sample_id, section_id, ion_mode,
                                           ion_mode_source,
                                           processing = c("processed_peak_lists", "profile_diff"),
                                           alignment_ppm = 10,
                                           peak_pick_snr = 6,
                                           min_detection_fraction = 0,
                                           progress_callback = NULL) {
  processing <- match.arg(processing)
  pair <- validate_imzml_ibd_pair(imzml_path, ibd_path)
  text_args <- list(sample_id = sample_id, section_id = section_id,
                    ion_mode = ion_mode, ion_mode_source = ion_mode_source)
  bad <- vapply(text_args, function(x) length(x) != 1L || is.na(x) || !nzchar(trimws(x)), logical(1))
  if (any(bad)) stop(paste(names(text_args)[bad], collapse = ", "), " must be supplied explicitly.", call. = FALSE)
  if (!ion_mode %in% c("positive", "negative")) stop("ion_mode must be positive or negative.", call. = FALSE)
  if (!is.finite(alignment_ppm) || alignment_ppm <= 0) stop("alignment_ppm must be positive.", call. = FALSE)
  if (!is.finite(min_detection_fraction) || min_detection_fraction < 0 || min_detection_fraction > 1) {
    stop("min_detection_fraction must be between zero and one.", call. = FALSE)
  }
  if (!requireNamespace("Cardinal", quietly = TRUE)) stop("Cardinal is required.", call. = FALSE)
  if (!is.null(progress_callback)) progress_callback(0.02, "Reading imzML/ibd")
  cardinal_imzml <- pair$imzml_path
  imz_stem <- tools::file_path_sans_ext(basename(pair$imzml_path))
  ibd_stem <- tools::file_path_sans_ext(basename(pair$ibd_path))
  if (!identical(tolower(imz_stem), tolower(ibd_stem))) {
    pairing_dir <- tempfile("SpatialOmicsMSI-cardinal-pair-")
    dir.create(pairing_dir)
    cardinal_imzml <- file.path(pairing_dir, "paired.imzML")
    cardinal_ibd <- file.path(pairing_dir, "paired.ibd")
    if (!file.copy(pair$imzml_path, cardinal_imzml, overwrite = FALSE, copy.mode = TRUE) ||
        !file.symlink(pair$ibd_path, cardinal_ibd)) {
      stop("Could not create temporary read-only same-basename links for Cardinal.", call. = FALSE)
    }
    Sys.chmod(cardinal_imzml, mode = "0444")
  }
  object <- Cardinal::readMSIData(cardinal_imzml)
  coordinates <- as.data.frame(Cardinal::coord(object))
  required_columns(coordinates, c("x", "y"), "Cardinal coordinates")
  if (any(!is.finite(coordinates$x) | !is.finite(coordinates$y)) || anyDuplicated(coordinates[c("x", "y")])) {
    stop("Coordinates must be finite and unique.", call. = FALSE)
  }
  experiment <- Cardinal::experimentData(object)
  spectrum_type <- paste(as.character(experiment$spectrumType), collapse = "; ")
  if (nzchar(spectrum_type) && !grepl("MS1", spectrum_type, fixed = TRUE)) {
    stop("Only MS1 MSI input is supported; observed: ", spectrum_type, call. = FALSE)
  }
  raw_mz <- Cardinal::mz(object)
  variable_axis <- is_list_like_spectra(raw_mz) || inherits(raw_mz, "matter_list")
  if (!variable_axis) stop("This adapter requires per-spectrum variable m/z arrays.", call. = FALSE)
  metadata_representation <- if (isTRUE(Cardinal::isCentroided(object))) "centroided_cv" else "not_centroided_or_cv_absent"
  if (!is.null(progress_callback)) progress_callback(0.12, "Spectral processing")
  processed <- if (processing == "profile_diff") {
    picked <- Cardinal::peakPick(object, method = "diff", SNR = peak_pick_snr)
    Cardinal::peakAlign(picked, tolerance = alignment_ppm, units = "ppm")
  } else {
    Cardinal::peakAlign(object, tolerance = alignment_ppm, units = "ppm")
  }
  aligned <- cardinal_aligned_matrix(processed, nrow(coordinates))
  detection <- rowMeans(is.finite(aligned$values) & aligned$values > 0)
  keep <- detection >= min_detection_fraction
  if (!any(keep)) stop("No aligned features passed the detection-frequency filter.", call. = FALSE)
  intensity <- aligned$values[keep, , drop = FALSE]
  mz_values <- aligned$mz[keep]
  detection <- detection[keep]
  column_names <- stable_mz_column_names(mz_values)
  pixel_feature_matrix <- data.frame(
    pixel_id = seq_len(nrow(coordinates)), sample_id = sample_id, section_id = section_id,
    x = coordinates$x, y = coordinates$y, t(intensity), check.names = FALSE,
    stringsAsFactors = FALSE
  )
  names(pixel_feature_matrix)[-(1:5)] <- column_names
  feature_metadata <- data.frame(
    feature_id = sprintf("feature_%05d", seq_along(mz_values)), column_name = column_names,
    mz = mz_values, detection_fraction = detection, ion_mode = ion_mode,
    stringsAsFactors = FALSE
  )
  coordinates_out <- pixel_feature_matrix[c("pixel_id", "sample_id", "section_id", "x", "y")]
  raw_arrays <- as.list(Cardinal::intensity(object));raw_lengths <- lengths(raw_arrays)
  raw_tic <- vapply(raw_arrays,sum,numeric(1),na.rm=TRUE)
  raw_peak_count <- vapply(raw_arrays,function(v)sum(is.finite(v)&v>0),numeric(1))
  parameters <- list(
    input_type = "variable_mz_imzML", processing = processing,
    alignment_method = "Cardinal::peakAlign", alignment_ppm = alignment_ppm,
    peak_pick_method = if (processing == "profile_diff") "Cardinal::peakPick(method=diff)" else "not_applied",
    peak_pick_snr = if (processing == "profile_diff") peak_pick_snr else NA_real_,
    min_detection_fraction = min_detection_fraction, ion_mode = ion_mode,
    ion_mode_source = ion_mode_source, spectrum_representation_source = "imzML_CV_via_Cardinal"
  )
  parameters$temporary_same_basename_pair <- !identical(cardinal_imzml, pair$imzml_path)
  qc_summary <- data.frame(
    spectra_count = nrow(coordinates), aligned_feature_count = length(mz_values),
    raw_array_length_min = min(raw_lengths), raw_array_length_median = stats::median(raw_lengths),
    raw_array_length_max = max(raw_lengths), duplicate_coordinates = 0L,
    missing_coordinates = 0L, mz_min = min(mz_values), mz_max = max(mz_values),
    zero_fraction = mean(intensity == 0), stringsAsFactors = FALSE
  )
  provenance <- make_pipeline_manifest(c(pair$imzml_path, pair$ibd_path),
    input_type = parameters$input_type, ion_mode = ion_mode, ppm = alignment_ppm,
    parameters = parameters)
  provenance <- rbind(provenance, data.frame(
    record_type = "spectral_metadata",
    key = c("ms_level", "spectrum_representation", "mz_layout", "polarity_source", "imzml_ibd_uuid_match"),
    value = c(ifelse(nzchar(spectrum_type), spectrum_type, "not_reported"), metadata_representation,
              "per_spectrum_variable", ion_mode_source, as.character(pair$uuid_match)),
    size_bytes = NA_real_, md5 = NA_character_, stringsAsFactors = FALSE))
  if (!is.null(progress_callback)) progress_callback(1, "Variable-m/z processing complete")
  pixel_qc<-data.frame(coordinates_out,raw_tic=raw_tic,raw_peak_count=raw_peak_count)
  list(pixel_feature_matrix = pixel_feature_matrix, coordinates = coordinates_out,
       feature_metadata = feature_metadata, qc_summary = qc_summary,
       parameters = parameters, provenance = provenance,pixel_qc=pixel_qc)
}

read_metaspace_transform <- function(path) {
  if (!requireNamespace("jsonlite", quietly = TRUE)) stop("jsonlite is required.", call. = FALSE)
  value <- jsonlite::fromJSON(path, simplifyVector = FALSE)
  transform <- value$data$rawOpticalImage$transform
  if (is.null(transform)) transform <- value$transform
  matrix <- do.call(rbind, lapply(transform, as.numeric))
  if (!identical(dim(matrix), c(3L, 3L)) || any(!is.finite(matrix))) {
    stop("The platform transform must be a finite 3 x 3 matrix.", call. = FALSE)
  }
  list(transform = matrix,
       optical_url = value$data$rawOpticalImage$url %||% NA_character_,
       direction = "relative_MSI_coordinates_to_optical_pixels")
}

`%||%` <- function(x, y) if (is.null(x)) y else x

apply_metaspace_transform <- function(coordinates, transform,
                                      origin = c(min(coordinates$x), min(coordinates$y))) {
  required_columns(coordinates, c("x", "y"), "MSI coordinates")
  if (!identical(dim(transform), c(3L, 3L))) stop("transform must be 3 x 3.", call. = FALSE)
  relative <- cbind(coordinates$x - origin[1], coordinates$y - origin[2], 1)
  mapped <- relative %*% t(transform)
  data.frame(coordinates, relative_x = relative[, 1], relative_y = relative[, 2],
             optical_x = mapped[, 1] / mapped[, 3], optical_y = mapped[, 2] / mapped[, 3],
             check.names = FALSE)
}

registration_mask_diagnostics <- function(measurement, tissue) {
  if (!is.matrix(measurement) || !is.matrix(tissue) || !identical(dim(measurement), dim(tissue))) {
    stop("measurement and tissue must be logical matrices with identical dimensions.", call. = FALSE)
  }
  measurement <- measurement != FALSE; tissue <- tissue != FALSE
  intersection <- sum(measurement & tissue)
  boundary <- function(mask) {
    up <- rbind(FALSE, mask[-nrow(mask),,drop=FALSE]); down <- rbind(mask[-1,,drop=FALSE],FALSE)
    left <- cbind(FALSE,mask[,-ncol(mask),drop=FALSE]); right <- cbind(mask[,-1,drop=FALSE],FALSE)
    mask & !(up & down & left & right)
  }
  a <- which(boundary(measurement),arr.ind=TRUE); b <- which(boundary(tissue),arr.ind=TRUE)
  median_distance <- p95_distance <- NA_real_
  if(nrow(a)&&nrow(b)){
    nearest <- function(from,to,chunk=500L) unlist(lapply(split(seq_len(nrow(from)),ceiling(seq_len(nrow(from))/chunk)),function(ii){
      dx<-outer(from[ii,1],to[,1],"-");dy<-outer(from[ii,2],to[,2],"-");sqrt(apply(dx*dx+dy*dy,1,min))
    }),use.names=FALSE)
    symmetric<-c(nearest(a,b),nearest(b,a))
    median_distance<-stats::median(symmetric);p95_distance<-unname(stats::quantile(symmetric,.95))
  }
  c(overlap=intersection/max(1,sum(measurement)),
    dice=2*intersection/max(1,sum(measurement)+sum(tissue)),
    boundary_distance_median_px=median_distance,boundary_distance_p95_px=p95_distance)
}

project_msi_measurement_mask <- function(coordinates, transform, optical_dim) {
  origin<-c(min(coordinates$x),min(coordinates$y));gw<-max(coordinates$x)-origin[1]+1L;gh<-max(coordinates$y)-origin[2]+1L
  grid<-matrix(FALSE,gh,gw);grid[cbind(coordinates$y-origin[2]+1L,coordinates$x-origin[1]+1L)]<-TRUE
  inverse<-solve(transform);corners<-rbind(c(0,0,1),c(gw,0,1),c(0,gh,1),c(gw,gh,1))%*%t(transform)
  xr<-max(0,floor(min(corners[,1]))):min(optical_dim[2]-1L,ceiling(max(corners[,1])));yr<-max(0,floor(min(corners[,2]))):min(optical_dim[1]-1L,ceiling(max(corners[,2])))
  out<-matrix(FALSE,optical_dim[1],optical_dim[2]);for(y in yr){sx<-inverse[1,1]*(xr+.5)+inverse[1,2]*(y+.5)+inverse[1,3];sy<-inverse[2,1]*(xr+.5)+inverse[2,2]*(y+.5)+inverse[2,3];ix<-floor(sx);iy<-floor(sy);ok<-ix>=0&ix<gw&iy>=0&iy<gh;if(any(ok))out[y+1L,xr[ok]+1L]<-grid[cbind(iy[ok]+1L,ix[ok]+1L)]}
  out
}

flood_connected <- function(allowed, seed) {
  current<-matrix(FALSE,nrow(allowed),ncol(allowed));current[seed[1],seed[2]]<-TRUE
  repeat{grown<-allowed&(current|rbind(FALSE,current[-nrow(current),,drop=FALSE])|rbind(current[-1,,drop=FALSE],FALSE)|cbind(FALSE,current[,-ncol(current),drop=FALSE])|cbind(current[,-1,drop=FALSE],FALSE));if(sum(grown)==sum(current))break;current<-grown};current
}

matched_filled_optical_component <- function(tissue_raw, measurement) {
  rr<-which(measurement,arr.ind=TRUE);center<-c(stats::median(rr[,1]),stats::median(rr[,2]));cand<-which(tissue_raw,arr.ind=TRUE)
  seed<-cand[which.min((cand[,1]-center[1])^2+(cand[,2]-center[2])^2),]
  component<-flood_connected(tissue_raw,seed)
  outside_allowed<-!component;edge<-matrix(FALSE,nrow(component),ncol(component));edge[c(1,nrow(edge)),]<-outside_allowed[c(1,nrow(edge)),];edge[,c(1,ncol(edge))]<-outside_allowed[,c(1,ncol(edge))]
  seeds<-which(edge,arr.ind=TRUE);outside<-matrix(FALSE,nrow(edge),ncol(edge));outside[seeds]<-TRUE
  repeat{grown<-outside_allowed&(outside|rbind(FALSE,outside[-nrow(outside),,drop=FALSE])|rbind(outside[-1,,drop=FALSE],FALSE)|cbind(FALSE,outside[,-ncol(outside),drop=FALSE])|cbind(outside[,-1,drop=FALSE],FALSE));if(sum(grown)==sum(outside))break;outside<-grown}
  component|(!component&!outside)
}

build_msi_tissue_mask <- function(coordinates, tic, peak_count,
                                  method = c("kmeans_log_tic_peak_count", "score_quantile"),
                                  score_quantile = 0.50, seed = 20260808) {
  method<-match.arg(method);n<-nrow(coordinates)
  if(length(tic)!=n||length(peak_count)!=n)stop("tic and peak_count must have one value per coordinate.",call.=FALSE)
  metrics<-scale(cbind(log1p(tic),log1p(peak_count)));score<-rowSums(metrics)
  if(method=="kmeans_log_tic_peak_count"){set.seed(seed);km<-stats::kmeans(metrics,2,nstart=50);candidate<-km$cluster==which.max(tapply(score,km$cluster,mean));threshold<-NA_real_}else{threshold<-unname(stats::quantile(score,score_quantile));candidate<-score>=threshold}
  px<-data.frame(pixel_id=seq_len(n),x=coordinates$x,y=coordinates$y)
  cleaned<-cleanup_foreground_components(px,candidate,method="largest_component",connectivity="4")
  tissue<-cleaned$keep_rows
  xr<-range(coordinates$x);yr<-range(coordinates$y);grid<-matrix(FALSE,diff(yr)+1L,diff(xr)+1L);grid[cbind(coordinates$y-yr[1]+1L,coordinates$x-xr[1]+1L)]<-tissue
  allowed<-!grid;edge<-matrix(FALSE,nrow(grid),ncol(grid));edge[c(1,nrow(edge)),]<-allowed[c(1,nrow(edge)),];edge[,c(1,ncol(edge))]<-allowed[,c(1,ncol(edge))];outside<-matrix(FALSE,nrow(grid),ncol(grid));outside[which(edge,arr.ind=TRUE)]<-TRUE
  repeat{grown<-allowed&(outside|rbind(FALSE,outside[-nrow(outside),,drop=FALSE])|rbind(outside[-1,,drop=FALSE],FALSE)|cbind(FALSE,outside[,-ncol(outside),drop=FALSE])|cbind(outside[,-1,drop=FALSE],FALSE));if(sum(grown)==sum(outside))break;outside<-grown}
  holes<-!grid&!outside
  boundary_contact<-c(left=any(coordinates$x[tissue]==xr[1]),right=any(coordinates$x[tissue]==xr[2]),bottom=any(coordinates$y[tissue]==yr[1]),top=any(coordinates$y[tissue]==yr[2]))
  list(mask=data.frame(coordinates,tissue_status=ifelse(tissue,"tissue","unclassified/background"),candidate=candidate,tissue=tissue),
       diagnostics=data.frame(method=method,score_quantile=if(method=="score_quantile")score_quantile else NA_real_,score_threshold=threshold,seed=seed,total_scan_pixels=n,candidate_tissue_pixels=sum(candidate),final_tissue_pixels=sum(tissue),background_pixels=n-sum(tissue),components=cleaned$stats$foreground_components,internal_hole_pixels=sum(holes),touch_left=boundary_contact[1],touch_right=boundary_contact[2],touch_bottom=boundary_contact[3],touch_top=boundary_contact[4]))
}

register_metaspace_optical <- function(coordinates, transform_json, optical_path,
                                       attribution_path = NULL, tic = NULL,
                                       representative_ions = NULL,
                                       output_dir = NULL,
                                       tissue_threshold = 165 / 255) {
  if (!file.exists(optical_path)) stop("Optical brightfield image is missing.", call. = FALSE)
  platform <- read_metaspace_transform(transform_json)
  origin <- c(min(coordinates$x), min(coordinates$y))
  registered <- apply_metaspace_transform(coordinates, platform$transform, origin)
  attribution <- if (!is.null(attribution_path) && file.exists(attribution_path) &&
                     requireNamespace("jsonlite", quietly = TRUE)) jsonlite::fromJSON(attribution_path) else NULL
  diagnostics <- data.frame(n_pixels = nrow(registered),
    optical_x_min = min(registered$optical_x), optical_x_max = max(registered$optical_x),
    optical_y_min = min(registered$optical_y), optical_y_max = max(registered$optical_y),
    y_axis_inversion = FALSE, overlap = NA_real_, dice = NA_real_,
    boundary_distance_median_px = NA_real_, boundary_distance_p95_px = NA_real_,
    stringsAsFactors = FALSE)
  output_files <- character()
  intermediates <- NULL
  if (requireNamespace("jpeg", quietly = TRUE)) {
    image <- jpeg::readJPEG(optical_path)
    gray <- if (length(dim(image)) == 3L) apply(image, c(1, 2), mean) else image
    h <- nrow(gray); w <- ncol(gray)
    measurement <- project_msi_measurement_mask(coordinates,platform$transform,c(h,w))
    tissue_raw <- gray < tissue_threshold
    pad<-50L;rows<-max(1,min(which(measurement,arr.ind=TRUE)[,1])-pad):min(h,max(which(measurement,arr.ind=TRUE)[,1])+pad);cols<-max(1,min(which(measurement,arr.ind=TRUE)[,2])-pad):min(w,max(which(measurement,arr.ind=TRUE)[,2])+pad)
    tissue<-matrix(FALSE,h,w);tissue[rows,cols]<-matched_filled_optical_component(tissue_raw[rows,cols],measurement[rows,cols])
    intermediates<-list(raster_dimensions=c(height=h,width=w),measurement_mask=measurement,
      optical_binary_mask=tissue_raw,matched_component_mask=tissue,
      threshold=tissue_threshold,morphology="4-connected component nearest transformed-mask center; binary hole filling")
    mask_diagnostics <- registration_mask_diagnostics(measurement, tissue)
    diagnostics[names(mask_diagnostics)] <- as.list(mask_diagnostics)
    if (!is.null(output_dir)) {
      dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
      overlay <- function(values, name, title) {
        path <- file.path(output_dir, name); grDevices::png(path, 1400, 1000, res=130)
        graphics::plot.new(); graphics::plot.window(c(0,w),c(h,0),asp=1)
        graphics::rasterImage(image,0,h,w,0)
        if (is.null(values)) graphics::points(registered$optical_x,registered$optical_y,pch=15,cex=.22,col=grDevices::adjustcolor("#00D7FF",.4)) else {
          limits <- stats::quantile(values[is.finite(values)],c(.01,.99),na.rm=TRUE)
          z <- pmin(pmax(values,limits[1]),limits[2]); pal <- grDevices::hcl.colors(256,"Inferno")
          idx <- 1L+pmin(255L,pmax(0L,as.integer((z-limits[1])/max(diff(limits),.Machine$double.eps)*255)))
          graphics::points(registered$optical_x,registered$optical_y,pch=15,cex=.22,col=grDevices::adjustcolor(pal[idx],.65))
        }
        graphics::box();graphics::title(title);grDevices::dev.off();path
      }
      output_files <- c(measurement_mask_overlay=overlay(NULL,"measurement_mask_overlay.png","Platform-registered MSI measurement mask"))
      if (!is.null(tic)) output_files <- c(output_files,tic_overlay=overlay(tic,"tic_overlay.png","TIC on optical brightfield"))
      if (!is.null(representative_ions)) for (name in names(representative_ions)) {
        output_files <- c(output_files,stats::setNames(overlay(representative_ions[[name]],paste0("ion_",make.names(name),"_overlay.png"),paste("Ion",name,"on optical brightfield")),name))
      }
    }
  }
  list(registered_coordinates = registered, transform = platform$transform,
       transform_direction = platform$direction, coordinate_origin = origin,
       coordinate_convention = "x right, y down; no implicit y-axis inversion",
       optical_path = normalizePath(optical_path), attribution = attribution,
       diagnostics = diagnostics, output_files = output_files,
       registration_intermediates = intermediates)
}

decode_mzml_binary <- function(block) {
  if (!requireNamespace("base64enc", quietly = TRUE)) stop("base64enc is required.", call. = FALSE)
  binary <- sub(".*<binary[^>]*>", "", block)
  binary <- sub("</binary>.*", "", binary)
  raw <- base64enc::base64decode(gsub("[[:space:]]", "", binary))
  if (grepl("MS:1000574", block, fixed = TRUE)) raw <- memDecompress(raw, type = "gzip")
  size <- if (grepl("MS:1000523", block, fixed = TRUE)) 8L else 4L
  readBin(raw, what = "double", n = length(raw) %/% size, size = size, endian = "little")
}

mzml_cv_value <- function(block, accession) {
  hit <- regmatches(block, regexpr(paste0('<cvParam[^>]*accession="', accession, '"[^>]*/?>'), block))
  if (!length(hit) || !nzchar(hit)) return(NA_real_)
  value <- sub('.*value="([^"]*)".*', "\\1", hit)
  suppressWarnings(as.numeric(value))
}

read_mzml_fragment_spectra <- function(mzml_path, precursor_target = NULL,
                                       precursor_ppm_tolerance = 10) {
  mzml_path <- normalizePath(mzml_path, mustWork = TRUE)
  con <- file(mzml_path, open = "rt"); on.exit(close(con), add = TRUE)
  scans <- list(); peaks <- list(); block <- character(); inside <- FALSE; index <- 0L
  repeat {
    lines <- readLines(con, n = 200L, warn = FALSE)
    if (!length(lines)) break
    for (line in lines) {
      if (grepl("<spectrum[ >]", line)) { inside <- TRUE; block <- line } else if (inside) block <- c(block, line)
      if (inside && grepl("</spectrum>", line, fixed = TRUE)) {
        text <- paste(block, collapse = "")
        inside <- FALSE; block <- character()
        if (!grepl('accession="MS:1000511"[^>]*value="2"', text)) next
        precursor <- mzml_cv_value(text, "MS:1000744")
        if (!is.null(precursor_target) && (is.na(precursor) ||
            abs(precursor - precursor_target) / precursor_target * 1e6 > precursor_ppm_tolerance)) next
        index <- index + 1L
        id <- sub('.*<spectrum[^>]*id="([^"]+)".*', "\\1", text)
        arrays <- regmatches(text, gregexpr("<binaryDataArray[[:space:][:print:]]*?</binaryDataArray>", text, perl = TRUE))[[1]]
        mz_block <- arrays[grepl("MS:1000514", arrays, fixed = TRUE)][1]
        int_block <- arrays[grepl("MS:1000515", arrays, fixed = TRUE)][1]
        mz <- decode_mzml_binary(mz_block); intensity <- decode_mzml_binary(int_block)
        n <- min(length(mz), length(intensity)); mz <- mz[seq_len(n)]; intensity <- intensity[seq_len(n)]
        bp <- if (n) max(intensity, na.rm = TRUE) else NA_real_
        scans[[index]] <- data.frame(scan_id = id, scan_index = index, ms_level = 2L,
          precursor_mz = precursor, retention_time = mzml_cv_value(text, "MS:1000016"),
          collision_energy = mzml_cv_value(text, "MS:1000045"),
          isolation_target_mz = mzml_cv_value(text, "MS:1000827"),
          isolation_lower_offset = mzml_cv_value(text, "MS:1000828"),
          isolation_upper_offset = mzml_cv_value(text, "MS:1000829"),
          total_ion_intensity = sum(intensity, na.rm = TRUE), base_peak_intensity = bp,
          spectrum_representation = if (grepl("MS:1000127", text, fixed = TRUE)) "centroided" else if (grepl("MS:1000128", text, fixed = TRUE)) "profile" else "not_reported",
          stringsAsFactors = FALSE)
        peaks[[index]] <- data.frame(scan_id = id, scan_index = index,
          fragment_mz = mz, fragment_intensity = intensity,
          peak_source = if (grepl("MS:1000127", text, fixed = TRUE)) "vendor_centroid" else "mzML_array",
          extraction_method = "none", extraction_parameters = "none",
          software_version = paste0("SpatialOmicsMSI ", utils::packageVersion("SpatialOmicsMSI")),
          stringsAsFactors = FALSE)
      }
    }
  }
  list(precursor_scan_metadata = if (length(scans)) do.call(rbind, scans) else data.frame(),
       fragment_peak_table = if (length(peaks)) do.call(rbind, peaks) else data.frame(),
       parameters = list(precursor_target = precursor_target,
                         precursor_ppm_tolerance = precursor_ppm_tolerance,
                         chemical_identity_inferred = FALSE))
}
