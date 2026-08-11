options(shiny.maxRequestSize = 5 * 1024^3)
suppressPackageStartupMessages({
  library(shiny)
  library(SpatialOmicsMSI)
  library(ggplot2)
  library(DT)
})

project_root <- Sys.getenv("SPATIALOMICS_PROJECT_ROOT", unset = "")
if (!nzchar(project_root)) {
  candidate <- normalizePath(file.path(getwd(), "../../.."), mustWork = FALSE)
  project_root <- if (dir.exists(file.path(candidate, "data_raw"))) candidate else getwd()
}
session_root <- Sys.getenv("SPATIALOMICS_SESSION_ROOT", unset = tempdir())

`%or%` <- function(x, y) if (is.null(x) || length(x) == 0L || is.na(x) || !nzchar(x)) y else x
present_path <- function(x) length(x) == 1L && !is.na(x) && nzchar(trimws(x))
require_columns <- function(data, columns, label) {
  missing <- setdiff(columns, names(data))
  if (length(missing)) stop(label, " is missing required column(s): ", paste(missing, collapse = ", "), call. = FALSE)
}
normalize_feature_names <- function(columns) {
  output <- columns; numeric <- suppressWarnings(is.finite(as.numeric(columns)))
  output[numeric] <- paste0("mz_", columns[numeric]); gsub("-", "_", output, fixed = TRUE)
}

example_datasets <- list(
  brain01 = list(
    label = "Brain01 registration example",
    msi_input_type = "imzml", imzml_path = file.path(project_root,
      "data_raw/mouse_brain_he_msi/metaspace_brain01/Brain01_Bregma-1-46_centroid.imzML"),
    ibd_path = file.path(project_root,
      "data_raw/mouse_brain_he_msi/metaspace_brain01/Brain01_Bregma-1-46_centroid.ibd"),
    optical_path = file.path(project_root,
      "data_raw/mouse_brain_he_msi/metaspace_brain01/optical_brightfield.jpg"),
    transform_path = file.path(project_root,
      "data_raw/mouse_brain_he_msi/metaspace_brain01/optical_transform_api.json"),
    attribution_path = file.path(project_root,
      "data_raw/mouse_brain_he_msi/metaspace_brain01/attribution_license.json"),
    sample_id = "Brain01", section_id = "Brain01_Bregma-1-46", subject_id = "Brain01",
    ion_mode = "positive", ion_mode_source = "METASPACE_dataset_metadata",
    spectral_processing = "processed_peak_lists", alignment_ppm = 10,
    peak_pick_snr = 6, min_detection_fraction = 0, make_tissue_mask = FALSE,
    lcms_relationship = "none",
    source_note = "METASPACE 2016-09-22_11h16m17s; CC BY 4.0; LAVIGNE Régis."
  ),
  omix = list(
    label = "OMIX016317 full-brain MSI example",
    msi_input_type = "imzml", imzml_path = file.path(project_root,
      "data_raw/full_tissue_mouse_brain/OMIX016317/OMIX016317-02.imzML"),
    ibd_path = file.path(project_root,
      "data_raw/full_tissue_mouse_brain/OMIX016317/OMIX016317-01.ibd"),
    sample_id = "mbrain1_neg100", section_id = "mbrain1_neg100", subject_id = "mbrain1",
    ion_mode = "negative", ion_mode_source = "OMIX_record_and_imzML_audit",
    spectral_processing = "profile_diff", alignment_ppm = 10,
    peak_pick_snr = 6, min_detection_fraction = .10, make_tissue_mask = TRUE,
    lcms_relationship = "none",
    source_note = "OMIX016317 mbrain1_neg100; independent complete-brain MSI technical example."
  ),
  lcms = list(
    label = "MSV000090179 LC-MS/MS example", msi_input_type = "none",
    lcms_path = file.path(project_root,
      "data_raw/mouse_brain_lcms/msv000090179/pos_mouse_female_brain_12w_1.mzML"),
    sample_id = "mouse_brain_12w_rep1", section_id = "not_applicable",
    subject_id = "mouse_brain_12w_rep1", ion_mode = "positive",
    ion_mode_source = "MassIVE_sample_metadata", spectral_processing = "processed_peak_lists",
    alignment_ppm = 10, peak_pick_snr = 6, min_detection_fraction = 0,
    make_tissue_mask = FALSE, lcms_relationship = "external_reference",
    source_note = "MassIVE MSV000090179; CC0; female mouse brain, 12 weeks, replicate 1."
  ),
  msiflow = list(
    label = "MSIflow MSI + LC-MS/MS study example", msi_input_type = "csv",
    csv_path = file.path(project_root,
      "results/real_data/UPEC_12/pixel_feature_matrix.csv"),
    labels_path = file.path(project_root,
      "results/real_data/UPEC_12/spatial_domain_labels.csv"),
    lcms_path = file.path(project_root,
      "data_raw/msiflow_lcms/massive_mzml/POS/Lumos_ISAS_Ecoli_Mouse_Urinarybladder_Sample_2_rep1_Pos_1.mzML"),
    sample_id = "UPEC_12", section_id = "UPEC_12", subject_id = "UPEC_12",
    ion_mode = "positive", ion_mode_source = "MSIflow study files",
    spectral_processing = "processed_peak_lists", alignment_ppm = 10,
    peak_pick_snr = 6, min_detection_fraction = 0, make_tissue_mask = FALSE,
    lcms_relationship = "study_matched",
    source_note = paste(
      "Previously validated MSIflow UPEC_12 pixel-feature output plus an LC-MS/MS file from the associated study.",
      "Treat them as study-level corroborating evidence unless specimen-level pairing is independently verified.")
  )
)

vicari_preset <- function(sample_id, section_id) {
  st_prefix <- file.path(project_root, "data_raw/vicari_sma/st_processed", section_id)
  list(
    label = paste0("Vicari 2024 ", sample_id, " MSI + Visium + H&E (registration required)"),
    msi_input_type = "csv",
    csv_path = file.path(project_root, "data_raw/vicari_sma/msi_raw",
      paste0("ITO_FMP10_", sample_id, ".csv")),
    sample_id = sample_id, section_id = section_id, subject_id = sample_id,
    ion_mode = "positive", ion_mode_source = "Vicari et al. Figshare metadata (FMP-10)",
    spectral_processing = "processed_peak_lists", alignment_ppm = 10,
    peak_pick_snr = 6, min_detection_fraction = 0, make_tissue_mask = TRUE,
    lcms_relationship = "none",
    st_expression_path = paste0(st_prefix, "_RNA_filtered_feature_bc_matrix.h5"),
    st_positions_path = paste0(st_prefix, "_RNA_tissue_positions_list.csv"),
    histology_path = paste0(st_prefix, "_RNA_tissue_hires_image.png"),
    st_region_path = paste0(st_prefix, "_RNA_region.csv"),
    st_lesion_path = paste0(st_prefix, "_RNA_lesion.csv"),
    st_loupe_path = paste0(st_prefix, "_RNA_RegionLoupe.csv"),
    source_note = paste(
      "Vicari et al. 2024 matched tissue section; MSI SCiLS export, Visium expression,",
      "H&E and manual annotations downloaded from Figshare."
    )
  )
}
example_datasets$vicari_mpd1 <- vicari_preset("mPD1", "V11L12-038_D1")
example_datasets$vicari_mpd3 <- vicari_preset("mPD3", "V11L12-109_A1")
example_datasets$vicari_mpd4 <- vicari_preset("mPD4", "V11L12-109_B1")
example_datasets$kasarla_kidney <- list(
  label = "Kasarla 2025 kidney MSI + LMD-LC-MS/MS study bundle",
  msi_input_type = "csv",
  csv_path = file.path(project_root,
    "data_raw/kasarla2025/processed/kasarla_kidney_pixel_feature_matrix.rds"),
  sample_id = "Kasarla_kidney_Dglu_washed",
  section_id = "2023-12-19_13h02m39s",
  subject_id = "Kasarla_Dglu_mouse",
  ion_mode = "negative",
  ion_mode_source = "METASPACE dataset 2023-12-19_13h02m39s",
  spectral_processing = "processed_peak_lists", alignment_ppm = 3,
  peak_pick_snr = 6, min_detection_fraction = 0, make_tissue_mask = FALSE,
  lcms_relationship = "study_matched",
  source_note = paste(
    "Kasarla et al. 2025 kidney study bundle. Pixel-level MSI matrix reconstructed from",
    "236 HMDB annotation ion images at 10% FDR for METASPACE 2023-12-19_13h02m39s;",
    "the downstream validation uses the published cortex, medulla and renal-pelvis",
    "MALDI/LMD measurements and the 12 study-matched LMD-LC-MS/MS file manifest."
  )
)

preset_spec <- function(key) {
  if (!key %in% names(example_datasets)) stop("Unknown example dataset.", call. = FALSE)
  example_datasets[[key]]
}

example_workflow <- function(spec) {
  has_msi <- !identical(spec$msi_input_type, "none")
  has_lcms <- present_path(spec$lcms_path)
  if (has_msi && has_lcms) "matched" else if (has_msi) "spatial" else "lcms"
}
spatial_example_datasets <- Filter(function(spec) !identical(spec$msi_input_type, "none"), example_datasets)
example_dataset_choices <- list(
  "End-to-end spatial MSI studies" = setNames(names(spatial_example_datasets),
    vapply(spatial_example_datasets, `[[`, character(1), "label")),
  "Module-only LC-MS/MS evidence" = setNames("lcms", example_datasets$lcms$label)
)

cardinal_pair_path <- function(pair) {
  same_stem <- identical(tolower(tools::file_path_sans_ext(basename(pair$imzml_path))),
                           tolower(tools::file_path_sans_ext(basename(pair$ibd_path))))
  if (same_stem) return(list(imzml_path = pair$imzml_path, temporary = FALSE))
  directory <- tempfile("SpatialOmicsMSI-cardinal-pair-"); dir.create(directory)
  imzml <- file.path(directory, "paired.imzML"); ibd <- file.path(directory, "paired.ibd")
  if (!file.copy(pair$imzml_path, imzml, overwrite = FALSE, copy.mode = TRUE) || !file.symlink(pair$ibd_path, ibd)) {
    stop("Could not create a temporary read-only same-basename pair for Cardinal.", call. = FALSE)
  }
  Sys.chmod(imzml, "0444")
  list(imzml_path = imzml, temporary = TRUE)
}

inspect_imzml_input <- function(imzml_path, ibd_path) {
  pair <- validate_imzml_ibd_pair(imzml_path, ibd_path)
  if (!requireNamespace("Cardinal", quietly = TRUE)) stop("Cardinal is required to inspect imzML metadata.", call. = FALSE)
  readable <- cardinal_pair_path(pair)
  object <- Cardinal::readMSIData(readable$imzml_path)
  coordinates <- as.data.frame(Cardinal::coord(object))
  require_columns(coordinates, c("x", "y"), "imzML coordinates")
  if (!nrow(coordinates) || any(!is.finite(coordinates$x) | !is.finite(coordinates$y))) {
    stop("imzML coordinates must be finite and non-empty.", call. = FALSE)
  }
  if (anyDuplicated(coordinates[c("x", "y")])) stop("imzML contains duplicate x/y coordinates.", call. = FALSE)
  experiment <- Cardinal::experimentData(object)
  spectrum_type <- paste(as.character(experiment$spectrumType), collapse = "; ")
  if (nzchar(spectrum_type) && !grepl("MS1", spectrum_type, fixed = TRUE)) {
    stop("The MSI workflow requires MS1 spectra; metadata reports: ", spectrum_type, call. = FALSE)
  }
  raw_mz <- Cardinal::mz(object)
  variable_axis <- is.list(raw_mz) || inherits(raw_mz, "matter_list")
  cv_text <- paste(readLines(pair$imzml_path, n = 1000L, warn = FALSE), collapse = " ")
  positive_cv <- grepl("MS:1000130", cv_text, fixed = TRUE)
  negative_cv <- grepl("MS:1000129", cv_text, fixed = TRUE)
  polarity_cv <- if (xor(positive_cv, negative_cv)) if (positive_cv) "positive" else "negative" else "not reported or ambiguous"
  list(pair = pair, coordinates = coordinates,
       summary = data.frame(
         item = c("Spectra/pixels", "x range", "y range", "MS level", "Spectrum representation",
                  "m/z layout", "Polarity interpretation"),
         value = c(nrow(coordinates), paste(range(coordinates$x), collapse = " to "),
                   paste(range(coordinates$y), collapse = " to "), spectrum_type %or% "not reported",
                   if (isTRUE(Cardinal::isCentroided(object))) "centroided (CV)" else "profile or CV not reported",
                   if (variable_axis) "per-spectrum variable" else "shared axis", polarity_cv),
         stringsAsFactors = FALSE),
       mz_layout = if (variable_axis) "variable" else "shared", polarity_cv = polarity_cv)
}

load_shared_axis_pair <- function(imzml_path, ibd_path, sample_id, section_id,
                                  ion_mode, ion_mode_source) {
  pair <- validate_imzml_ibd_pair(imzml_path, ibd_path)
  readable <- cardinal_pair_path(pair)
  output <- load_centroided_msi_features(readable$imzml_path, sample_id, section_id, ion_mode, ion_mode_source)
  output$parameters$temporary_same_basename_pair <- readable$temporary
  output
}

read_pixel_feature_csv <- function(path, sample_id, section_id, subject_id) {
  data <- if (identical(tolower(tools::file_ext(path)), "rds")) {
    readRDS(path)
  } else {
    utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  }
  if (!is.data.frame(data)) data <- as.data.frame(data, check.names = FALSE)
  require_columns(data, c("x", "y"), "Processed pixel-by-feature table")
  if (!nrow(data) || any(!is.finite(data$x) | !is.finite(data$y)) || anyDuplicated(data[c("x", "y")])) {
    stop("Pixel table x/y coordinates must be finite, non-empty and unique.", call. = FALSE)
  }
  metadata <- c("pixel_id", "sample_id", "section_id", "subject_id", "x", "y", "domain_id", "domain_label")
  candidates <- setdiff(names(data), metadata)
  numeric_feature <- vapply(data[candidates], is.numeric, logical(1))
  features <- candidates[numeric_feature]
  if (!length(features)) stop("Pixel table must contain at least one numeric feature column in addition to x/y.", call. = FALSE)
  normalized <- normalize_feature_names(features)
  normalized <- make.unique(normalized, sep = "__")
  values <- as.data.frame(data[features], check.names = FALSE)
  names(values) <- normalized
  pixel_id <- if ("pixel_id" %in% names(data)) data$pixel_id else seq_len(nrow(data))
  if (anyNA(pixel_id) || anyDuplicated(pixel_id)) stop("Pixel-table pixel_id must be non-missing and unique.", call. = FALSE)
  pixel <- data.frame(pixel_id = pixel_id, sample_id = sample_id, section_id = section_id,
                      x = data$x, y = data$y, values, check.names = FALSE)
  mz <- suppressWarnings(as.numeric(sub("^mz_", "", sub("__.*$", "", normalized))))
  feature_metadata <- data.frame(feature_id = sprintf("feature_%05d", seq_along(features)),
    column_name = normalized, mz = mz, source_column = features, ion_mode = NA_character_)
  coordinates <- data.frame(pixel_id = pixel_id, sample_id = sample_id, section_id = section_id,
                            subject_id = subject_id, x = data$x, y = data$y)
  matrix <- as.matrix(values); storage.mode(matrix) <- "double"
  pixel_qc <- data.frame(coordinates, raw_tic = rowSums(matrix, na.rm = TRUE),
                         raw_peak_count = rowSums(is.finite(matrix) & matrix > 0))
  list(pixel_feature_matrix = pixel, coordinates = coordinates,
       feature_metadata = feature_metadata,
       qc_summary = data.frame(spectra_count = nrow(data), feature_count = length(features),
         x_min = min(data$x), x_max = max(data$x), y_min = min(data$y), y_max = max(data$y),
         missing_intensity = sum(!is.finite(matrix)), zero_fraction = mean(matrix == 0, na.rm = TRUE)),
       parameters = list(input_type = "processed_pixel_feature_csv"),
       provenance = make_pipeline_manifest(path, input_type = "processed_pixel_feature_csv",
         parameters = list(sample_id = sample_id, section_id = section_id, subject_id = subject_id)),
       pixel_qc = pixel_qc)
}

validate_label_file <- function(path) {
  labels <- utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  require_columns(labels, c("x", "y"), "ROI/domain labels")
  label_col <- intersect(c("domain_id", "roi_id", "label", "domain_label", "roi_label"), names(labels))[1]
  if (is.na(label_col)) stop("ROI/domain CSV needs one label column: domain_id, roi_id, label, domain_label or roi_label.", call. = FALSE)
  key_cols <- intersect(c("sample_id", "section_id"), names(labels))
  if (anyDuplicated(labels[c(key_cols, "x", "y")])) stop("ROI/domain labels are not unique by their coordinate key.", call. = FALSE)
  list(data = labels, label_column = label_col)
}

validate_analysis_spec <- function(spec, inspect = TRUE) {
  errors <- character(); notes <- character(); metadata <- NULL
  has_msi <- spec$msi_input_type %in% c("imzml", "csv")
  has_lcms <- present_path(spec$lcms_path)
  lcms_relationship <- spec$lcms_relationship %or% "none"
  if (!has_msi) errors <- c(errors, "Provide MSI data. LC-MS/MS is optional downstream evidence and cannot start an analysis by itself.")
  if (has_lcms && identical(lcms_relationship, "none")) {
    errors <- c(errors, "Declare whether the LC-MS/MS evidence is specimen-matched, study-matched, or an external reference.")
  }
  if (!has_lcms && !identical(lcms_relationship, "none")) {
    notes <- c(notes, "An LC-MS/MS relationship was declared, but no LC-MS/MS mzML file was supplied.")
  }
  if (has_msi) {
    for (field in c("sample_id", "section_id", "subject_id", "ion_mode", "ion_mode_source")) {
      if (!present_path(spec[[field]])) errors <- c(errors, paste(field, "is required for MSI analysis."))
    }
    if (!spec$ion_mode %in% c("positive", "negative")) errors <- c(errors, "Ion mode must be positive or negative.")
  }
  required_files <- character()
  if (identical(spec$msi_input_type, "imzml")) required_files <- c(required_files, spec$imzml_path, spec$ibd_path)
  if (identical(spec$msi_input_type, "csv")) required_files <- c(required_files, spec$csv_path)
  optional_files <- unlist(spec[c("optical_path", "transform_path", "labels_path", "lcms_path", "attribution_path")], use.names = FALSE)
  supplied <- c(required_files, optional_files[vapply(optional_files, present_path, logical(1))])
  missing <- supplied[!file.exists(supplied)]
  if (length(missing)) errors <- c(errors, paste(length(missing), "supplied file(s) do not exist; see Technical details."))
  if (!length(errors) && identical(spec$msi_input_type, "imzml")) {
    result <- tryCatch(if (inspect) inspect_imzml_input(spec$imzml_path, spec$ibd_path)
                       else list(pair = validate_imzml_ibd_pair(spec$imzml_path, spec$ibd_path), mz_layout = NA_character_),
                       error = function(e) e)
    if (inherits(result, "error")) errors <- c(errors, conditionMessage(result)) else {
      metadata <- result; notes <- c(notes, "imzML/ibd pairing and coordinate metadata passed.")
      if (inspect && result$polarity_cv %in% c("positive", "negative") && !identical(result$polarity_cv, spec$ion_mode)) {
        errors <- c(errors, paste0("Explicit ion mode ('", spec$ion_mode, "') conflicts with imzML CV polarity ('", result$polarity_cv, "')."))
      }
      if (inspect && !result$polarity_cv %in% c("positive", "negative")) notes <- c(notes, "imzML CV polarity was not definitive; the explicit ion-mode source is retained in provenance.")
    }
  }
  if (!length(errors) && identical(spec$msi_input_type, "csv")) {
    result <- tryCatch(read_pixel_feature_csv(spec$csv_path, spec$sample_id, spec$section_id, spec$subject_id), error = function(e) e)
    if (inherits(result, "error")) errors <- c(errors, conditionMessage(result)) else {
      metadata <- list(coordinates = result$coordinates, mz_layout = "table",
        summary = data.frame(item = c("Pixels", "Features", "x range", "y range", "Polarity interpretation"),
          value = c(nrow(result$coordinates), nrow(result$feature_metadata),
            paste(range(result$coordinates$x), collapse = " to "), paste(range(result$coordinates$y), collapse = " to "),
            "supplied explicitly by the user")))
      notes <- c(notes, "Processed pixel-by-feature table structure passed.")
    }
  }
  if (present_path(spec$optical_path) != present_path(spec$transform_path)) {
    notes <- c(notes, "Registration is disabled until both an optical JPEG and a transform JSON are supplied.")
  }
  if (!length(errors) && present_path(spec$optical_path) && present_path(spec$transform_path)) {
    if (!tolower(tools::file_ext(spec$optical_path)) %in% c("jpg", "jpeg")) errors <- c(errors, "Optical input must currently be JPEG.")
    optical_ok <- tryCatch({ jpeg::readJPEG(spec$optical_path); TRUE }, error = function(e) e)
    transform_ok <- tryCatch({ read_metaspace_transform(spec$transform_path); TRUE }, error = function(e) e)
    if (inherits(optical_ok, "error")) errors <- c(errors, paste("Optical image:", conditionMessage(optical_ok)))
    if (inherits(transform_ok, "error")) errors <- c(errors, paste("Registration transform:", conditionMessage(transform_ok)))
  }
  labels <- NULL
  if (!length(errors) && present_path(spec$labels_path)) {
    labels <- tryCatch(validate_label_file(spec$labels_path), error = function(e) e)
    if (inherits(labels, "error")) errors <- c(errors, conditionMessage(labels))
    else if (!is.null(metadata$coordinates)) {
      relevant <- labels$data
      if ("sample_id" %in% names(relevant)) relevant <- relevant[as.character(relevant$sample_id) == spec$sample_id, , drop = FALSE]
      if ("section_id" %in% names(relevant)) relevant <- relevant[as.character(relevant$section_id) == spec$section_id, , drop = FALSE]
      pixel_key <- paste(metadata$coordinates$x, metadata$coordinates$y, sep = "\r")
      label_key <- paste(relevant$x, relevant$y, sep = "\r")
      unknown <- sum(!label_key %in% pixel_key)
      if (unknown) errors <- c(errors, paste(unknown, "ROI/domain coordinate(s) do not occur in the MSI input."))
      else notes <- c(notes, paste(nrow(relevant), "ROI/domain coordinate(s) matched the MSI coordinate system."))
    }
  }
  if (!length(errors) && has_lcms) {
    if (tolower(tools::file_ext(spec$lcms_path)) != "mzml") errors <- c(errors, "LC-MS/MS input must be mzML.")
    else {
      connection <- file(spec$lcms_path, open = "rb"); on.exit(close(connection), add = TRUE)
      header <- readChar(connection, nchars = min(file.info(spec$lcms_path)$size, 2e6), useBytes = TRUE)
      if (!grepl("<mzML", header, fixed = TRUE)) errors <- c(errors, "LC-MS/MS input does not contain an mzML document header.")
    }
  }
  capabilities <- list(msi = has_msi && !length(errors), processing = has_msi && !length(errors),
    registration = has_msi && present_path(spec$optical_path) && present_path(spec$transform_path) && !length(errors),
    labels = has_msi && present_path(spec$labels_path) && !length(errors),
    domains = has_msi && !length(errors), moran = has_msi && !length(errors),
    lcms = has_lcms && !length(errors), download = !length(errors))
  list(valid = !length(errors), errors = unique(errors), notes = unique(notes), metadata = metadata,
       labels = if (inherits(labels, "error")) NULL else labels, capabilities = capabilities,
       files = supplied)
}

resolve_label_mapping <- function(labels_info, coordinates, sample_id, section_id) {
  labels <- labels_info$data
  if ("sample_id" %in% names(labels)) labels <- labels[as.character(labels$sample_id) == sample_id, , drop = FALSE]
  if ("section_id" %in% names(labels)) labels <- labels[as.character(labels$section_id) == section_id, , drop = FALSE]
  label_key <- paste(labels$x, labels$y, sep = "\r")
  pixel_key <- paste(coordinates$x, coordinates$y, sep = "\r")
  index <- match(pixel_key, label_key)
  raw <- rep("-1", length(index)); raw[!is.na(index)] <- as.character(labels[[labels_info$label_column]][index[!is.na(index)]])
  data.frame(coordinates, domain_id = raw,
    domain_label = ifelse(raw == "-1", "unclassified/background", paste0("user-supplied domain/ROI ", raw)),
    label_source = "user-supplied ROI/domain CSV", stringsAsFactors = FALSE)
}

workflow_card <- function(title, description, button_id, button_label) {
  column(4, tags$div(class = "well", style = "min-height:210px",
    h3(title), tags$p(description), actionButton(button_id, button_label, class = "btn-primary")))
}

ui <- navbarPage(
  id = "workflow_tabs", title = "SpatialOmicsMSI",
  header = tagList(
    tags$style(HTML(paste(
      ".workflow-status{margin:8px 12px}.workflow-status .label{font-size:90%;margin-right:5px}.well h3{margin-top:4px}",
      ".registration-workspace{display:flex;align-items:flex-start;gap:18px;margin:12px 0 18px}",
      ".registration-image-card,.registration-control-card{background:#fff;border:1px solid #ddd;border-radius:6px;padding:14px}",
      ".registration-image-card{flex:1 1 58%;min-width:0}.registration-control-card{flex:1 1 42%;min-width:300px}",
      ".responsive-image-output{height:min(56vh,520px)!important;min-height:280px}",
      ".responsive-image-output img{width:100%!important;height:100%!important;object-fit:contain;object-position:center}",
      ".shiny-notification-panel{z-index:1050!important}.shiny-notification{opacity:1!important}",
      ".server-busy-card{display:none;position:fixed;right:22px;bottom:22px;z-index:1045;align-items:center;gap:11px;max-width:340px;padding:12px 16px;background:#fff;border:1px solid #bce8f1;border-left:5px solid #337ab7;border-radius:7px;box-shadow:0 3px 16px rgba(0,0,0,.18);color:#245269}",
      ".shiny-busy .server-busy-card{display:flex}.server-busy-spinner{width:24px;height:24px;flex:0 0 24px;border:3px solid #d9edf7;border-top-color:#337ab7;border-radius:50%;animation:server-busy-spin .75s linear infinite}",
      ".server-busy-copy{line-height:1.25}.server-busy-copy strong{display:block}.server-busy-copy small{color:#666}",
      ".result-card{border:1px solid #ddd;border-radius:7px;padding:14px;margin:12px 0;background:#fff;overflow:hidden}",
      ".matched-layout{display:grid;grid-template-columns:minmax(0,3fr) minmax(320px,2fr);gap:16px;align-items:start}",
      ".matched-image{height:min(62vh,650px)!important;min-height:320px;overflow:hidden}.matched-image img{width:100%!important;height:100%!important;object-fit:contain!important;display:block}",
      ".matched-tables{min-width:0}.matched-tables .datatables{overflow-x:auto}",
      "@keyframes server-busy-spin{to{transform:rotate(360deg)}}",
      "@media(max-width:900px){.registration-workspace{display:block}.registration-control-card{min-width:0;margin-top:12px}.responsive-image-output{height:42vh!important}.matched-layout{display:block}.matched-tables{margin-top:12px}}"
    ))),
    tags$div(class = "server-busy-card", role = "status", `aria-live` = "polite",
      tags$div(class = "server-busy-spinner"),
      tags$div(class = "server-busy-copy",
        tags$strong("Working…"), tags$small("The page will update automatically when this step finishes."))),
    tags$div(class = "workflow-status", uiOutput("workflow_status"))
  ),
  tabPanel("Home", value = "home",
    h2("Spatial MSI analysis"),
    tags$p("Run the spatial MSI workflow from input through tissue-level results. LC-MS/MS is used only as downstream chemical-identification evidence for MSI candidates."),
    fluidRow(
      workflow_card("Spatial MSI", "Process spatial spectra, inspect QC, detect domains and niches, define regions, run spatial statistics, validate candidate identities with optional LC-MS/MS evidence, and map results back to tissue.", "choose_spatial", "Start analysis")
    ),
    hr(), h3("Current session"), uiOutput("workflow_next_step"), uiOutput("workflow_progress_home")
  ),
  tabPanel("Data setup", value = "data",
    tags$div(class = "alert alert-info",
      "Provide your MSI input here. Module-specific example buttons are located beside the functions they actually validate."),
    hr(),
    tags$div(class = "well",
      h3("Choose an example dataset"),
      tags$p("End-to-end studies contain pixel-level MSI and start at Step 1. Module-only evidence datasets open the specific downstream function they can validate."),
      selectInput("example_dataset", "Example dataset", choices = example_dataset_choices),
      actionButton("load_example_dataset", "Load selected example", class = "btn-success")),
    hr(),
    radioButtons("data_origin", "Input location", inline = TRUE,
      choices = c("Server files" = "server", "Browser upload" = "upload"), selected = "server"),
    selectInput("msi_input_type", "Primary MSI input", choices = c("Paired imzML + ibd" = "imzml",
      "Processed pixel × feature CSV" = "csv")),
    conditionalPanel("input.data_origin == 'upload' && input.msi_input_type == 'imzml'",
      fileInput("upload_imzml", "imzML", accept = c(".imzML", ".imzml")), fileInput("upload_ibd", "ibd", accept = c(".ibd", ".IBD"))),
    conditionalPanel("input.data_origin == 'upload' && input.msi_input_type == 'csv'", fileInput("upload_csv", "Pixel × feature CSV", accept = ".csv")),
    conditionalPanel("input.data_origin == 'upload'",
      fileInput("upload_optical", "Optional optical/H&E/brightfield JPEG", accept = c(".jpg", ".jpeg")),
      fileInput("upload_transform", "Optional registration transform JSON", accept = ".json"),
      fileInput("upload_labels", "Optional ROI/domain CSV", accept = ".csv"),
      fileInput("upload_lcms", "Optional LC-MS/MS evidence mzML", accept = c(".mzML", ".mzml"))),
    tags$details(tags$summary("Technical details — server paths"),
      conditionalPanel("input.data_origin == 'server' && input.msi_input_type == 'imzml'",
        textInput("server_imzml", "imzML path"), textInput("server_ibd", "ibd path")),
      conditionalPanel("input.data_origin == 'server' && input.msi_input_type == 'csv'",
        textInput("server_csv", "Pixel × feature table path (CSV or internal RDS cache)")),
      conditionalPanel("input.data_origin == 'server'",
        textInput("server_optical", "Optional optical JPEG path"), textInput("server_transform", "Optional transform JSON path"),
        textInput("server_labels", "Optional ROI/domain CSV path"),
        textInput("server_lcms", "Optional LC-MS/MS evidence mzML path"),
        textInput("server_attribution", "Optional attribution metadata path"))),
    selectInput("lcms_relationship", "LC-MS/MS evidence relationship",
      c("No LC-MS/MS evidence supplied" = "none",
        "Same specimen / matched aliquot" = "specimen_matched",
        "Same study, specimen pairing not established" = "study_matched",
        "External reference sample or library" = "external_reference")),
    fluidRow(column(4, textInput("sample_id", "Sample ID")), column(4, textInput("section_id", "Section ID")),
      column(4, textInput("subject_id", "Biological subject ID"))),
    fluidRow(column(4, selectInput("ion_mode", "Ion mode", c("Positive" = "positive", "Negative" = "negative"))),
      column(8, textInput("ion_mode_source", "Ion-mode metadata source", placeholder = "e.g. acquisition record; inferred from adducts"))),
    tags$details(tags$summary("Processing parameters"),
      selectInput("spectral_processing", "Spectral processing", c("Align stored processed peak lists" = "processed_peak_lists", "Profile peakPick(diff) then align" = "profile_diff")),
      numericInput("alignment_ppm", "Alignment tolerance (ppm)", 10, min = .01),
      numericInput("peak_pick_snr", "Profile peak-picking SNR", 6, min = 0),
      numericInput("min_detection_fraction", "Minimum detection fraction", 0, min = 0, max = 1, step = .01),
      checkboxInput("make_tissue_mask", "Build transparent TIC/peak-count tissue mask", FALSE)),
    actionButton("validate_input", "Validate input", class = "btn-primary"),
    uiOutput("validation_card"), uiOutput("continue_processing"),
    tags$details(tags$summary("Technical details — validated files and metadata"), DTOutput("technical_files"), DTOutput("input_metadata"))
  ),
  navbarMenu("Analysis workflow",
  tabPanel("1. Input & provenance", value = "spatial_provenance", uiOutput("spatial_gate"), uiOutput("module_overview"), DTOutput("provenance_table"), textOutput("source_note"), uiOutput("next_from_provenance")),
  tabPanel("2. Processing & QC", value = "spatial_processing", uiOutput("spatial_gate"),
    uiOutput("active_dataset_card"),
    checkboxInput("tic_normalize", "TIC normalization", TRUE),
    checkboxInput("log1p_transform", "log10(x + 1) transformation", TRUE),
    uiOutput("processing_action"), uiOutput("processing_gate"),
    tags$details(tags$summary("Example datasets for MSI import, processing and QC"),
      actionButton("use_brain01_qc", "View Brain01 QC example", class = "btn-info"),
      actionButton("use_omix_qc", "View OMIX016317 QC example", class = "btn-info"),
      actionButton("use_msiflow_qc", "Run MSIflow UPEC_12 QC example", class = "btn-info")),
    uiOutput("qc_example_results"),
    DTOutput("qc_table"), DTOutput("feature_preview"), uiOutput("tissue_mask_results"),
    hr(), h3("Spatial feature screening (optional)"),
    tags$p("Compute Moran's I once for every m/z feature. The cached result can screen features before domain detection and later check spatial coherence of statistical hits."),
    fluidRow(
      column(4, selectInput("neighbor_method", "Adjacency", c("Queen" = "queen", "Rook" = "rook"))),
      column(4, numericInput("permutations", "Two-sided permutations", 499, min = 99, step = 100)),
      column(4, br(), uiOutput("moran_action"))),
    uiOutput("moran_workload_warning"),
    uiOutput("moran_gate"), uiOutput("moran_screening_ui"), uiOutput("next_from_processing")),
  tabPanel("3. Registration (optional)", value = "spatial_registration", uiOutput("spatial_gate"),
    uiOutput("registration_workspace"), uiOutput("registration_gate"), uiOutput("registration_results"),
    tags$div(class = "well", strong("Example dataset: "),
      "Brain01 validates optical-image registration with a supplied transform. ",
      actionButton("use_brain01_registration", "View Brain01 registration example", class = "btn-info")),
    uiOutput("registration_example_results"), uiOutput("next_from_registration")),
  tabPanel("4. Optional pre-analysis ROI", value = "spatial_roi", uiOutput("spatial_gate"),
    h3("Choose how much tissue enters domain detection"),
    tags$div(class = "alert alert-info",
      "Recommended default: use the complete tissue, detect molecular domains, then choose domains as the downstream ROI. Use a prior ROI only when the target area is already known independently."),
    radioButtons("roi_path", "Input area for domain detection",
      choices = c("Complete tissue (recommended; choose domains afterward)" = "whole_tissue",
        "Known anatomical area from registered H&E" = "histology",
        "Known area drawn directly on the MSI coordinate map" = "msi"),
      selected = "whole_tissue"),
    conditionalPanel("input.roi_path == 'whole_tissue'",
      tags$div(class = "alert alert-success",
        "No drawing is required. Continue to Domains & niches; the complete tissue mask will be analyzed.")),
    conditionalPanel("input.roi_path == 'histology'",
      tags$div(class = "well",
        h4("Path A — anatomical / histology ROI"),
        tags$p("Import polygon vertices drawn in H&E pixel coordinates. A validated H&E→MSI control-point registration is required before the polygon can be transferred."),
        fileInput("histology_roi_polygon", "H&E polygon CSV (roi_id, x, y; optional vertex_order)", accept = ".csv"),
        uiOutput("histology_roi_action"), uiOutput("histology_roi_gate"))),
    conditionalPanel("input.roi_path == 'msi'",
      tags$div(class = "well",
        h4("Direct MSI-coordinate ROI (optional)"),
        tags$p("This selects pixels by their x/y positions on the MSI map. Use it when no reliable H&E registration is available but a region is already known from the MSI image or acquisition coordinates."),
        fluidRow(
          column(4, textInput("roi_name", "ROI name", "roi_01")),
          column(8, br(), actionButton("add_brushed_roi", "Add brushed rectangle", class = "btn-primary"),
            actionButton("use_full_tissue_roi", "Use complete tissue"),
            actionButton("clear_rois", "Clear ROI"))
        ),
        fileInput("roi_table_upload", "Optional MSI coordinate table (roi_id, x_min, x_max, y_min, y_max)", accept = ".csv"))),
    conditionalPanel("input.roi_path != 'whole_tissue'",
      plotOutput("roi_selection_plot", height = 560, brush = brushOpts("roi_brush", resetOnNew = TRUE))),
    uiOutput("roi_gate"), DTOutput("roi_summary"), uiOutput("roi_download"), uiOutput("next_from_roi")),
  tabPanel("5. Domains & niches", value = "spatial_structure", uiOutput("spatial_gate"), h3("Metabolic domains"), fileInput("domain_csv_runtime", "Add/replace domain CSV", accept = ".csv"),
    tags$div(class = "alert alert-info",
      "Domains classify pixels by their own MSI molecular profiles. With no prior ROI, domains are fitted across the tissue mask; with a prior ROI, they are fitted only inside that field."),
    tags$details(tags$summary("Example datasets for domains and niches"),
      tags$p("OMIX016317 validates full-tissue exploratory domains; MSIflow UPEC_12 provides previously generated study domain labels."),
      actionButton("use_omix_domains", "View OMIX016317 domain example", class = "btn-info"),
      actionButton("use_msiflow_domains", "Run MSIflow domain example", class = "btn-info")),
    uiOutput("domain_example_results"),
    numericInput("domain_k", "Exploratory domain count k", 4, min = 2, max = 12),
    numericInput("domain_seed", "Random seed", 20260808, min = 1), numericInput("domain_pcs", "PCA components", 10, min = 2, max = 30),
    selectInput("domain_feature_source", "Features used for domain detection",
      c("All processed m/z features" = "all", "Moran FDR-significant features" = "moran_fdr", "Top Moran-ranked features" = "moran_top")),
    fluidRow(
      column(6, numericInput("domain_moran_fdr", "Moran FDR threshold", 0.05, min = 0, max = 1, step = 0.01)),
      column(6, numericInput("domain_moran_top_n", "Top Moran features", 200, min = 2, step = 10))),
    uiOutput("domain_action"), tags$div(class = "alert alert-warning",
      "Pseudoreplication warning: pixels are not biological replicates. Generated domains are descriptive and data-driven, never anatomical ROI."),
    uiOutput("domain_gate"), plotOutput("domain_plot", height = 600), DTOutput("domain_counts"), DTOutput("domain_features"), uiOutput("domain_download"),
    uiOutput("domain_compare_selector"), uiOutput("domain_compare_summary"),
    plotOutput("domain_compare_plot", height = 520), DTOutput("domain_compare_features"),
    uiOutput("domain_roi_selector"), uiOutput("domain_roi_gate"),
    hr(), h3("Local domain-composition niches"),
    tags$div(class = "alert alert-info",
      "Optional: niches group pixels by the proportions of surrounding domains, revealing domain interiors, mixed boundaries and transition environments. They do not replace ROI selection and are not cell types, cellular niches, or independent biological replicates."),
    fluidRow(
      column(4, selectInput("niche_neighborhood", "Neighborhood", c("Radius (recommended)" = "radius", "k nearest neighbors" = "knn"))),
      column(4, numericInput("niche_radius", "Radius (coordinate units)", 2.01, min = .Machine$double.eps, step = 0.5)),
      column(4, numericInput("niche_k", "k neighbors", 10, min = 1, step = 1))
    ),
    fluidRow(
      column(3, numericInput("niche_min_neighbors", "Minimum neighbors", 5, min = 1, step = 1)),
      column(3, numericInput("niche_count", "Number of niches", 4, min = 2, step = 1)),
      column(3, selectInput("niche_transform", "Composition transform",
        c("Hellinger (recommended)" = "hellinger", "None" = "none", "CLR" = "clr"))),
      column(3, checkboxInput("niche_include_self", "Include focal pixel", TRUE))
    ),
    fluidRow(
      column(6, selectInput("niche_domain_alignment", "Domain-label provenance",
        c("Domains fitted jointly" = "joint", "Domains aligned across fields" = "aligned"))),
      column(6, numericInput("niche_seed", "Random seed", 20260810, min = 1, step = 1))
    ),
    uiOutput("niche_action"), uiOutput("niche_gate"),
    plotOutput("niche_plot", height = 600), plotOutput("niche_ambiguity_plot", height = 600),
    uiOutput("niche_compare_selector"), uiOutput("niche_compare_summary"),
    plotOutput("niche_compare_plot", height = 520), DTOutput("niche_compare_composition"),
    DTOutput("niche_compare_features"),
    DTOutput("niche_exclusion_summary"), uiOutput("niche_download"), uiOutput("next_from_structure")),
  tabPanel("6. ROI summaries & statistics", value = "spatial_statistics", uiOutput("spatial_gate"),
    h3("ROI/subregion summaries"),
    tags$div(class = "alert alert-warning",
      "Subregions are technical spatial summaries, not biological replicates. Population inference requires independently sampled subjects."),
    fluidRow(
      column(4, numericInput("sampling_grid_size", "Grid divisions per axis", 5, min = 2, step = 1)),
      column(4, numericInput("sampling_min_pixels", "Minimum pixels per subregion", 30, min = 1, step = 1)),
      column(4, selectInput("sampling_grid_scope", "Grid scope", c("Within each ROI" = "roi", "Global field" = "global")))
    ),
    checkboxInput("acknowledge_exploratory_roi",
      "I understand that statistics using data-driven domains from the same features are exploratory, not independent confirmation.", FALSE),
    uiOutput("sampling_action"), uiOutput("sampling_gate"), uiOutput("sampling_results_ui"),
    hr(), h3("Import statistical result and map it back"),
    tags$div(class = "alert alert-info",
      "Upload a MetaboAnalyst VIP/differential/PCA/PLS-DA CSV. Score rows must retain the exported Sample identifiers for spatial back-mapping."),
    fileInput("metabo_result_csv", "MetaboAnalyst result CSV", accept = ".csv"),
    uiOutput("metabo_result_controls"), uiOutput("metabo_result_gate"),
    uiOutput("metabo_results_ui"), uiOutput("next_from_statistics")),
  tabPanel("7. Spatial validation", value = "spatial_validation", uiOutput("spatial_gate"),
    h3("Differential-result spatial coherence check"),
    tags$p("This reuses the Moran's I result computed after preprocessing. It does not recompute Moran's I and is not an independent validation when Moran-filtered features were used to create the domains."),
    tags$details(tags$summary("Technical details — load a previously generated result"),
      textInput("moran_result_dir", "Result directory", Sys.getenv("SPATIALOMICS_OMIX_MORAN_DIR", "")), actionButton("load_moran", "Load result")),
    uiOutput("moran_validation_gate"), uiOutput("moran_concordance_ui"),
    tags$details(class = "result-card", tags$summary("Optional advanced check — compare domains with an independent label map"),
      tags$p("Use this only when you have a second coordinate-aligned map created independently of the current MSI clustering—for example manual anatomy, pathology labels, or another modality. It measures agreement; it is not another ROI-selection step."),
      fileInput("corroboration_csv", "Independent label CSV (x, y and a label column)", accept = ".csv"),
      selectInput("corroboration_mapping", "Mapping rule",
        c("Mutual best match" = "mutual_best", "One-way best match" = "one_way")),
      fluidRow(
        column(6, numericInput("corroboration_fraction", "Minimum conditional overlap", 0.5, min = 0, max = 1, step = 0.05)),
        column(6, numericInput("corroboration_count", "Minimum shared pixels", 10, min = 1, step = 1))),
      uiOutput("corroboration_action"), uiOutput("corroboration_gate"), uiOutput("corroboration_results_ui")),
    uiOutput("next_from_validation")),
  tabPanel("8. Matched transcriptomics & H&E", value = "matched_omics",
    uiOutput("matched_omics_gate"),
    tags$p("These data support or interpret MSI findings. They do not replace the MSI analysis and are only enabled when the selected tissue section has matched files."),
    tags$div(class = "matched-layout",
      tags$div(class = "result-card", h4("Matched H&E image"),
        tags$div(class = "matched-image", imageOutput("matched_histology", height = "100%"))),
      tags$div(class = "matched-tables",
        tags$div(class = "result-card", h4("Matched files"), DTOutput("matched_file_inventory")),
        tags$div(class = "result-card", h4("Annotation summary"), DTOutput("matched_annotation_summary")))),
    uiOutput("next_from_matched")),
  tabPanel("9. LC-MS/MS evidence (optional)", value = "lcms_evidence", uiOutput("lcms_workflow_gate"),
    tags$div(class = "alert alert-info",
      "This is not a standalone LC-MS/MS analysis. It evaluates chemical-identification evidence for features discovered in the current MSI analysis."),
    tags$div(class = "well", strong("Example evidence: "),
      actionButton("use_msv_evidence", "Run MSV000090179 external-reference example", class = "btn-info"),
      actionButton("use_msiflow_evidence", "View MSIflow study-matched evidence", class = "btn-info"),
      actionButton("use_kasarla_validation", "Run Kasarla 2025 kidney ROI validation", class = "btn-success"),
      tags$p("MSV supplies external chemical evidence; MSIflow is study-matched; Kasarla compares MALDI-MSI and LMD-LC-MS/MS across cortex, medulla and renal pelvis of serial kidney sections.")),
    uiOutput("lcms_example_results"),
    uiOutput("kasarla_validation_results"),
    uiOutput("lcms_relationship_notice"),
    numericInput("msi_mz", "Candidate MSI m/z", 775.55261535, step = .00000001),
    numericInput("lcms_precursor", "LC-MS/MS precursor target", 775.550137928655, step = .000000000001), uiOutput("lcms_action"),
    tags$p("Precursor-level match and an unassigned fragment spectrum are shown. Chemical identity is never inferred without user-supplied diagnostic-ion definitions."),
    uiOutput("lcms_gate"), DTOutput("precursor_table"), plotOutput("fragment_plot", height = 520),
    hr(), h3("MSI–LC-MS feature-table matching"),
    tags$p("Both CSV files require an mz column. Optional shared ion_mode and log2fc columns add polarity filtering and direction agreement. Matching is one-to-one."),
    fluidRow(column(6, fileInput("msi_feature_csv", "MSI feature CSV", accept = ".csv")),
      column(6, fileInput("lcms_feature_csv", "LC-MS feature CSV", accept = ".csv"))),
    fluidRow(column(6, numericInput("feature_match_ppm", "m/z tolerance (ppm)", 5, min = 0, step = 0.5)),
      column(6, selectInput("assignment_method", "Assignment", c("Optimal (recommended)" = "optimal", "Greedy" = "greedy")))),
    actionButton("run_feature_matching", "Match feature tables", class = "btn-primary"),
    uiOutput("feature_match_gate"), DTOutput("feature_match_table"), uiOutput("feature_match_download"),
    hr(), h3("CCS / ion-mobility evidence"),
    tags$p("Observed CSV: candidate_id, observed_ccs, observed_source. Reference CSV: candidate_id, reference_ccs, reference_source. CCS supports a candidate but is not identity proof."),
    fluidRow(column(6, fileInput("observed_ccs_csv", "Observed CCS CSV", accept = ".csv")),
      column(6, fileInput("reference_ccs_csv", "Reference CCS CSV", accept = ".csv"))),
    numericInput("ccs_tolerance", "CCS tolerance (%)", 2, min = 0.01, step = 0.1),
    actionButton("run_ccs_validation", "Validate CCS evidence", class = "btn-primary"),
    uiOutput("ccs_gate"), DTOutput("ccs_table"), uiOutput("ccs_download"), uiOutput("next_from_lcms"))
  ),
  tabPanel("10. Results & export", value = "results", h2("Session results"),
    uiOutput("workflow_progress_results"), uiOutput("result_inventory"),
    tags$p("Each session writes to a unique temporary directory; source data are read-only."),
    uiOutput("download_gate"), downloadButton("download_bundle", "Download session bundle"),
    tags$details(tags$summary("Technical details — temporary session location"), textOutput("session_path")))
)

server <- function(input, output, session) {
  state <- reactiveValues(valid = FALSE, spec = NULL, validation = NULL, processed = NULL,
    analysis_matrix = NULL, tissue_gate = NULL, tissue_mask = NULL, registration = NULL,
    registration_points = NULL, rois = NULL, roi_source = NULL, domains = NULL, domain_source = NULL, domain_counts = NULL, domain_features = NULL, neighbors = NULL,
    niches = NULL,
    sampling = NULL, metabo_input = NULL, functional_peak_table = NULL,
    functional_peak_list = NULL, metabo_result = NULL,
    metabo_result_type = NULL, mapped_scores = NULL,
    corroboration = NULL, corroborated_pixels = NULL, moran = NULL, lcms = NULL,
    feature_matches = NULL, ccs_evidence = NULL, example_note = NULL,
    runtime_lcms_path = NULL, runtime_lcms_relationship = NULL,
    preset_key = NULL, qc_example = NULL, registration_example = NULL, domain_example = NULL,
    lcms_example = NULL, kasarla_validation = NULL, kasarla_bundle = NULL,
    lipid_export = NULL, lipid_name_list = NULL,
    kasarla_raw_msi_available = FALSE, kasarla_pixel_msi_available = FALSE)
  session_dir <- file.path(session_root, paste0("SpatialOmicsMSI-session-", Sys.getpid(), "-", substr(session$token, 1, 8)))
  dir.create(session_dir, recursive = TRUE, showWarnings = FALSE)
  output$session_path <- renderText(session_dir)

  observeEvent(input$choose_spatial, {
    updateNavbarPage(session, "workflow_tabs", selected = "data")
  })

  upload_path <- function(value) if (is.null(value)) "" else value$datapath
  current_spec <- reactive({
    server <- identical(input$data_origin, "server")
    list(msi_input_type = input$msi_input_type,
      imzml_path = if (server) input$server_imzml %or% "" else upload_path(input$upload_imzml),
      ibd_path = if (server) input$server_ibd %or% "" else upload_path(input$upload_ibd),
      csv_path = if (server) input$server_csv %or% "" else upload_path(input$upload_csv),
      optical_path = if (server) input$server_optical %or% "" else upload_path(input$upload_optical),
      transform_path = if (server) input$server_transform %or% "" else upload_path(input$upload_transform),
      labels_path = if (server) input$server_labels %or% "" else upload_path(input$upload_labels),
      lcms_path = if (server) input$server_lcms %or% "" else upload_path(input$upload_lcms),
      attribution_path = if (server) input$server_attribution %or% "" else "",
      sample_id = input$sample_id %or% "", section_id = input$section_id %or% "",
      subject_id = input$subject_id %or% "", ion_mode = input$ion_mode %or% "",
      ion_mode_source = input$ion_mode_source %or% "", lcms_relationship = input$lcms_relationship %or% "none",
      spectral_processing = input$spectral_processing,
      alignment_ppm = input$alignment_ppm, peak_pick_snr = input$peak_pick_snr,
      min_detection_fraction = input$min_detection_fraction, make_tissue_mask = isTRUE(input$make_tissue_mask))
  })

  matched_paths <- function(spec = state$spec) {
    if (is.null(spec)) return(character())
    fields <- c("st_expression_path", "st_positions_path", "histology_path",
      "st_region_path", "st_lesion_path", "st_loupe_path")
    paths <- unlist(spec[fields], use.names = TRUE)
    if (!length(paths)) return(character())
    paths[vapply(paths, present_path, logical(1)) & file.exists(paths)]
  }

  reset_analysis <- function() {
    state$valid <- FALSE; state$validation <- NULL; state$processed <- NULL; state$analysis_matrix <- NULL
    state$tissue_gate <- NULL; state$tissue_mask <- NULL; state$registration <- NULL; state$registration_points <- NULL
    state$rois <- NULL; state$roi_source <- NULL; state$domains <- NULL; state$domain_source <- NULL; state$domain_counts <- NULL; state$domain_features <- NULL
    state$niches <- NULL
    state$sampling <- NULL; state$metabo_input <- NULL; state$functional_peak_table <- NULL
    state$functional_peak_list <- NULL; state$metabo_result <- NULL
    state$lipid_export <- NULL; state$lipid_name_list <- NULL
    state$metabo_result_type <- NULL; state$mapped_scores <- NULL
    state$corroboration <- NULL; state$corroborated_pixels <- NULL
    state$neighbors <- NULL; state$moran <- NULL; state$lcms <- NULL
    state$feature_matches <- NULL; state$ccs_evidence <- NULL
    state$runtime_lcms_path <- NULL; state$runtime_lcms_relationship <- NULL
    state$kasarla_validation <- NULL; state$kasarla_bundle <- NULL
    state$kasarla_raw_msi_available <- FALSE; state$kasarla_pixel_msi_available <- FALSE
  }

  observeEvent(input$load_example_dataset, {
    key <- input$example_dataset
    spec <- preset_spec(key)
    reset_analysis()
    state$preset_key <- key
    state$spec <- spec
    state$example_note <- spec$source_note
    if (identical(example_workflow(spec), "lcms")) {
      load_msv_validation()
      updateNavbarPage(session, "workflow_tabs", selected = "lcms_evidence")
      return()
    }
    updateRadioButtons(session, "data_origin", selected = "server")
    updateSelectInput(session, "msi_input_type", selected = spec$msi_input_type)
    updateTextInput(session, "server_csv", value = spec$csv_path %or% "")
    updateTextInput(session, "server_imzml", value = spec$imzml_path %or% "")
    updateTextInput(session, "server_ibd", value = spec$ibd_path %or% "")
    updateTextInput(session, "sample_id", value = spec$sample_id)
    updateTextInput(session, "section_id", value = spec$section_id)
    updateTextInput(session, "subject_id", value = spec$subject_id)
    updateSelectInput(session, "ion_mode", selected = spec$ion_mode)
    updateTextInput(session, "ion_mode_source", value = spec$ion_mode_source)
    updateSelectInput(session, "lcms_relationship", selected = spec$lcms_relationship %or% "none")
    state$validation <- withProgress(message = paste("Loading", spec$label), value = 0, {
      incProgress(.15, detail = "Reading the pixel-level MSI table")
      out <- validate_analysis_spec(spec, inspect = FALSE)
      incProgress(.70, detail = "Checking coordinates and feature metadata")
      if (identical(key, "kasarla_kidney")) {
        incProgress(.10, detail = "Attaching matched LMD-LC-MS/MS evidence")
        load_kasarla_validation()
      }
      incProgress(.05, detail = "Example dataset is ready")
      out
    })
    state$valid <- state$validation$valid
    if (isTRUE(state$valid)) {
      write.csv(data.frame(workflow = "spatial_msi_with_optional_matched_omics",
        preset = key, sample_id = spec$sample_id, section_id = spec$section_id,
        matched_companion_files = length(matched_paths(spec)), stringsAsFactors = FALSE),
        file.path(session_dir, "workflow_manifest.csv"), row.names = FALSE)
      updateNavbarPage(session, "workflow_tabs", selected = "spatial_provenance")
    }
  })

  active_lcms_path <- reactive(state$runtime_lcms_path %or% if (is.null(state$spec)) "" else state$spec$lcms_path %or% "")
  active_lcms_relationship <- reactive(state$runtime_lcms_relationship %or%
    if (is.null(state$spec)) "none" else state$spec$lcms_relationship %or% "none")
  load_msv_validation <- function() {
    path <- example_datasets$lcms$lcms_path
    precursor <- input$lcms_precursor %or% 775.550137928655
    result <- withProgress(message = "Reading MSV000090179 example", value = .1, {
      result <- read_mzml_fragment_spectra(path, precursor, 10); incProgress(.9); result
    })
    result$example_key <- "msv"
    result$relationship <- "external_reference"
    state$lcms_example <- result
    state$lcms <- result
    state$runtime_lcms_path <- path
    state$runtime_lcms_relationship <- "external_reference"
  }
  observeEvent(input$use_msv_evidence, load_msv_validation())
  observeEvent(input$use_msiflow_evidence, {
    scans <- file.path(project_root, "results/real_data/lcms_validation/Sample_2_rep1_Pos_1/target_precursor_ms2_scans.csv")
    fragments <- file.path(project_root, "results/real_data/lcms_validation/Sample_2_rep1_Pos_1/fragment_centroid_table.csv")
    validate(need(file.exists(scans) && file.exists(fragments), "Saved MSIflow evidence example files are unavailable."))
    result <- list(precursor_scan_metadata = read.csv(scans, check.names = FALSE),
      fragment_peak_table = read.csv(fragments, check.names = FALSE))
    result$example_key <- "msiflow"
    result$relationship <- "study_matched"
    state$lcms_example <- result
    state$lcms <- result
    state$runtime_lcms_path <- example_datasets$msiflow$lcms_path
    state$runtime_lcms_relationship <- "study_matched"
  })
  load_kasarla_validation <- function() {
    directory <- file.path(project_root, "data_raw/kasarla2025/processed")
    filenames <- c("kasarla_kidney_maldi_roi_matrix.csv",
      "kasarla_kidney_lmd_lcms_roi_matrix.csv", "kasarla_kidney_feature_mapping.csv")
    paths <- file.path(directory, filenames)
    if (!all(file.exists(paths))) {
      directory <- system.file("extdata", "kasarla2025", package = "SpatialOmicsMSI")
      paths <- file.path(directory, filenames)
    }
    validate(need(all(file.exists(paths)), "Kasarla 2025 processed validation tables are unavailable."))
    msi <- read.csv(paths[1], check.names = FALSE, stringsAsFactors = FALSE)
    lcm <- read.csv(paths[2], check.names = FALSE, stringsAsFactors = FALSE)
    mapping <- read.csv(paths[3], check.names = FALSE, stringsAsFactors = FALSE)
    comparison <- compare_msi_lcm_quantification(msi, lcm, mapping,
      id_columns = c("roi_id", "section_id"), min_pairs = 3, method = "spearman")
    comparison$metabolite_name <- mapping$metabolite_name[match(comparison$msi_feature, mapping$msi_feature)]
    state$kasarla_validation <- comparison
    manifest_candidates <- c(
      file.path(project_root, "inst/extdata/kasarla2025/kasarla_kidney_bundle_manifest.csv"),
      system.file("extdata", "kasarla2025", "kasarla_kidney_bundle_manifest.csv", package = "SpatialOmicsMSI"))
    manifest_path <- manifest_candidates[file.exists(manifest_candidates)][1]
    validate(need(length(manifest_path) == 1L && !is.na(manifest_path),
      "Kasarla kidney bundle manifest is unavailable."))
    state$kasarla_bundle <- read.csv(manifest_path, check.names = FALSE, stringsAsFactors = FALSE)
    msi_directory <- file.path(project_root, "data_raw/kasarla2025/metaspace/2023-12-19_13h02m39s")
    state$kasarla_raw_msi_available <- dir.exists(msi_directory) &&
      length(list.files(msi_directory, pattern = "[.]imzML$", ignore.case = TRUE)) == 1L &&
      length(list.files(msi_directory, pattern = "[.]ibd$", ignore.case = TRUE)) == 1L
    state$kasarla_pixel_msi_available <- state$kasarla_raw_msi_available || any(file.exists(file.path(
      project_root, "data_raw/kasarla2025/processed",
      c("kasarla_kidney_pixel_feature_matrix.rds", "kasarla_kidney_pixel_feature_matrix.csv.gz"))))
    write.csv(comparison, file.path(session_dir, "kasarla2025_kidney_maldi_lmd_lcms_concordance.csv"), row.names = FALSE)
  }
  observeEvent(input$use_kasarla_validation, load_kasarla_validation())
  observeEvent(input$open_kasarla_validation, {
    load_kasarla_validation()
    updateNavbarPage(session, "workflow_tabs", selected = "lcms_evidence")
  })

  shiny_asset_path <- function(name) {
    source_path <- file.path(project_root, "inst/shiny/spatial_pipeline/www", name)
    if (file.exists(source_path)) return(source_path)
    installed_path <- system.file("shiny", "spatial_pipeline", "www", name, package = "SpatialOmicsMSI")
    installed_path
  }
  activate_msiflow_msi_example <- function() {
    spec <- preset_spec("msiflow")
    state$spec <- spec
    state$validation <- validate_analysis_spec(spec, inspect = FALSE)
    state$valid <- state$validation$valid
    state$processed <- read_pixel_feature_csv(spec$csv_path, spec$sample_id, spec$section_id, spec$subject_id)
    feature_names <- state$processed$feature_metadata$column_name
    matrix <- as.matrix(state$processed$pixel_feature_matrix[, feature_names, drop = FALSE])
    storage.mode(matrix) <- "double"
    tic <- rowSums(matrix, na.rm = TRUE)
    state$analysis_matrix <- log10(matrix / pmax(tic, .Machine$double.eps) + 1)
    state$tic <- tic
    state$tissue_mask <- rep(TRUE, nrow(matrix))
    state$example_note <- spec$source_note
    invisible(state$processed)
  }
  observeEvent(input$use_brain01_qc, {
    state$qc_example <- list(key = "brain01", image = shiny_asset_path("brain01_tic_overlay.png"), table = NULL,
      interpretation = "Brain01 centroid MSI: TIC distribution shown in registered tissue coordinates.")
  })
  observeEvent(input$use_omix_qc, {
    state$qc_example <- list(key = "omix", image = shiny_asset_path("omix016317_full_field_tissue_gate.png"), table = NULL,
      interpretation = "OMIX016317: full acquisition field and tissue-gating QC.")
  })
  observeEvent(input$use_msiflow_qc, {
    activate_msiflow_msi_example()
    state$qc_example <- list(key = "msiflow", image = NULL, table = state$processed$qc_summary,
      interpretation = "MSIflow UPEC_12 was loaded and processed directly in this module using TIC normalization and log10(x + 1).")
  })
  observeEvent(input$use_brain01_registration, {
    state$registration_example <- list(key = "brain01",
      image = shiny_asset_path("brain01_registered_mask_overlay.png"),
      interpretation = "Brain01 optical image and MSI measurement mask after applying the supplied transform.")
    state$registration <- list(
      diagnostics = data.frame(example = "Brain01", transform = "supplied 3 x 3 METASPACE transform",
        interpretation = "real-data registration example", stringsAsFactors = FALSE),
      output_files = c(measurement_mask_overlay = shiny_asset_path("brain01_registered_mask_overlay.png")))
  })
  observeEvent(input$use_omix_domains, {
    state$domain_example <- list(key = "omix", image = shiny_asset_path("omix016317_tissue_only_domains.png"),
      pixels = NULL, interpretation = "OMIX016317 exploratory metabolic domains within the tissue gate.")
  })
  observeEvent(input$use_msiflow_domains, {
    path <- file.path(project_root, "results/real_data/UPEC_12/spatial_domain_labels.csv")
    validate(need(file.exists(path), "Saved MSIflow domain example is unavailable."))
    activate_msiflow_msi_example()
    pixels <- read.csv(path, check.names = FALSE)
    state$domains <- pixels
    state$domain_source <- "data_driven"
    state$domain_counts <- as.data.frame(table(pixels$domain_label)); names(state$domain_counts) <- c("domain", "pixel_count")
    state$domain_example <- list(key = "msiflow", image = NULL, pixels = pixels,
      interpretation = "MSIflow UPEC_12 data-driven domains; labels are descriptive rather than anatomical ground truth.")
  })

  output$qc_example_status <- renderUI({ req(state$qc_example)
    tags$div(class = "alert alert-success", state$qc_example$interpretation) })
  output$qc_example_results <- renderUI({
    if (is.null(state$qc_example)) return(NULL)
    tagList(uiOutput("qc_example_status"),
      if (!is.null(state$qc_example$image)) imageOutput("qc_example_image", height = "600px"),
      if (!is.null(state$qc_example$table)) DTOutput("qc_example_table"))
  })
  output$qc_example_image <- renderImage({ req(state$qc_example$image)
    list(src = state$qc_example$image, contentType = "image/png", alt = state$qc_example$interpretation)
  }, deleteFile = FALSE)
  output$qc_example_table <- renderDT({ req(state$qc_example$table)
    datatable(state$qc_example$table, options = list(dom = "t", scrollX = TRUE)) })
  output$registration_example_status <- renderUI({ req(state$registration_example)
    tags$div(class = "alert alert-success", state$registration_example$interpretation) })
  output$registration_example_results <- renderUI({
    if (is.null(state$registration_example)) return(NULL)
    tagList(uiOutput("registration_example_status"),
      if (!is.null(state$registration_example$image)) imageOutput("registration_example_image", height = "650px"))
  })
  output$registration_example_image <- renderImage({ req(state$registration_example$image)
    list(src = state$registration_example$image, contentType = "image/png", alt = state$registration_example$interpretation)
  }, deleteFile = FALSE)
  output$domain_example_status <- renderUI({ req(state$domain_example)
    tags$div(class = "alert alert-success", state$domain_example$interpretation) })
  output$domain_example_results <- renderUI({
    if (is.null(state$domain_example)) return(NULL)
    tagList(uiOutput("domain_example_status"),
      if (!is.null(state$domain_example$image)) imageOutput("domain_example_image", height = "600px"),
      if (!is.null(state$domain_example$pixels)) tagList(
        plotOutput("domain_example_plot", height = 600), DTOutput("domain_example_table")))
  })
  output$domain_example_image <- renderImage({ req(state$domain_example$image)
    list(src = state$domain_example$image, contentType = "image/png", alt = state$domain_example$interpretation)
  }, deleteFile = FALSE)
  output$domain_example_plot <- renderPlot({ req(state$domain_example$pixels); d <- state$domain_example$pixels
    ggplot(d, aes(x, y, color = domain_label)) + geom_point(size = .6) + coord_equal() + scale_y_reverse() +
      theme_minimal() + labs(title = "MSIflow UPEC_12 domain example", color = "Domain") })
  output$domain_example_table <- renderDT({ req(state$domain_example$pixels); d <- state$domain_example$pixels
    counts <- as.data.frame(table(d$domain_label)); names(counts) <- c("domain", "pixel_count")
    datatable(counts, options = list(dom = "t")) })
  output$lcms_example_status <- renderUI({ req(state$lcms_example)
    tags$div(class = "alert alert-success", paste("Example loaded directly in this module. Evidence relationship:",
      state$lcms_example$relationship)) })
  output$lcms_example_results <- renderUI({
    if (is.null(state$lcms_example)) return(NULL)
    tagList(uiOutput("lcms_example_status"), DTOutput("lcms_example_table"),
      plotOutput("lcms_example_plot", height = 500))
  })
  output$lcms_example_table <- renderDT({ req(state$lcms_example)
    datatable(state$lcms_example$precursor_scan_metadata, options = list(pageLength = 10, scrollX = TRUE)) })
  output$lcms_example_plot <- renderPlot({ req(state$lcms_example); d <- state$lcms_example$fragment_peak_table
    validate(need(nrow(d) > 0L, "No fragment peaks are available for this example target."))
    mz_col <- if ("fragment_mz" %in% names(d)) "fragment_mz" else "mz"
    intensity_col <- if ("fragment_intensity" %in% names(d)) "fragment_intensity" else "intensity"
    ggplot(d, aes(x = .data[[mz_col]], y = .data[[intensity_col]])) +
      geom_segment(aes(xend = .data[[mz_col]], yend = 0)) + theme_minimal() +
      labs(title = "Example LC-MS/MS fragment evidence", x = "fragment m/z", y = "intensity")
  })
  output$kasarla_validation_results <- renderUI({
    if (is.null(state$kasarla_validation)) return(NULL)
    d <- state$kasarla_validation
    tagList(tags$div(class = "alert alert-success", sprintf(
      "Kasarla 2025 loaded: %d matched metabolites; %d show positive three-region MALDI/LMD Spearman concordance. Supporting tables are CC BY-NC 4.0.",
      nrow(d), sum(is.finite(d$correlation) & d$correlation > 0))),
      tags$div(class = if (isTRUE(state$kasarla_pixel_msi_available)) "alert alert-success" else "alert alert-warning",
        strong("Selected MSI: "), "METASPACE 2023-12-19_13h02m39s · 231218_SK_labeledkidney_HEXANE_NH4OHwash_NEDC_MALDI_40um. ",
        if (isTRUE(state$kasarla_raw_msi_available))
          "The raw imzML/ibd pair is installed."
        else if (isTRUE(state$kasarla_pixel_msi_available))
          "The 45,623-pixel × 236-feature METASPACE matrix is installed, so the spatial workflow starts at Step 1; raw-spectrum peak detection is the only unavailable preprocessing branch."
        else "Only the published ROI validation tables are installed; the pixel-level spatial workflow is unavailable."),
      fluidRow(
        column(5, tags$div(class = "result-card", plotOutput("kasarla_validation_plot", height = 430))),
        column(7, tags$div(class = "result-card", DTOutput("kasarla_validation_table")))),
      tags$div(class = "result-card", h4("Matched kidney bundle inventory"),
        tags$p("One MALDI-MSI input plus 12 LMD-LC-MS/MS files; the table records the exact repository identifiers."),
        DTOutput("kasarla_bundle_table")))
  })
  output$kasarla_validation_plot <- renderPlot({
    req(state$kasarla_validation); d <- state$kasarla_validation
    ggplot(d, aes(x = correlation)) + geom_histogram(binwidth = .25, boundary = 0, color = "white") +
      geom_vline(xintercept = 0, linetype = 2, color = "#B2182B") + theme_minimal() +
      labs(title = "Kidney ROI concordance", subtitle = "Cortex · medulla · renal pelvis", x = "Spearman correlation", y = "Metabolites")
  })
  output$kasarla_validation_table <- renderDT({
    req(state$kasarla_validation)
    datatable(state$kasarla_validation[order(state$kasarla_validation$correlation, decreasing = TRUE), ],
      options = list(pageLength = 10, scrollX = TRUE))
  })
  output$kasarla_bundle_table <- renderDT({
    req(state$kasarla_bundle)
    datatable(state$kasarla_bundle, options = list(pageLength = 13, scrollX = TRUE, dom = "tip"))
  })

  output$active_dataset_card <- renderUI({
    if (!isTRUE(state$valid) || is.null(state$spec))
      return(tags$div(class = "alert alert-warning", "Choose and validate a dataset on Data setup first."))
    tags$div(class = "alert alert-info",
      strong("Active dataset: "), state$spec$sample_id,
      " · section: ", state$spec$section_id,
      " · input: ", if (identical(state$spec$msi_input_type, "csv")) basename(state$spec$csv_path) else basename(state$spec$imzml_path),
      if (!is.null(state$preset_key)) tags$span(" · preset: ", state$preset_key))
  })

  observeEvent(input$validate_input, {
    candidate <- current_spec()
    if (!is.null(state$preset_key)) {
      preset <- preset_spec(state$preset_key)
      if (identical(candidate$sample_id, preset$sample_id) &&
          identical(candidate$section_id, preset$section_id)) {
        companion_fields <- c("st_expression_path", "st_positions_path", "histology_path",
          "st_region_path", "st_lesion_path", "st_loupe_path", "source_note")
        candidate[companion_fields] <- preset[companion_fields]
      } else state$preset_key <- NULL
    }
    reset_analysis(); state$spec <- candidate
    state$validation <- withProgress(message = "Checking input metadata", value = .2, {
      out <- validate_analysis_spec(state$spec, inspect = TRUE); incProgress(.8); out
    })
    state$valid <- state$validation$valid
    if (isTRUE(state$valid)) write.csv(data.frame(
      workflow = "spatial_msi_with_optional_lcms_evidence",
      has_msi = isTRUE(state$validation$capabilities$msi),
      has_lcms = isTRUE(state$validation$capabilities$lcms),
      lcms_relationship = state$spec$lcms_relationship,
      stringsAsFactors = FALSE
    ), file.path(session_dir, "workflow_manifest.csv"), row.names = FALSE)
  })

  output$workflow_status <- renderUI({
    tagList(tags$span(class = "label label-primary", "Spatial MSI"),
      tags$span(class = paste("label", if (isTRUE(state$valid)) "label-success" else "label-default"),
        if (isTRUE(state$valid)) "input validated" else "input not validated"),
      if (!is.null(state$processed)) tags$span(class = "label label-success", "MSI processed"),
      if (!is.null(state$domains)) tags$span(class = "label label-success", "domains ready"),
      if (!is.null(state$niches)) tags$span(class = "label label-success", "niches ready"),
      if (length(matched_paths())) tags$span(
        class = paste("label", if (!is.null(state$registration)) "label-success" else "label-warning"),
        if (!is.null(state$registration)) "matched omics registered" else "matched files · registration required"),
      if (!is.null(state$lcms)) tags$span(class = "label label-success", "LC-MS/MS read"),
      if (!is.null(state$feature_matches)) tags$span(class = "label label-success", "features matched"))
  })
  progress_ui <- reactive({
    items <- data.frame(
      stage = c("Input validated", "MSI processed", "ROI selected", "Domains available", "Spatial niches", "Statistical export", "Spatial back-map", "Matched omics registered", "LC-MS/MS evidence", "MSI–LC-MS/MS matches"),
      complete = c(isTRUE(state$valid), !is.null(state$processed),
        !is.null(state$rois) && any(!is.na(state$rois$roi_id)), !is.null(state$domains), !is.null(state$niches),
        !is.null(state$metabo_input), !is.null(state$mapped_scores),
        length(matched_paths()) > 0L && !is.null(state$registration),
        !is.null(state$lcms), !is.null(state$feature_matches)),
      stringsAsFactors = FALSE)
    tags$ul(lapply(seq_len(nrow(items)), function(i) tags$li(
      tags$span(class = paste("label", if (items$complete[i]) "label-success" else "label-default"),
        if (items$complete[i]) "complete" else "pending"), " ", items$stage[i])))
  })
  output$workflow_progress_home <- renderUI(progress_ui())
  output$workflow_progress_results <- renderUI(progress_ui())
  output$workflow_next_step <- renderUI({
    if (!isTRUE(state$valid)) return(tags$div(class = "alert alert-warning", "Next: provide MSI input and validate it."))
    if (is.null(state$processed))
      return(tags$div(class = "alert alert-info", "Next: run MSI processing and review QC."))
    if (is.null(state$domains))
      return(tags$div(class = "alert alert-info", "Next: optionally define a prior ROI, or detect metabolic domains across the complete tissue."))
    if (is.null(state$rois) || !any(!is.na(state$rois$roi_id)))
      return(tags$div(class = "alert alert-info", "Next: select one or more detected domains as the downstream ROI, or define another ROI."))
    tags$div(class = "alert alert-success", "The selected workflow has usable results. Review Results & export for the session bundle.")
  })
  output$spatial_gate <- renderUI({
    if (!isTRUE(state$valid) || !isTRUE(state$validation$capabilities$msi))
      tags$div(class = "alert alert-danger", "Spatial MSI modules require a validated imzML/ibd pair or pixel-by-feature CSV.")
  })
  output$lcms_workflow_gate <- renderUI({
    if (!isTRUE(state$valid) || is.null(state$processed))
      return(tags$div(class = "alert alert-danger", "Process the MSI data before evaluating LC-MS/MS evidence."))
    if (!present_path(active_lcms_path()))
      tags$div(class = "alert alert-warning", "No LC-MS/MS mzML evidence file was supplied. Feature-table and CCS evidence can still be added, but their provenance must be reported.")
  })
  output$lcms_relationship_notice <- renderUI({
    req(state$spec)
    labels <- c(none = "No LC-MS/MS relationship was declared.",
      specimen_matched = "Specimen-matched evidence may support identity and within-design concordance, subject to the study design.",
      study_matched = "Study-matched evidence supports chemical identity only; specimen-level abundance concordance must not be claimed.",
      external_reference = "External-reference evidence supports chemical identity only and is not biological validation of the MSI cohort.")
    relationship <- active_lcms_relationship()
    tags$div(class = paste("alert", if (relationship == "specimen_matched") "alert-success" else "alert-warning"), labels[[relationship]])
  })
  output$result_inventory <- renderUI({
    files <- list.files(session_dir, recursive = TRUE)
    if (!length(files)) return(tags$div(class = "alert alert-secondary", "No result files have been created in this session."))
    tagList(h3("Generated files"), tags$ul(lapply(files, tags$li)))
  })

  output$matched_omics_gate <- renderUI({
    paths <- matched_paths()
    if (!isTRUE(state$valid)) return(tags$div(class = "alert alert-danger", "Validate an MSI dataset first."))
    if (!length(paths)) return(tags$div(class = "alert alert-warning",
      "The current dataset has no declared matched Visium/H&E companion files. Continue to optional LC-MS/MS evidence or export."))
    if (is.null(state$registration)) return(tags$div(class = "alert alert-warning", sprintf(
      "%d same-section companion files are attached to %s (%s), but MSI and Visium/H&E are still in different coordinate systems. Complete and validate registration before spatial integration.",
      length(paths), state$spec$sample_id, state$spec$section_id)))
    tags$div(class = "alert alert-success", sprintf(
      "%d companion files are attached and a registration has been fitted/validated for %s (%s).",
      length(paths), state$spec$sample_id, state$spec$section_id))
  })
  output$matched_file_inventory <- renderDT({
    paths <- matched_paths(); req(length(paths))
    datatable(data.frame(role = names(paths), file = basename(paths), bytes = file.info(paths)$size,
      stringsAsFactors = FALSE), options = list(dom = "t", scrollX = TRUE))
  })
  output$matched_histology <- renderImage({
    req(state$spec, present_path(state$spec$histology_path), file.exists(state$spec$histology_path))
    list(src = state$spec$histology_path, contentType = "image/png",
      alt = paste("Matched H&E image for", state$spec$sample_id))
  }, deleteFile = FALSE)
  output$matched_annotation_summary <- renderDT({
    req(state$spec)
    annotation_paths <- unlist(state$spec[c("st_region_path", "st_lesion_path", "st_loupe_path")], use.names = TRUE)
    annotation_paths <- annotation_paths[file.exists(annotation_paths)]
    req(length(annotation_paths))
    summary <- do.call(rbind, lapply(names(annotation_paths), function(role) {
      d <- read.csv(annotation_paths[[role]], stringsAsFactors = FALSE, check.names = FALSE)
      label_col <- setdiff(names(d), c("Barcode", "barcode"))[1]
      data.frame(annotation = role, spots = nrow(d), labels = length(unique(d[[label_col]])),
        missing_labels = sum(is.na(d[[label_col]]) | !nzchar(as.character(d[[label_col]]))),
        stringsAsFactors = FALSE)
    }))
    datatable(summary, options = list(dom = "t"))
  })

  next_button <- function(id, label) actionButton(id, label, class = "btn-success")
  output$next_from_provenance <- renderUI(if (isTRUE(state$valid)) next_button("go_processing", "Next: processing & QC") else NULL)
  output$next_from_processing <- renderUI(if (!is.null(state$processed)) next_button("go_registration", "Next: registration (optional)") else NULL)
  output$next_from_registration <- renderUI(if (!is.null(state$processed)) next_button("go_roi", "Next: optional prior ROI") else NULL)
  output$next_from_roi <- renderUI(if (!is.null(state$processed)) next_button("go_structure", "Next: detect metabolic domains") else NULL)
  output$next_from_structure <- renderUI(if (!is.null(state$domains) && !is.null(state$rois) && any(!is.na(state$rois$roi_id))) next_button("go_statistics", "Next: ROI summaries & statistics") else NULL)
  output$next_from_statistics <- renderUI(if (!is.null(state$processed)) next_button("go_validation", "Next: spatial validation") else NULL)
  output$next_from_validation <- renderUI(if (!is.null(state$processed)) next_button("go_matched", "Next: matched transcriptomics & H&E") else NULL)
  output$next_from_matched <- renderUI(if (isTRUE(state$valid)) next_button("go_lcms", "Next: optional LC-MS/MS evidence") else NULL)
  output$next_from_lcms <- renderUI(if (isTRUE(state$valid)) next_button("go_results", "Finish: results & export") else NULL)
  navigation <- c(go_processing = "spatial_processing", go_registration = "spatial_registration",
    go_roi = "spatial_roi", go_structure = "spatial_structure", go_statistics = "spatial_statistics",
    go_validation = "spatial_validation", go_matched = "matched_omics",
    go_lcms = "lcms_evidence", go_results = "results")
  lapply(names(navigation), function(id) observeEvent(input[[id]], {
    updateNavbarPage(session, "workflow_tabs", selected = navigation[[id]])
  }, ignoreInit = TRUE))

  output$validation_card <- renderUI({
    if (is.null(state$validation)) return(tags$div(class = "alert alert-secondary", "Choose inputs and validate them before continuing."))
    if (!state$validation$valid) return(tags$div(class = "alert alert-danger", strong("Input needs attention."),
      tags$ul(lapply(state$validation$errors, tags$li))))
    caps <- names(Filter(isTRUE, state$validation$capabilities))
    tags$div(class = "alert alert-success", strong("Input is ready."),
      tags$p("File pairing, format and coordinate metadata passed the available checks."),
      tags$p("Available workflow modules: ", paste(caps, collapse = ", "), "."),
      if (length(state$validation$notes)) tags$ul(lapply(state$validation$notes, tags$li)))
  })
  output$continue_processing <- renderUI(if (isTRUE(state$valid))
    actionButton("continue_button", if (isTRUE(state$validation$capabilities$msi)) "Continue to MSI processing" else "Continue to LC-MS/MS evidence",
      class = "btn-success") else NULL)
  observeEvent(input$continue_button, updateNavbarPage(session, "workflow_tabs",
    selected = if (isTRUE(state$validation$capabilities$msi)) "spatial_processing" else "lcms_evidence"))
  output$technical_files <- renderDT({ req(state$validation); files <- state$validation$files
    if (!length(files)) return(datatable(data.frame(message = "No files"), options = list(dom = "t")))
    roles <- names(files); if (is.null(roles) || !length(roles) || all(!nzchar(roles))) roles <- paste0("file_", seq_along(files))
    datatable(data.frame(role = roles, path = files,
      bytes = file.info(files)$size, md5 = unname(tools::md5sum(files))), options = list(pageLength = 5, scrollX = TRUE)) })
  output$input_metadata <- renderDT({ req(state$validation); d <- state$validation$metadata$summary
    if (is.null(d)) d <- data.frame(item = "Input type", value = state$spec$msi_input_type)
    datatable(d, options = list(dom = "t")) })
  output$module_overview <- renderUI({
    if (!isTRUE(state$valid)) return(tags$div(class = "alert alert-warning", "Validate input on the Start page first."))
    c <- state$validation$capabilities
    tags$div(class = "alert alert-info", paste(names(c), ifelse(unlist(c), "available", "not supplied"), collapse = " · "))
  })
  output$source_note <- renderText(state$example_note %or% "User-supplied dataset; attribution and study design remain the user's responsibility.")
  output$provenance_table <- renderDT({ req(state$valid); s <- state$spec
    datatable(data.frame(field = c("sample_id", "section_id", "subject_id", "ion_mode", "ion_mode_source", "input_type", "lcms_relationship"),
      value = c(s$sample_id, s$section_id, s$subject_id, s$ion_mode, s$ion_mode_source, s$msi_input_type,
        s$lcms_relationship %or% "none")), options = list(dom = "t")) })

  output$processing_action <- renderUI(if (isTRUE(state$valid) && isTRUE(state$validation$capabilities$processing))
    actionButton("run_processing", "Run processing", class = "btn-primary") else tags$div(class = "alert alert-secondary", "Processing is available after validating MSI input."))
  observeEvent(input$run_processing, {
    req(state$valid); s <- state$spec
    withProgress(message = "Preparing MSI", value = 0, {
      setProgress(value = .05, detail = "Reading the MSI input")
      if (identical(s$msi_input_type, "csv")) {
        state$processed <- read_pixel_feature_csv(s$csv_path, s$sample_id, s$section_id, s$subject_id)
      } else {
        layout <- state$validation$metadata$mz_layout
        if (identical(layout, "shared")) {
          state$processed <- load_shared_axis_pair(s$imzml_path, s$ibd_path, s$sample_id, s$section_id,
            s$ion_mode, s$ion_mode_source)
        } else {
          state$processed <- load_variable_mz_msi_features(s$imzml_path, s$ibd_path,
            s$sample_id, s$section_id, s$ion_mode, s$ion_mode_source,
            processing = s$spectral_processing, alignment_ppm = s$alignment_ppm,
            peak_pick_snr = s$peak_pick_snr, min_detection_fraction = s$min_detection_fraction,
            progress_callback = function(value, detail)
              setProgress(value = min(.45, .05 + .40 * value), detail = detail))
        }
      }
      setProgress(value = .45, detail = "Normalizing the pixel-by-feature matrix")
      feature_names <- state$processed$feature_metadata$column_name
      matrix <- as.matrix(state$processed$pixel_feature_matrix[, feature_names, drop = FALSE]); storage.mode(matrix) <- "double"
      tic <- rowSums(matrix, na.rm = TRUE); analysis <- matrix
      if (isTRUE(input$tic_normalize)) analysis <- analysis / pmax(tic, .Machine$double.eps)
      if (isTRUE(input$log1p_transform)) analysis <- log10(analysis + 1)
      state$analysis_matrix <- analysis; state$tic <- tic
      setProgress(value = .70, detail = "Preparing the tissue mask and QC summaries")
      if (isTRUE(s$make_tissue_mask)) {
        pq <- state$processed$pixel_qc
        if (is.null(pq)) pq <- data.frame(state$processed$coordinates, raw_tic = tic,
          raw_peak_count = rowSums(is.finite(matrix) & matrix > 0))
        state$tissue_gate <- build_msi_tissue_mask(state$processed$coordinates, pq$raw_tic, pq$raw_peak_count,
          "kmeans_log_tic_peak_count", seed = 20260808)
        state$tissue_mask <- state$tissue_gate$mask$tissue
        write.csv(state$tissue_gate$mask, file.path(session_dir, "tissue_mask.csv"), row.names = FALSE)
        write.csv(state$tissue_gate$diagnostics, file.path(session_dir, "tissue_mask_parameters.csv"), row.names = FALSE)
      } else state$tissue_mask <- rep(TRUE, nrow(analysis))
      setProgress(value = .85, detail = "Saving compact session metadata")
      write.csv(state$processed$qc_summary, file.path(session_dir, "qc_summary.csv"), row.names = FALSE)
      write.csv(state$processed$feature_metadata, file.path(session_dir, "feature_metadata.csv"), row.names = FALSE)
      write.csv(state$processed$coordinates, file.path(session_dir, "coordinates.csv"), row.names = FALSE)
      write.csv(state$processed$provenance, file.path(session_dir, "provenance_manifest.csv"), row.names = FALSE)
      write.csv(data.frame(tic_normalization = isTRUE(input$tic_normalize),
        log_transform = if (isTRUE(input$log1p_transform)) "log10(x + 1)" else "none",
        alignment_ppm = s$alignment_ppm, peak_pick_snr = s$peak_pick_snr,
        min_detection_fraction = s$min_detection_fraction), file.path(session_dir, "processing_parameters.csv"), row.names = FALSE)
      setProgress(value = 1, detail = "Processing complete")
    })
  })
  output$processing_gate <- renderUI(if (is.null(state$processed)) tags$div(class = "alert alert-warning", "No processed MSI is available yet.")
    else tags$div(class = "alert alert-success", "Processing completed. The full matrix remains in this R session; the browser receives previews only."))
  output$qc_table <- renderDT({ req(state$processed); datatable(state$processed$qc_summary, options = list(dom = "t", scrollX = TRUE)) })
  output$feature_preview <- renderDT({ req(state$processed); datatable(head(state$processed$feature_metadata, 25), options = list(pageLength = 10)) })
  output$tissue_mask_results <- renderUI({
    if (is.null(state$tissue_gate)) return(NULL)
    tagList(
      tags$div(class = "result-card", plotOutput("tissue_mask_plot", height = 480)),
      DTOutput("tissue_mask_diagnostics"), uiOutput("tissue_download")
    )
  })
  output$tissue_mask_plot <- renderPlot({ req(state$tissue_gate); d <- state$tissue_gate$mask
    ggplot(d, aes(x, y, color = tissue_status)) + geom_point(size = .7) + coord_equal() + scale_y_reverse() +
      scale_color_manual(values = c("tissue" = "#2166AC", "unclassified/background" = "#D9D9D9")) + theme_minimal() +
      labs(title = "Full acquisition field and reproducible tissue gate") })
  output$tissue_mask_diagnostics <- renderDT({ req(state$tissue_gate); datatable(state$tissue_gate$diagnostics, options = list(dom = "t", scrollX = TRUE)) })
  output$tissue_download <- renderUI(if (!is.null(state$tissue_gate)) downloadButton("download_tissue_mask", "Download tissue mask") else NULL)
  output$download_tissue_mask <- downloadHandler(filename = function() "tissue_mask.csv", content = function(file) {
    req(state$tissue_gate); write.csv(state$tissue_gate$mask, file, row.names = FALSE) })

  output$registration_workspace <- renderUI({
    if (!isTRUE(state$valid))
      return(tags$div(class = "alert alert-warning", "Validate a dataset before registration."))
    if (is.null(state$processed))
      return(tags$div(class = "alert alert-warning", "Process MSI before registration."))
    has_histology <- present_path(state$spec$histology_path) && file.exists(state$spec$histology_path)
    official_ready <- isTRUE(state$validation$capabilities$registration)
    methods <- c("Fit affine transform from control points" = "control_points")
    if (official_ready) methods <- c("Validate supplied official 3 × 3 transform" = "official", methods)
    controls <- tags$div(class = "registration-control-card",
      h3("Registration"),
      tags$p(if (has_histology)
        "The matched H&E is available. Choose how its coordinates should be connected to the MSI field."
        else "No matched H&E is attached. An optical image and registration information are required."),
      selectInput("registration_method", "Registration path", choices = methods,
        selected = if (official_ready) "official" else "control_points"),
      conditionalPanel("input.registration_method == 'control_points'",
        fileInput("registration_control_points", "Control-point CSV",
          accept = ".csv", placeholder = "histology_x, histology_y, msi_x, msi_y"),
        tags$p(class = "help-block",
          "Provide at least three finite, non-collinear point pairs. Optional section_id supports section-specific fits.")),
      uiOutput("registration_action"))
    if (!has_histology) return(controls)
    tags$div(class = "registration-workspace",
      tags$div(class = "registration-image-card",
        h3(paste("Matched H&E —", state$spec$sample_id)),
        tags$div(class = "responsive-image-output", imageOutput("registration_histology_preview", height = "100%"))),
      controls)
  })
  output$registration_action <- renderUI({
    req(state$processed)
    method <- input$registration_method %or% if (isTRUE(state$validation$capabilities$registration)) "official" else "control_points"
    if (identical(method, "official")) {
      if (!isTRUE(state$validation$capabilities$registration))
        return(tags$div(class = "alert alert-warning", "A compatible optical image and 3 × 3 transform JSON are required."))
      return(actionButton("run_registration", "Validate supplied transform", class = "btn-primary"))
    }
    if (!present_path(state$spec$histology_path) || !file.exists(state$spec$histology_path))
      return(tags$div(class = "alert alert-warning", "Attach an H&E/optical image before fitting control points."))
    if (is.null(input$registration_control_points))
      return(tags$div(class = "alert alert-info", "Upload a control-point CSV to enable affine fitting."))
    actionButton("run_control_registration", "Fit and validate affine transform", class = "btn-primary")
  })
  observeEvent(input$run_registration, {
    req(state$processed); s <- state$spec
    state$registration <- register_metaspace_optical(state$processed$coordinates, s$transform_path,
      s$optical_path, if (present_path(s$attribution_path)) s$attribution_path else NULL, tic = state$tic,
      representative_ions = { m <- state$analysis_matrix; ii <- head(order(apply(m, 2, var), decreasing = TRUE), 3)
        stats::setNames(lapply(ii, function(j) m[, j]), format(state$processed$feature_metadata$mz[ii], digits = 9)) },
      output_dir = file.path(session_dir, "live_registration"))
    state$registration$type <- "official_transform_validation"
    state$registration$diagnostics$diagnostic_source <- "live run using user-supplied optical image, transform and current MSI"
    write.csv(state$registration$registered_coordinates, file.path(session_dir, "registered_coordinates.csv"), row.names = FALSE)
  })
  observeEvent(input$run_control_registration, {
    req(state$processed, input$registration_control_points)
    points <- read.csv(input$registration_control_points$datapath, check.names = FALSE, stringsAsFactors = FALSE)
    require_columns(points, c("histology_x", "histology_y", "msi_x", "msi_y"), "Control-point CSV")
    section_column <- if ("section_id" %in% names(points)) "section_id" else NULL
    fit <- fit_histology_msi_registration(points, section_column = section_column)
    predicted <- transform_histology_coordinates(points, fit,
      x_column = "histology_x", y_column = "histology_y",
      output_columns = c("fitted_msi_x", "fitted_msi_y"), section_column = section_column)
    predicted$error <- sqrt((predicted$msi_x - predicted$fitted_msi_x)^2 +
      (predicted$msi_y - predicted$fitted_msi_y)^2)
    state$registration <- list(type = "control_point_affine", fit = fit,
      diagnostics = registration_diagnostics(fit))
    state$registration_points <- predicted
    write.csv(state$registration$diagnostics, file.path(session_dir, "control_point_registration_diagnostics.csv"), row.names = FALSE)
    write.csv(predicted, file.path(session_dir, "control_point_registration_fit.csv"), row.names = FALSE)
  })
  output$registration_gate <- renderUI({
    if (is.null(state$registration)) return(NULL)
    if (identical(state$registration$type, "control_point_affine"))
      return(tags$div(class = "alert alert-success",
        "Affine registration fitted from control points. Review RMSE, maximum error and residual segments before using histology-derived ROI."))
    tags$div(class = "alert alert-success",
      "The supplied official transform was applied and checked against independently derived optical and MSI masks.")
  })
  output$registration_diagnostics <- renderDT({ req(state$registration); datatable(state$registration$diagnostics, options = list(dom = "t")) })
  output$registration_results <- renderUI({
    if (is.null(state$registration)) return(NULL)
    tagList(h3("Registration diagnostics"), DTOutput("registration_diagnostics"),
      if (identical(state$registration$type, "control_point_affine"))
        plotOutput("control_registration_plot", height = 480)
      else tags$div(class = "responsive-image-output", imageOutput("registration_plot", height = "100%")))
  })
  output$registration_plot <- renderImage({ req(state$registration); req(!identical(state$registration$type, "control_point_affine")); src <- unname(state$registration$output_files["measurement_mask_overlay"])
    list(src = src, contentType = "image/png", alt = "Live registered measurement mask") }, deleteFile = FALSE)
  output$control_registration_plot <- renderPlot({
    req(state$registration_points); d <- state$registration_points
    ggplot(d) +
      geom_segment(aes(x = msi_x, y = msi_y, xend = fitted_msi_x, yend = fitted_msi_y),
        color = "#D73027", linewidth = .7, arrow = arrow(length = grid::unit(.12, "inches"))) +
      geom_point(aes(msi_x, msi_y), color = "#2166AC", size = 2.8) +
      geom_point(aes(fitted_msi_x, fitted_msi_y), color = "#D73027", shape = 4, size = 3) +
      coord_equal() + scale_y_reverse() + theme_minimal() +
      labs(title = "Observed versus affine-fitted MSI control points",
        subtitle = "Blue: observed MSI point; red cross: fitted point; segment: residual",
        x = "MSI x", y = "MSI y")
  })
  output$registration_histology_preview <- renderImage({
    req(state$spec, present_path(state$spec$histology_path), file.exists(state$spec$histology_path))
    list(src = state$spec$histology_path, contentType = "image/png",
      alt = paste("Matched H&E image for", state$spec$sample_id))
  }, deleteFile = FALSE)

  initialize_rois <- function() {
    req(state$processed)
    if (is.null(state$rois)) {
      state$rois <- data.frame(state$processed$coordinates, roi_id = NA_character_,
        stringsAsFactors = FALSE)
    }
    state$rois
  }
  invalidate_after_roi <- function() {
    state$domains <- NULL; state$domain_source <- NULL; state$domain_counts <- NULL
    state$domain_features <- NULL; state$niches <- NULL; state$sampling <- NULL
    state$metabo_input <- NULL; state$functional_peak_table <- NULL; state$functional_peak_list <- NULL; state$mapped_scores <- NULL; state$moran <- NULL
  }
  output$histology_roi_action <- renderUI({
    if (is.null(state$registration) || !identical(state$registration$type, "control_point_affine"))
      return(tags$div(class = "alert alert-warning",
        "Complete the control-point H&E→MSI registration in Step 3 before transferring an anatomical ROI."))
    if (is.null(input$histology_roi_polygon))
      return(tags$div(class = "alert alert-info", "Upload an H&E polygon CSV to continue."))
    actionButton("transfer_histology_roi", "Transfer H&E ROI to MSI", class = "btn-primary")
  })
  observeEvent(input$transfer_histology_roi, {
    req(state$processed, state$registration, input$histology_roi_polygon)
    validate(need(identical(state$registration$type, "control_point_affine"),
      "Histology ROI transfer requires a fitted H&E→MSI control-point registration."))
    vertices <- read.csv(input$histology_roi_polygon$datapath, check.names = FALSE, stringsAsFactors = FALSE)
    require_columns(vertices, c("roi_id", "x", "y"), "H&E polygon CSV")
    if (!"vertex_order" %in% names(vertices)) {
      vertices$vertex_order <- ave(seq_len(nrow(vertices)), vertices$roi_id, FUN = seq_along)
    }
    section_column <- state$registration$fit$section_column
    if (!is.null(section_column)) {
      require_columns(vertices, section_column, "H&E polygon CSV")
    }
    transformed <- transform_histology_coordinates(vertices, state$registration$fit,
      x_column = "x", y_column = "y", output_columns = c("x", "y"),
      section_column = section_column)
    selected <- select_rois(state$processed$coordinates,
      selection_mode = "manual", manual_method = "polygon",
      polygon_vertices = transformed, section_column = section_column,
      roi_column = "roi_id", overlap = "error")
    state$rois <- selected$annotated_pixels
    if (!is.null(state$tissue_mask)) state$rois$roi_id[!state$tissue_mask] <- NA_character_
    validate(need(any(!is.na(state$rois$roi_id)),
      "The transferred H&E polygon contains no current MSI tissue pixels."))
    state$rois$roi_selection_source <- "histology_polygon_transferred_after_control_point_registration"
    state$roi_source <- "histology_registered"
    invalidate_after_roi()
    write.csv(transformed, file.path(session_dir, "histology_roi_vertices_in_msi_coordinates.csv"), row.names = FALSE)
    write.csv(state$rois, file.path(session_dir, "selected_rois.csv"), row.names = FALSE)
  })
  output$histology_roi_gate <- renderUI({
    if (!identical(state$roi_source, "histology_registered")) return(NULL)
    tags$div(class = "alert alert-success",
      "The anatomical ROI was transformed from H&E coordinates into MSI coordinates using the validated control-point fit.")
  })
  observeEvent(input$add_brushed_roi, {
    req(state$processed)
    brush <- input$roi_brush
    validate(need(!is.null(brush), "Brush a rectangle on the plot first."),
      need(present_path(input$roi_name), "Provide a non-empty ROI name."))
    rois <- initialize_rois()
    selected <- rois$x >= min(brush$xmin, brush$xmax) & rois$x <= max(brush$xmin, brush$xmax) &
      rois$y >= min(brush$ymin, brush$ymax) & rois$y <= max(brush$ymin, brush$ymax)
    if (!is.null(state$tissue_mask)) selected <- selected & state$tissue_mask
    validate(need(any(selected), "The brushed rectangle contains no tissue pixels."))
    rois$roi_id[selected] <- trimws(input$roi_name)
    rois$roi_selection_source <- "direct_msi_brush"
    state$rois <- rois; state$roi_source <- "msi_coordinate"; invalidate_after_roi()
    write.csv(state$rois, file.path(session_dir, "selected_rois.csv"), row.names = FALSE)
  })
  observeEvent(input$use_full_tissue_roi, {
    rois <- initialize_rois()
    selected <- if (is.null(state$tissue_mask)) rep(TRUE, nrow(rois)) else state$tissue_mask
    rois$roi_id <- NA_character_; rois$roi_id[selected] <- "full_tissue"
    rois$roi_selection_source <- "full_msi_tissue_mask"
    state$rois <- rois; state$roi_source <- "msi_coordinate"; invalidate_after_roi()
    write.csv(state$rois, file.path(session_dir, "selected_rois.csv"), row.names = FALSE)
  })
  observeEvent(input$clear_rois, {
    req(state$processed); state$rois <- NULL; state$roi_source <- NULL; invalidate_after_roi()
  })
  observeEvent(input$roi_table_upload, {
    req(state$processed)
    table <- read.csv(input$roi_table_upload$datapath, check.names = FALSE, stringsAsFactors = FALSE)
    pixel <- data.frame(state$processed$coordinates, check.names = FALSE)
    state$rois <- define_rois(pixel, mode = "coordinate", roi_table = table)
    if (!is.null(state$tissue_mask)) state$rois$roi_id[!state$tissue_mask] <- NA_character_
    validate(need(any(!is.na(state$rois$roi_id)), "The ROI table did not select any tissue pixels."))
    state$rois$roi_selection_source <- "imported_msi_coordinate_table"
    state$roi_source <- "msi_coordinate"; invalidate_after_roi()
    write.csv(state$rois, file.path(session_dir, "selected_rois.csv"), row.names = FALSE)
  })
  output$roi_selection_plot <- renderPlot({
    req(state$processed)
    d <- state$processed$coordinates
    d$roi_id <- if (is.null(state$rois)) NA_character_ else state$rois$roi_id
    d$tissue <- if (is.null(state$tissue_mask)) TRUE else state$tissue_mask
    ggplot(d, aes(x, y)) +
      geom_point(data = d[!d$tissue, , drop = FALSE], color = "grey90", size = .55) +
      geom_point(data = d[d$tissue, , drop = FALSE], aes(color = roi_id), size = .7) +
      coord_equal() + scale_y_reverse() + theme_minimal() +
      labs(title = if (identical(input$roi_path, "histology"))
        "MSI pixels receiving the transferred anatomical ROI" else if (identical(input$roi_path, "msi"))
          "Brush the MSI field to define a known prior ROI" else "Complete tissue will enter domain detection",
        color = "Selected ROI")
  })
  output$roi_gate <- renderUI({
    if (is.null(state$processed)) return(tags$div(class = "alert alert-warning", "Process MSI before selecting ROI."))
    selected <- if (is.null(state$rois)) 0L else sum(!is.na(state$rois$roi_id) & nzchar(state$rois$roi_id))
    if (!selected) tags$div(class = "alert alert-info", "No prior ROI selected. Domain detection will use the complete tissue mask; you can select domains as ROI afterward.")
    else tags$div(class = "alert alert-success", sprintf(
      "Prior ROI contains %d pixels (source: %s). Domain detection will be restricted to this field.",
      selected, state$roi_source %or% "unspecified"))
  })
  output$roi_summary <- renderDT({
    req(state$rois); d <- state$rois[!is.na(state$rois$roi_id), , drop = FALSE]; req(nrow(d))
    counts <- as.data.frame(table(d$roi_id), stringsAsFactors = FALSE); names(counts) <- c("roi_id", "pixel_count")
    datatable(counts, options = list(dom = "t"))
  })
  output$roi_download <- renderUI(if (!is.null(state$rois)) downloadButton("download_rois", "Download selected ROI pixels") else NULL)
  output$download_rois <- downloadHandler(filename = function() "selected_rois.csv", content = function(file) {
    req(state$rois); write.csv(state$rois, file, row.names = FALSE)
  })

  domain_analysis_keep <- function() {
    req(state$processed)
    keep <- if (is.null(state$tissue_mask)) rep(TRUE, nrow(state$processed$coordinates)) else state$tissue_mask
    if (!is.null(state$rois) && any(!is.na(state$rois$roi_id) & nzchar(state$rois$roi_id))) {
      keep <- keep & !is.na(state$rois$roi_id) & nzchar(state$rois$roi_id)
    }
    keep
  }
  restrict_domains_to_field <- function(domains) {
    keep <- domain_analysis_keep()
    domains$domain_id[!keep] <- "-1"
    domains$domain_label[!keep] <- "unclassified/background"
    domains
  }
  observeEvent(input$domain_csv_runtime, {
    req(state$processed)
    info <- validate_label_file(input$domain_csv_runtime$datapath)
    state$domains <- restrict_domains_to_field(resolve_label_mapping(
      info, state$processed$coordinates, state$spec$sample_id, state$spec$section_id))
    state$domain_source <- "external"
    state$niches <- NULL
    state$sampling <- NULL; state$metabo_input <- NULL; state$functional_peak_table <- NULL; state$functional_peak_list <- NULL; state$mapped_scores <- NULL
    state$domain_counts <- as.data.frame(table(state$domains$domain_label)); names(state$domain_counts) <- c("domain", "pixel_count")
  })
  output$domain_action <- renderUI({
    if (!isTRUE(state$valid) || !isTRUE(state$validation$capabilities$domains))
      return(tags$div(class = "alert alert-secondary", "Domain analysis is disabled until MSI input is validated."))
    if (is.null(state$processed)) return(tags$div(class = "alert alert-warning", "Process MSI first."))
    tagList(actionButton("load_spec_labels", "Use validated ROI/domain labels"),
      actionButton("generate_domains", "Generate exploratory metabolic domains", class = "btn-primary"))
  })
  observeEvent(input$load_spec_labels, {
    validate(need(isTRUE(state$validation$capabilities$labels), "No validated ROI/domain CSV was supplied."))
    state$domains <- restrict_domains_to_field(resolve_label_mapping(
      state$validation$labels, state$processed$coordinates, state$spec$sample_id, state$spec$section_id))
    state$domain_source <- "external"
    state$niches <- NULL
    state$sampling <- NULL; state$metabo_input <- NULL; state$functional_peak_table <- NULL; state$functional_peak_list <- NULL; state$mapped_scores <- NULL
    state$domain_counts <- as.data.frame(table(state$domains$domain_label)); names(state$domain_counts) <- c("domain", "pixel_count")
  })
  observeEvent(input$generate_domains, {
    req(state$processed, state$analysis_matrix)
    keep_tissue <- domain_analysis_keep()
    feature_index <- seq_len(ncol(state$analysis_matrix))
    if (!identical(input$domain_feature_source, "all")) {
      validate(need(!is.null(state$moran), "Compute Moran's I in Processing & QC before using Moran-screened features."))
      selected_features <- if (identical(input$domain_feature_source, "moran_fdr")) {
        state$moran$feature[is.finite(state$moran$adj_p_value) & state$moran$adj_p_value <= input$domain_moran_fdr]
      } else head(state$moran$feature, input$domain_moran_top_n)
      feature_index <- match(selected_features, state$processed$feature_metadata$column_name)
      feature_index <- unique(feature_index[!is.na(feature_index)])
      validate(need(length(feature_index) >= 2L, "The selected Moran filter retained fewer than two available features."))
    }
    matrix <- state$analysis_matrix[keep_tissue, feature_index, drop = FALSE]; keep <- apply(matrix, 2, var) > 0
    validate(need(sum(keep) >= 2, "At least two non-zero-variance features are required."))
    pc <- prcomp(matrix[, keep, drop = FALSE], center = TRUE, scale. = TRUE); npc <- min(input$domain_pcs, ncol(pc$x))
    set.seed(input$domain_seed); km <- kmeans(pc$x[, seq_len(npc), drop = FALSE], centers = input$domain_k, nstart = 25)
    ids <- rep(-1L, nrow(state$analysis_matrix)); ids[keep_tissue] <- km$cluster - 1L
    state$domains <- data.frame(state$processed$coordinates, tissue_status = ifelse(keep_tissue, "tissue", "unclassified/background"),
      domain_id = ids, domain_label = ifelse(ids == -1L, "unclassified/background", paste0("data-driven exploratory metabolic domain ", ids)))
    state$domain_source <- "data_driven"
    state$niches <- NULL
    state$sampling <- NULL; state$metabo_input <- NULL; state$functional_peak_table <- NULL; state$functional_peak_list <- NULL; state$mapped_scores <- NULL
    state$domain_counts <- as.data.frame(table(state$domains$domain_label)); names(state$domain_counts) <- c("domain", "pixel_count")
    between <- vapply(seq_len(ncol(matrix)), function(j) var(tapply(matrix[, j], km$cluster, mean)), numeric(1)); order <- head(order(between, decreasing = TRUE), 25)
    retained_metadata <- state$processed$feature_metadata[feature_index, , drop = FALSE]
    state$domain_features <- data.frame(feature = retained_metadata$column_name[order],
      mz = retained_metadata$mz[order], between_domain_variance = between[order])
    write.csv(state$domains, file.path(session_dir, "exploratory_metabolic_domains.csv"), row.names = FALSE)
    write.csv(data.frame(k = input$domain_k, seed = input$domain_seed, n_pcs = npc,
      feature_source = input$domain_feature_source, features_used = sum(keep),
      inference = "descriptive_only_no_pixel_replicates"), file.path(session_dir, "domain_parameters.csv"), row.names = FALSE)
  })
  output$domain_gate <- renderUI(if (is.null(state$domains)) tags$div(class = "alert alert-warning", "No labels or generated domains are currently loaded.") else NULL)
  output$domain_plot <- renderPlot({ req(state$domains); ggplot(state$domains, aes(x, y, color = domain_label)) + geom_point(size = .7) + coord_equal() + scale_y_reverse() + theme_minimal() })
  output$domain_counts <- renderDT({ req(state$domain_counts); datatable(state$domain_counts, options = list(dom = "t")) })
  output$domain_features <- renderDT({ req(state$domain_features); datatable(state$domain_features, options = list(pageLength = 10)) })
  output$domain_download <- renderUI(if (!is.null(state$domains)) downloadButton("download_domains", "Download domain CSV") else NULL)
  output$download_domains <- downloadHandler(filename = function() "spatial_domains.csv", content = function(file) { req(state$domains); write.csv(state$domains, file, row.names = FALSE) })

  comparison_feature_table <- function(labels, selected, label_column) {
    validate(need(length(selected) == 2L, "Select exactly two groups to compare."))
    idx_a <- which(as.character(labels[[label_column]]) == selected[1])
    idx_b <- which(as.character(labels[[label_column]]) == selected[2])
    validate(need(length(idx_a) > 0L && length(idx_b) > 0L, "Both selected groups must contain pixels."))
    matrix <- state$analysis_matrix
    mean_a <- colMeans(matrix[idx_a, , drop = FALSE], na.rm = TRUE)
    mean_b <- colMeans(matrix[idx_b, , drop = FALSE], na.rm = TRUE)
    sd_a <- apply(matrix[idx_a, , drop = FALSE], 2, sd, na.rm = TRUE)
    sd_b <- apply(matrix[idx_b, , drop = FALSE], 2, sd, na.rm = TRUE)
    pooled <- sqrt((sd_a^2 + sd_b^2) / 2)
    standardized <- (mean_a - mean_b) / pooled
    standardized[!is.finite(standardized)] <- 0
    out <- data.frame(
      feature = state$processed$feature_metadata$column_name,
      mz = state$processed$feature_metadata$mz,
      group_a = selected[1], group_b = selected[2],
      mean_a = mean_a, mean_b = mean_b,
      difference_a_minus_b = mean_a - mean_b,
      standardized_difference = standardized,
      stringsAsFactors = FALSE)
    out[order(abs(out$standardized_difference), decreasing = TRUE), , drop = FALSE]
  }
  comparison_contact_summary <- function(labels, selected, label_column) {
    validate(need(length(selected) == 2L, "Select exactly two groups to compare."))
    x <- as.numeric(labels$x); y <- as.numeric(labels$y)
    positive_step <- function(z) {
      d <- diff(sort(unique(z))); d <- d[is.finite(d) & d > 0]
      if (length(d)) min(d) else 1
    }
    gx <- round((x - min(x)) / positive_step(x))
    gy <- round((y - min(y)) / positive_step(y))
    graph <- build_spatial_neighbors(gx, gy, method = "queen")
    edges <- graph$undirected_edges
    group <- as.character(labels[[label_column]])
    cross <- if (nrow(edges)) {
      (group[edges$from] == selected[1] & group[edges$to] == selected[2]) |
        (group[edges$from] == selected[2] & group[edges$to] == selected[1])
    } else logical()
    touching <- if (nrow(edges)) group[edges$from] %in% selected & group[edges$to] %in% selected else logical()
    data.frame(group_a = selected[1], group_b = selected[2],
      pixels_a = sum(group == selected[1]), pixels_b = sum(group == selected[2]),
      direct_boundary_edges = sum(cross),
      boundary_fraction_among_selected_edges = if (sum(touching)) sum(cross) / sum(touching) else 0,
      stringsAsFactors = FALSE)
  }
  domain_comparison <- reactive({
    req(state$domains, input$domain_compare_labels)
    selected <- input$domain_compare_labels
    list(selected = selected,
      features = comparison_feature_table(state$domains, selected, "domain_label"),
      contact = comparison_contact_summary(state$domains, selected, "domain_label"))
  })
  output$domain_compare_selector <- renderUI({
    req(state$domains)
    labels <- sort(setdiff(unique(as.character(state$domains$domain_label)), "unclassified/background"))
    tags$div(class = "well", h4("Explore two domains"),
      tags$p("Choose two domains repeatedly without rerunning clustering. The map, boundary contact and descriptive molecular contrast update automatically."),
      checkboxGroupInput("domain_compare_labels", "Select exactly two domains", choices = labels))
  })
  output$domain_compare_summary <- renderUI({
    result <- domain_comparison(); d <- result$contact
    tags$div(class = "alert alert-info", sprintf(
      "%s: %d pixels; %s: %d pixels. Direct shared-boundary edges: %d.",
      d$group_a, d$pixels_a, d$group_b, d$pixels_b, d$direct_boundary_edges))
  })
  output$domain_compare_plot <- renderPlot({
    result <- domain_comparison(); d <- state$domains
    d$comparison <- ifelse(d$domain_label %in% result$selected, d$domain_label, "other/background")
    ggplot(d, aes(x, y, color = comparison)) + geom_point(size = .75) + coord_equal() +
      scale_y_reverse() + theme_minimal() + labs(title = "Selected domain pair", color = NULL)
  })
  output$domain_compare_features <- renderDT({
    result <- domain_comparison()
    datatable(head(result$features, 100), options = list(pageLength = 10, scrollX = TRUE))
  })

  output$domain_roi_selector <- renderUI({
    req(state$domains)
    labels <- sort(unique(as.character(state$domains$domain_label)))
    labels <- labels[!is.na(labels) & nzchar(labels) & labels != "unclassified/background"]
    req(length(labels))
    tags$div(class = "well",
      h4("Data-driven ROI from detected domains"),
      tags$p("Select one or more domains after reviewing the map. This creates the ROI used for subregion summaries and downstream statistics; it does not refit the domains."),
      checkboxGroupInput("domain_roi_labels", "Domains to retain", choices = labels),
      radioButtons("domain_roi_mode", "When several domains are selected", inline = TRUE,
        choices = c("Keep as separate ROI groups" = "separate", "Combine into one ROI" = "combined"),
        selected = "separate"),
      textInput("combined_domain_roi_name", "Combined ROI name", "selected_domains"),
      actionButton("use_domains_as_roi", "Use selected domains as downstream ROI", class = "btn-success"))
  })
  observeEvent(input$use_domains_as_roi, {
    req(state$processed, state$domains)
    selected_labels <- input$domain_roi_labels
    validate(need(length(selected_labels) > 0L, "Select at least one non-background domain."))
    selected <- as.character(state$domains$domain_label) %in% selected_labels
    validate(need(any(selected), "The selected domains contain no pixels."))
    roi_id <- rep(NA_character_, nrow(state$domains))
    if (identical(input$domain_roi_mode, "combined")) {
      validate(need(present_path(input$combined_domain_roi_name), "Provide a non-empty combined ROI name."))
      roi_id[selected] <- trimws(input$combined_domain_roi_name)
    } else {
      roi_id[selected] <- as.character(state$domains$domain_label[selected])
    }
    state$rois <- data.frame(state$processed$coordinates, roi_id = roi_id,
      roi_selection_source = "selected_data_driven_metabolic_domains",
      stringsAsFactors = FALSE)
    state$roi_source <- "data_driven_domains"
    state$sampling <- NULL; state$metabo_input <- NULL
    state$functional_peak_table <- NULL; state$functional_peak_list <- NULL
    state$mapped_scores <- NULL; state$moran <- NULL
    write.csv(state$rois, file.path(session_dir, "selected_rois.csv"), row.names = FALSE)
  })
  output$domain_roi_gate <- renderUI({
    if (!identical(state$roi_source, "data_driven_domains")) return(NULL)
    n <- sum(!is.na(state$rois$roi_id) & nzchar(state$rois$roi_id))
    tags$div(class = "alert alert-success", sprintf(
      "The selected metabolic domain(s) now define the downstream ROI (%d pixels). Domain and optional niche results were preserved.", n))
  })

  output$niche_action <- renderUI({
    if (is.null(state$domains)) return(tags$div(class = "alert alert-warning", "Load or generate spatial domains first."))
    actionButton("run_niches", "Detect spatial niches", class = "btn-primary")
  })
  observeEvent(input$run_niches, {
    req(state$domains)
    domain_id <- as.character(state$domains$domain_id)
    keep <- !is.na(domain_id) & domain_id != "-1" &
      !is.na(state$domains$domain_label) & state$domains$domain_label != "unclassified/background"
    validate(need(sum(keep) >= input$niche_count, "Too few non-background domain pixels for the requested niche count."))
    niche_input <- state$domains[keep, , drop = FALSE]
    subject_column <- if ("subject_id" %in% names(niche_input)) "subject_id" else NULL
    section_column <- if ("section_id" %in% names(niche_input)) "section_id" else NULL
    state$niches <- withProgress(message = "Detecting spatial niches", value = .1, {
      result <- detect_spatial_niches(
        niche_input, domain_column = "domain_label", x_col = "x", y_col = "y",
        subject_column = subject_column, section_column = section_column,
        neighborhood = input$niche_neighborhood,
        radius = if (identical(input$niche_neighborhood, "radius")) input$niche_radius else NULL,
        k = input$niche_k, include_self = isTRUE(input$niche_include_self),
        min_neighbors = input$niche_min_neighbors,
        transform = input$niche_transform, k_niches = input$niche_count,
        domain_alignment = input$niche_domain_alignment, seed = input$niche_seed
      )
      incProgress(.9)
      result
    })
    write.csv(state$niches$matrix, file.path(session_dir, "spatial_niches.csv"), row.names = FALSE)
    write.csv(state$niches$exclusion_summary, file.path(session_dir, "spatial_niche_exclusion_summary.csv"), row.names = FALSE)
    write.csv(state$niches$centers, file.path(session_dir, "spatial_niche_transformed_centers.csv"), row.names = TRUE)
    write.csv(state$niches$composition_centers, file.path(session_dir, "spatial_niche_composition_centers.csv"), row.names = TRUE)
    write.csv(data.frame(
      neighborhood = input$niche_neighborhood,
      radius = if (identical(input$niche_neighborhood, "radius")) input$niche_radius else NA_real_,
      k = if (identical(input$niche_neighborhood, "knn")) input$niche_k else NA_integer_,
      min_neighbors = input$niche_min_neighbors, include_self = isTRUE(input$niche_include_self),
      transform = input$niche_transform, k_niches = input$niche_count,
      domain_alignment = input$niche_domain_alignment, seed = input$niche_seed,
      stringsAsFactors = FALSE
    ), file.path(session_dir, "spatial_niche_parameters.csv"), row.names = FALSE)
  })
  output$niche_gate <- renderUI({
    if (is.null(state$niches)) return(tags$div(class = "alert alert-secondary", "No niche result is available."))
    d <- state$niches$matrix
    tags$div(class = "alert alert-success", sprintf(
      "Assigned %d of %d domain pixels; %d were excluded for insufficient neighborhood support.",
      sum(d$eligible), nrow(d), sum(!d$eligible)
    ))
  })
  output$niche_plot <- renderPlot({
    req(state$niches); d <- state$niches$matrix
    ggplot(d, aes(x, y, color = niche_label)) + geom_point(size = .7) + coord_equal() +
      scale_y_reverse() + theme_minimal() +
      labs(title = "Local metabolic-domain composition niches", color = "Niche")
  })
  output$niche_ambiguity_plot <- renderPlot({
    req(state$niches); d <- state$niches$matrix
    ggplot(d, aes(x, y, color = ambiguity_ratio)) + geom_point(size = .7) + coord_equal() +
      scale_y_reverse() + scale_color_viridis_c(limits = c(0, 1), na.value = "#D9D9D9") + theme_minimal() +
      labs(title = "Niche assignment ambiguity", color = "Nearest / second")
  })
  niche_labels_full <- reactive({
    req(state$niches, state$processed)
    niche <- state$niches$matrix
    full <- data.frame(state$processed$coordinates, niche_label = NA_character_, stringsAsFactors = FALSE)
    if ("pixel_id" %in% names(niche) && "pixel_id" %in% names(full)) {
      index <- match(as.character(niche$pixel_id), as.character(full$pixel_id))
    } else {
      full_key <- paste(full$sample_id, full$section_id, full$x, full$y, sep = "\r")
      niche_key <- paste(niche$sample_id, niche$section_id, niche$x, niche$y, sep = "\r")
      index <- match(niche_key, full_key)
    }
    valid <- !is.na(index)
    full$niche_label[index[valid]] <- as.character(niche$niche_label[valid])
    full
  })
  niche_comparison <- reactive({
    req(input$niche_compare_labels)
    selected <- input$niche_compare_labels
    labels <- niche_labels_full()
    composition <- state$niches$composition_centers
    ids <- suppressWarnings(as.integer(sub("^niche_", "", selected)))
    composition_pair <- as.data.frame(composition[ids, , drop = FALSE], check.names = FALSE)
    composition_pair$niche <- selected
    composition_pair <- composition_pair[, c("niche", setdiff(names(composition_pair), "niche")), drop = FALSE]
    list(selected = selected,
      features = comparison_feature_table(labels, selected, "niche_label"),
      contact = comparison_contact_summary(labels, selected, "niche_label"),
      composition = composition_pair)
  })
  output$niche_compare_selector <- renderUI({
    req(state$niches)
    labels <- sort(unique(na.omit(as.character(state$niches$matrix$niche_label))))
    tags$div(class = "well", h4("Explore two niches"),
      tags$p("Choose two niches repeatedly. The comparison shows where they occur, whether they share a boundary, how their surrounding-domain compositions differ, and their descriptive molecular contrast."),
      checkboxGroupInput("niche_compare_labels", "Select exactly two niches", choices = labels))
  })
  output$niche_compare_summary <- renderUI({
    result <- niche_comparison(); d <- result$contact
    tags$div(class = "alert alert-info", sprintf(
      "%s: %d pixels; %s: %d pixels. Direct shared-boundary edges: %d. This describes spatial association, not cell–cell interaction.",
      d$group_a, d$pixels_a, d$group_b, d$pixels_b, d$direct_boundary_edges))
  })
  output$niche_compare_plot <- renderPlot({
    result <- niche_comparison(); d <- niche_labels_full()
    d$comparison <- ifelse(d$niche_label %in% result$selected, d$niche_label, "other/background")
    ggplot(d, aes(x, y, color = comparison)) + geom_point(size = .75) + coord_equal() +
      scale_y_reverse() + theme_minimal() + labs(title = "Selected niche pair", color = NULL)
  })
  output$niche_compare_composition <- renderDT({
    result <- niche_comparison()
    datatable(result$composition, options = list(dom = "t", scrollX = TRUE))
  })
  output$niche_compare_features <- renderDT({
    result <- niche_comparison()
    datatable(head(result$features, 100), options = list(pageLength = 10, scrollX = TRUE))
  })
  output$niche_exclusion_summary <- renderDT({
    req(state$niches); datatable(state$niches$exclusion_summary, options = list(dom = "t", scrollX = TRUE))
  })
  output$niche_download <- renderUI(if (!is.null(state$niches)) tagList(
    downloadButton("download_niches", "Download niche pixels"),
    downloadButton("download_niche_summary", "Download exclusion summary")
  ) else NULL)
  output$download_niches <- downloadHandler(filename = function() "spatial_niches.csv",
    content = function(file) { req(state$niches); write.csv(state$niches$matrix, file, row.names = FALSE) })
  output$download_niche_summary <- downloadHandler(filename = function() "spatial_niche_exclusion_summary.csv",
    content = function(file) { req(state$niches); write.csv(state$niches$exclusion_summary, file, row.names = FALSE) })

  analysis_pixel_matrix <- function() {
    req(state$processed, state$analysis_matrix)
    out <- data.frame(state$processed$coordinates, state$analysis_matrix, check.names = FALSE)
    names(out)[(ncol(out) - ncol(state$analysis_matrix) + 1L):ncol(out)] <- state$processed$feature_metadata$column_name
    out
  }
  make_functional_peak_table <- function(sample_matrix) {
    features <- grep("^mz_", names(sample_matrix), value = TRUE)
    validate(need(length(features) > 0L, "No m/z features are available for functional export."))
    mz <- suppressWarnings(as.numeric(sub("^mz_", "", sub("__.*$", "", features))))
    validate(need(all(is.finite(mz)), "Every functional-analysis feature must have a numeric m/z."))
    values <- t(as.matrix(sample_matrix[, features, drop = FALSE]))
    storage.mode(values) <- "double"
    table <- data.frame(feature = c("Sample", "Group", format(mz, scientific = FALSE, trim = TRUE)),
      stringsAsFactors = FALSE, check.names = FALSE)
    for (i in seq_len(nrow(sample_matrix))) {
      table[[as.character(sample_matrix$sample_id[i])]] <- c(
        as.character(sample_matrix$sample_id[i]), as.character(sample_matrix$roi_id[i]),
        format(values[, i], scientific = FALSE, trim = TRUE))
    }
    table
  }
  prepare_lipid_name_export <- function(annotation, peak_statistics, ppm = 5, p_cutoff = .05) {
    if (!is.data.frame(annotation) || !nrow(annotation))
      stop("The lipid annotation table is empty.", call. = FALSE)
    require_columns(peak_statistics, c("m.z", "p.value", "t.score"), "Differential peak statistics")
    first_column <- function(candidates, label, required = TRUE) {
      hit <- intersect(candidates, names(annotation))
      if (length(hit)) return(hit[1])
      if (required) stop("Lipid annotation CSV needs a ", label, " column. Accepted names: ",
        paste(candidates, collapse = ", "), ".", call. = FALSE)
      NA_character_
    }
    name_col <- first_column(c("lipid_name", "LipidName", "lipid", "Name", "name",
      "metabolite_name", "moleculeName", "moleculeNames"), "lipid name")
    mz_col <- first_column(c("mz", "m.z", "MZ", "Input Mass", "input_mass",
      "experimental_mz", "theoretical_mz"), "m/z", required = FALSE)
    feature_col <- first_column(c("column_name", "feature", "Feature", "msi_column_name"),
      "MSI feature", required = FALSE)
    id_col <- first_column(c("lipid_id", "LipidID", "LipidMaps_ID", "LIPID_MAPS_ID",
      "SwissLipids_ID", "HMDB_ID", "hmdb_id", "KEGG_ID", "moleculeIds"),
      "lipid identifier", required = FALSE)
    if (is.na(mz_col) && is.na(feature_col))
      stop("Lipid annotation CSV needs either an m/z column or an mz_* feature column.", call. = FALSE)
    annotation_mz <- if (!is.na(mz_col)) suppressWarnings(as.numeric(annotation[[mz_col]])) else {
      suppressWarnings(as.numeric(sub("__.*$", "", sub("^mz_", "", annotation[[feature_col]]))))
    }
    if (all(!is.finite(annotation_mz)))
      stop("No finite m/z values could be read from the lipid annotation CSV.", call. = FALSE)
    peak_mz <- suppressWarnings(as.numeric(peak_statistics$m.z))
    nearest <- vapply(annotation_mz, function(value) {
      if (!is.finite(value)) return(NA_integer_)
      which.min(abs(peak_mz - value))
    }, integer(1))
    matched_mz <- rep(NA_real_, length(nearest)); valid_index <- !is.na(nearest)
    matched_mz[valid_index] <- peak_mz[nearest[valid_index]]
    ppm_error <- abs(matched_mz - annotation_mz) / annotation_mz * 1e6
    matched <- is.finite(ppm_error) & ppm_error <= ppm
    raw_names <- trimws(as.character(annotation[[name_col]]))
    rows <- lapply(which(matched & !is.na(raw_names) & nzchar(raw_names)), function(i) {
      candidate_names <- trimws(unlist(strsplit(raw_names[i], "[;|]")))
      candidate_names <- unique(candidate_names[nzchar(candidate_names)])
      if (!length(candidate_names)) return(NULL)
      j <- nearest[i]
      data.frame(
        lipid_name = candidate_names,
        lipid_id = if (!is.na(id_col)) as.character(annotation[[id_col]][i]) else NA_character_,
        annotation_mz = annotation_mz[i], msi_mz = peak_mz[j], ppm_error = ppm_error[i],
        msi_feature = paste0("mz_", format(peak_mz[j], scientific = FALSE, trim = TRUE)),
        p_value = peak_statistics$p.value[j], effect_score = peak_statistics$t.score[j],
        selected_for_enrichment = is.finite(peak_statistics$p.value[j]) && peak_statistics$p.value[j] <= p_cutoff,
        mapping_evidence = "user-supplied lipid name mapped to MSI m/z within ppm tolerance",
        stringsAsFactors = FALSE
      )
    })
    rows <- Filter(Negate(is.null), rows)
    if (!length(rows)) stop("No lipid annotation matched the differential MSI features within the ppm tolerance.", call. = FALSE)
    audit <- unique(do.call(rbind, rows))
    audit <- audit[order(audit$p_value, audit$ppm_error, audit$lipid_name), , drop = FALSE]
    names <- unique(audit$lipid_name[audit$selected_for_enrichment])
    list(audit = audit, names = names, settings = data.frame(ppm = ppm, p_cutoff = p_cutoff))
  }
  output$sampling_action <- renderUI({
    if (is.null(state$rois) || !any(!is.na(state$rois$roi_id)))
      return(tags$div(class = "alert alert-warning", "Select ROI pixels first."))
    actionButton("run_sampling", "Create ROI subregion samples", class = "btn-primary")
  })
  observeEvent(input$run_sampling, {
    req(state$rois)
    pixels <- analysis_pixel_matrix()
    validate(need(nrow(pixels) == nrow(state$rois), "ROI and processed-pixel rows are not aligned."))
    pixels$roi_id <- as.character(state$rois$roi_id)
    state$sampling <- sample_subregions(
      pixels, grid_size = input$sampling_grid_size,
      min_pixels = input$sampling_min_pixels, roi_column = "roi_id",
      grid_scope = input$sampling_grid_scope
    )
    validate(need(nrow(state$sampling$sample_matrix) > 0L,
      "No subregion met the minimum-pixel requirement."))
    state$metabo_input <- make_metaboanalyst_data(state$sampling$sample_matrix, group_column = "roi_id")
    state$functional_peak_table <- make_functional_peak_table(state$sampling$sample_matrix)
    roi_groups <- unique(as.character(state$sampling$sample_matrix$roi_id))
    roi_groups <- roi_groups[!is.na(roi_groups) & nzchar(roi_groups)]
    state$functional_peak_list <- NULL
    state$lipid_export <- NULL; state$lipid_name_list <- NULL
    if (length(roi_groups) == 2L) {
      differential <- suppressWarnings(differential_region_analysis(
        state$sampling$sample_matrix, group_column = "roi_id"))
      state$functional_peak_list <- data.frame(
        m.z = differential$mz,
        p.value = ifelse(is.finite(differential$p_value), pmax(differential$p_value, .Machine$double.xmin), 1),
        t.score = ifelse(is.finite(differential$effect_size), differential$effect_size, 0),
        check.names = FALSE
      )
      state$functional_peak_list <- state$functional_peak_list[
        is.finite(state$functional_peak_list$m.z), , drop = FALSE]
    }
    state$mapped_scores <- NULL
    write.csv(state$sampling$sample_matrix, file.path(session_dir, "roi_subregion_sample_matrix.csv"), row.names = FALSE)
    write.csv(state$sampling$sample_mapping, file.path(session_dir, "roi_subregion_sample_mapping.csv"), row.names = FALSE)
    write.csv(state$metabo_input, file.path(session_dir, "metaboanalyst_statistical_analysis_samples_in_rows.csv"), row.names = FALSE)
    write.table(state$functional_peak_table,
      file.path(session_dir, "metaboanalyst_functional_peak_intensity_table.txt"),
      sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
    if (!is.null(state$functional_peak_list)) write.table(state$functional_peak_list,
      file.path(session_dir, "metaboanalyst_functional_peak_list_mz_pvalue_score.txt"),
      sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
    write.csv(data.frame(
      grid_size = input$sampling_grid_size, min_pixels = input$sampling_min_pixels,
      grid_scope = input$sampling_grid_scope,
      inference_warning = "subregions_are_not_independent_biological_replicates",
      stringsAsFactors = FALSE
    ), file.path(session_dir, "roi_sampling_parameters.csv"), row.names = FALSE)
  })
  output$sampling_gate <- renderUI({
    if (is.null(state$sampling)) return(tags$div(class = "alert alert-secondary", "No ROI samples are available."))
    tags$div(class = "alert alert-success", sprintf(
      "Created %d subregion samples across %d ROI/domain labels.",
      nrow(state$sampling$sample_matrix), length(unique(state$sampling$sample_matrix$roi_id))
    ))
  })
  output$sampling_results_ui <- renderUI(if (!is.null(state$sampling)) tagList(
    uiOutput("sampling_download"), uiOutput("lipid_export_ui"),
    tags$details(class = "result-card", tags$summary("Preview generated subregion tables"),
      h4("Sample matrix"), DTOutput("sample_matrix_preview"),
      h4("Pixel-to-sample mapping"), DTOutput("sample_mapping_preview"))
  ) else NULL)
  output$sample_matrix_preview <- renderDT({ req(state$sampling)
    datatable(head(state$sampling$sample_matrix, 100), options = list(pageLength = 10, scrollX = TRUE)) })
  output$sample_mapping_preview <- renderDT({ req(state$sampling)
    datatable(head(state$sampling$sample_mapping, 100), options = list(pageLength = 10, scrollX = TRUE)) })
  output$sampling_download <- renderUI(if (!is.null(state$sampling)) tagList(
    tags$div(class = "alert alert-info",
      tags$p(strong("For Statistical Analysis → one-factor:"),
        " upload the samples-in-rows CSV. Sample is column 1 and Group is column 2."),
      tags$p(strong("For Functional Analysis:"),
        " use the peak-intensity TXT under ‘A peak intensity table’, or use the m.z/p.value/t.score TXT under ‘A peak list profile’. Functional analysis requires exactly two ROI groups for the generated differential peak list.")),
    downloadButton("download_metabo_input", "Statistical Analysis table (.csv)"),
    if (length(unique(state$sampling$sample_matrix$roi_id)) >= 2L)
      downloadButton("download_functional_peak_table", "Functional peak-intensity table (.txt)")
    else tags$span(class = "text-warning",
      " Functional export disabled: define at least two named ROI groups."),
    if (!is.null(state$functional_peak_list))
      downloadButton("download_functional_peak_list", "Functional peak list: m.z + statistics (.txt)"),
    downloadButton("download_sample_mapping", "Download sample mapping")
  ) else NULL)
  output$download_metabo_input <- downloadHandler(filename = function() "metaboanalyst_statistical_analysis_samples_in_rows.csv",
    content = function(file) { req(state$metabo_input); write.csv(state$metabo_input, file, row.names = FALSE) })
  output$download_functional_peak_table <- downloadHandler(
    filename = function() "metaboanalyst_functional_peak_intensity_table.txt",
    content = function(file) {
      req(state$functional_peak_table)
      write.table(state$functional_peak_table, file, sep = "\t", row.names = FALSE,
        col.names = FALSE, quote = FALSE)
    })
  output$download_functional_peak_list <- downloadHandler(
    filename = function() "metaboanalyst_functional_peak_list_mz_pvalue_score.txt",
    content = function(file) {
      req(state$functional_peak_list)
      write.table(state$functional_peak_list, file, sep = "\t", row.names = FALSE,
        col.names = TRUE, quote = FALSE)
    })
  output$download_sample_mapping <- downloadHandler(filename = function() "roi_subregion_sample_mapping.csv",
    content = function(file) { req(state$sampling); write.csv(state$sampling$sample_mapping, file, row.names = FALSE) })

  output$lipid_export_ui <- renderUI({
    req(state$sampling)
    tags$div(class = "result-card",
      h4("MetaboAnalyst lipid enrichment export"),
      tags$p("MetaboAnalyst maps common lipid names. Supply an annotation table that connects the current MSI m/z features to lipid names; the app will not convert an unannotated m/z into a unique lipid identity."),
      fileInput("lipid_annotation_csv", "Lipid annotation CSV",
        accept = ".csv", placeholder = "Needs lipid_name plus mz (or an mz_* feature column)"),
      downloadButton("download_lipid_annotation_template", "Download current-feature annotation template (.csv)"),
      fluidRow(
        column(4, numericInput("lipid_match_ppm", "Annotation match tolerance (ppm)", 5, min = .Machine$double.eps, step = 1)),
        column(4, numericInput("lipid_p_cutoff", "Differential p-value cutoff", .05, min = 0, max = 1, step = .01)),
        column(4, br(), uiOutput("lipid_export_action"))),
      uiOutput("lipid_export_status"), uiOutput("lipid_export_downloads"))
  })
  output$lipid_export_action <- renderUI({
    if (is.null(state$functional_peak_list))
      return(tags$div(class = "alert alert-warning", "A differential lipid-name export requires exactly two ROI groups."))
    if (is.null(input$lipid_annotation_csv))
      return(tags$div(class = "alert alert-secondary", "Upload a lipid annotation CSV first."))
    actionButton("prepare_lipid_export", "Prepare lipid export", class = "btn-primary")
  })
  observeEvent(input$prepare_lipid_export, {
    req(state$functional_peak_list, input$lipid_annotation_csv)
    annotation <- read.csv(input$lipid_annotation_csv$datapath, check.names = FALSE, stringsAsFactors = FALSE)
    result <- prepare_lipid_name_export(annotation, state$functional_peak_list,
      ppm = input$lipid_match_ppm, p_cutoff = input$lipid_p_cutoff)
    state$lipid_export <- result$audit
    state$lipid_name_list <- result$names
    write.csv(result$audit, file.path(session_dir, "metaboanalyst_lipid_mapping_audit.csv"), row.names = FALSE)
    writeLines(result$names, file.path(session_dir, "metaboanalyst_lipid_name_list.txt"))
  })
  output$lipid_export_status <- renderUI({
    req(state$lipid_export)
    if (!length(state$lipid_name_list)) return(tags$div(class = "alert alert-warning",
      "Annotations matched, but no lipid name passes the selected p-value cutoff."))
    tags$div(class = "alert alert-success", sprintf(
      "Mapped %d annotation rows; %d unique lipid names are ready for MetaboAnalyst Enrichment Analysis → Lipidomics.",
      nrow(state$lipid_export), length(state$lipid_name_list)))
  })
  output$lipid_export_downloads <- renderUI(if (!is.null(state$lipid_export)) tagList(
    tags$div(class = "alert alert-info",
      tags$strong("Upload destination matters: "),
      "the one-column lipid-name TXT is for MetaboAnalyst Enrichment Analysis → metabolite/lipid list, then select the Lipidomics library. Do not upload it to Statistical Analysis / Data Processing, which requires a sample-by-feature matrix and rejects one-column input."),
    downloadButton("download_lipid_names", "MetaboAnalyst lipid-name list (.txt)"),
    downloadButton("download_lipid_audit", "Lipid mapping audit (.csv)")) else NULL)
  output$download_lipid_names <- downloadHandler(filename = function() "metaboanalyst_lipid_name_list.txt",
    content = function(file) { req(length(state$lipid_name_list)); writeLines(state$lipid_name_list, file) })
  output$download_lipid_audit <- downloadHandler(filename = function() "metaboanalyst_lipid_mapping_audit.csv",
    content = function(file) { req(state$lipid_export); write.csv(state$lipid_export, file, row.names = FALSE) })
  output$download_lipid_annotation_template <- downloadHandler(
    filename = function() "msi_lipid_annotation_template.csv",
    content = function(file) {
      req(state$processed)
      metadata <- state$processed$feature_metadata
      write.csv(data.frame(
        column_name = metadata$column_name, mz = metadata$mz,
        ion_mode = state$spec$ion_mode %or% NA_character_,
        lipid_name = "", lipid_id = "", annotation_source = "",
        stringsAsFactors = FALSE
      ), file, row.names = FALSE)
    })

  observeEvent(input$metabo_result_csv, {
    raw <- read.csv(input$metabo_result_csv$datapath, check.names = FALSE, stringsAsFactors = FALSE)
    state$metabo_result <- normalize_metaboanalyst_result(raw, source_name = input$metabo_result_csv$name %or% "")
    state$metabo_result_type <- detect_metaboanalyst_result(state$metabo_result)
    state$mapped_scores <- NULL
    write.csv(state$metabo_result, file.path(session_dir, "metaboanalyst_result_normalized.csv"), row.names = FALSE)
  })
  output$metabo_result_controls <- renderUI({
    req(state$metabo_result)
    scores <- names(state$metabo_result)[grepl("^(PC|Comp)[0-9]+$", names(state$metabo_result))]
    if (length(scores)) selectInput("metabo_score_column", "Score to map", scores) else
      numericInput("metabo_top_n", "Top features", 20, min = 1, step = 1)
  })
  observe({
    req(state$metabo_result, state$sampling)
    scores <- names(state$metabo_result)[grepl("^(PC|Comp)[0-9]+$", names(state$metabo_result))]
    if (!length(scores)) return()
    score_column <- input$metabo_score_column %or% scores[1]
    req(score_column %in% scores)
    mapped <- backmap_sample_scores(state$metabo_result, state$sampling$sample_mapping, score_column)
    pixels <- analysis_pixel_matrix()
    pixels$pixel_id <- as.character(pixels$pixel_id)
    mapped$pixel_id <- as.character(mapped$pixel_id)
    state$mapped_scores <- merge(pixels[, c("pixel_id", "x", "y"), drop = FALSE], mapped,
      by = "pixel_id", all = FALSE, sort = FALSE)
    write.csv(state$mapped_scores, file.path(session_dir, paste0("spatial_backmap_", score_column, ".csv")), row.names = FALSE)
  })
  output$metabo_result_gate <- renderUI({
    if (is.null(state$metabo_result)) return(tags$div(class = "alert alert-secondary", "No statistical result is loaded."))
    tags$div(class = "alert alert-success", paste("Detected result type:", state$metabo_result_type))
  })
  output$metabo_results_ui <- renderUI(if (!is.null(state$metabo_result)) tagList(
    uiOutput("metabo_result_download"),
    fluidRow(
      column(6, tags$div(class = "result-card", plotOutput("metabo_result_plot", height = 430))),
      column(6, tags$div(class = "result-card", plotOutput("metabo_spatial_plot", height = 430)))),
    tags$details(class = "result-card", tags$summary("Preview imported statistical result"),
      DTOutput("metabo_result_table"))
  ) else NULL)
  output$metabo_result_table <- renderDT({ req(state$metabo_result)
    datatable(head(state$metabo_result, 100), options = list(pageLength = 10, scrollX = TRUE)) })
  output$metabo_result_plot <- renderPlot({
    req(state$metabo_result)
    if (state$metabo_result_type %in% c("pca_scores", "plsda_scores")) {
      print(plot_metabo_score_plot(state$metabo_result,
        sample_matrix = if (is.null(state$sampling)) NULL else state$sampling$sample_matrix))
    } else if (state$metabo_result_type %in% c("vip", "differential")) {
      print(plot_metabo_feature_rank(state$metabo_result, top_n = input$metabo_top_n %or% 20))
    } else validate(need(FALSE, "Unsupported result structure."))
  })
  output$metabo_spatial_plot <- renderPlot({
    req(state$metabo_result)
    if (state$metabo_result_type %in% c("pca_scores", "plsda_scores")) {
      req(state$sampling)
      print(plot_sample_score_map(state$metabo_result, state$sampling$sample_mapping,
        analysis_pixel_matrix(), score_column = input$metabo_score_column))
    } else if (state$metabo_result_type %in% c("vip", "differential")) {
      print(plot_metabo_feature_view(state$metabo_result, analysis_pixel_matrix(),
        top_n = input$metabo_top_n %or% 20, transform = "identity"))
    }
  })
  output$metabo_result_download <- renderUI(if (!is.null(state$metabo_result)) tagList(
    downloadButton("download_normalized_metabo_result", "Download normalized result"),
    if (!is.null(state$mapped_scores)) downloadButton("download_mapped_scores", "Download spatial back-map")
  ) else NULL)
  output$download_normalized_metabo_result <- downloadHandler(filename = function() "metaboanalyst_result_normalized.csv",
    content = function(file) { req(state$metabo_result); write.csv(state$metabo_result, file, row.names = FALSE) })
  output$download_mapped_scores <- downloadHandler(filename = function() "spatial_backmapped_scores.csv",
    content = function(file) { req(state$mapped_scores); write.csv(state$mapped_scores, file, row.names = FALSE) })

  output$corroboration_action <- renderUI({
    if (is.null(state$domains)) return(tags$div(class = "alert alert-warning", "Load or generate spatial domains first."))
    if (is.null(input$corroboration_csv)) return(tags$div(class = "alert alert-secondary", "Supply an independent coordinate-aligned label CSV."))
    actionButton("run_corroboration", "Compare label maps", class = "btn-primary")
  })
  observeEvent(input$run_corroboration, {
    req(state$domains, input$corroboration_csv, state$processed)
    info <- validate_label_file(input$corroboration_csv$datapath)
    second <- resolve_label_mapping(info, state$processed$coordinates,
      state$spec$sample_id, state$spec$section_id)
    label_a <- as.character(state$domains$domain_id)
    label_b <- as.character(second$domain_id)
    label_a[label_a == "-1"] <- NA_character_
    label_b[label_b == "-1"] <- NA_character_
    result <- corroborate_cluster_labels(label_a, label_b,
      min_cooccurrence_fraction = input$corroboration_fraction,
      min_pair_count = input$corroboration_count,
      mapping = input$corroboration_mapping)
    state$corroboration <- result
    state$corroborated_pixels <- data.frame(
      state$processed$coordinates, label_a = label_a, label_b = label_b,
      combined_label = result$combined_label,
      corroborated = result$corroborated, stringsAsFactors = FALSE)
    write.csv(result$pair_summary, file.path(session_dir, "cross_modal_label_pairs.csv"), row.names = FALSE)
    write.csv(state$corroborated_pixels, file.path(session_dir, "cross_modal_corroborated_pixels.csv"), row.names = FALSE)
  })
  output$corroboration_gate <- renderUI({
    req(state$corroboration)
    tags$div(class = "alert alert-success",
      sprintf("Corroborated %.1f%% of pixels carrying both labels. %s",
        100 * state$corroboration$corroboration_rate_valid_pixels,
        state$corroboration$interpretation))
  })
  output$corroboration_results_ui <- renderUI(if (!is.null(state$corroboration)) tagList(
    uiOutput("corroboration_download"),
    fluidRow(
      column(6, tags$div(class = "result-card", plotOutput("corroboration_plot", height = 430))),
      column(6, tags$div(class = "result-card", DTOutput("corroboration_pairs"))))
  ) else NULL)
  output$corroboration_pairs <- renderDT({ req(state$corroboration)
    datatable(state$corroboration$pair_summary, options = list(pageLength = 10, scrollX = TRUE)) })
  output$corroboration_plot <- renderPlot({ req(state$corroborated_pixels); d <- state$corroborated_pixels
    ggplot(d, aes(x, y, color = combined_label)) + geom_point(size = .7) + coord_equal() +
      scale_y_reverse() + theme_minimal() + labs(title = "Cross-modal label corroboration", color = "Joint label") })
  output$corroboration_download <- renderUI(if (!is.null(state$corroboration))
    downloadButton("download_corroboration", "Download corroborated pixels") else NULL)
  output$download_corroboration <- downloadHandler(filename = function() "cross_modal_corroborated_pixels.csv",
    content = function(file) { req(state$corroborated_pixels); write.csv(state$corroborated_pixels, file, row.names = FALSE) })

  output$moran_action <- renderUI({
    if (!isTRUE(state$valid) || !isTRUE(state$validation$capabilities$moran)) return(tags$div(class = "alert alert-secondary", "Moran's I is disabled until MSI input is validated."))
    if (is.null(state$processed)) return(tags$div(class = "alert alert-warning", "Process MSI first."))
    actionButton("run_moran", if (is.null(state$moran)) "Compute Moran's I" else "Recompute Moran's I", class = "btn-primary")
  })
  moran_work_units <- reactive({
    if (is.null(state$analysis_matrix)) return(NA_real_)
    keep <- if (is.null(state$tissue_mask)) rep(TRUE, nrow(state$analysis_matrix)) else state$tissue_mask
    sum(keep) * ncol(state$analysis_matrix) * (input$permutations %or% 0)
  })
  output$moran_workload_warning <- renderUI({
    units <- moran_work_units()
    if (!is.finite(units)) return(NULL)
    message <- sprintf(
      "Requested workload: %s tissue pixels × %s features × %s permutations = %.2f billion pixel-feature permutations.",
      format(sum(if (is.null(state$tissue_mask)) rep(TRUE, nrow(state$analysis_matrix)) else state$tissue_mask), big.mark = ","),
      format(ncol(state$analysis_matrix), big.mark = ","), format(input$permutations, big.mark = ","), units / 1e9)
    if (units <= 2e8) return(tags$div(class = "alert alert-info", message))
    tags$div(class = "alert alert-danger",
      tags$strong("This Moran permutation job may take a very long time in an interactive Shiny session."),
      tags$p(message),
      tags$p("The statistical method has not been changed. Reduce permutations only if scientifically acceptable, or run the full job offline."),
      checkboxInput("confirm_large_moran", "I understand the cost and want to run this full permutation job", FALSE))
  })
  observeEvent(input$run_moran, {
    req(state$processed, state$analysis_matrix)
    keep <- if (is.null(state$tissue_mask)) rep(TRUE, nrow(state$analysis_matrix)) else state$tissue_mask
    validate(need(sum(keep) >= 3L, "At least three tissue pixels are required for Moran analysis."))
    coordinates <- state$processed$coordinates[keep, , drop = FALSE]; matrix <- state$analysis_matrix[keep, , drop = FALSE]
    work_units <- nrow(matrix) * ncol(matrix) * input$permutations
    validate(need(work_units <= 2e8 || isTRUE(input$confirm_large_moran),
      "This is a large permutation job. Review the workload warning and explicitly confirm it before running."))
    withProgress(message = "Computing Moran's I", value = 0, {
      setProgress(value = .10, detail = "Building the reusable spatial-neighbor graph")
      threshold <- if (input$neighbor_method == "queen") sqrt(2) else 1
      state$neighbors <- build_spatial_neighbors(coordinates$x, coordinates$y, input$neighbor_method, threshold)
      setProgress(value = .35, detail = sprintf("Scoring %d m/z features", ncol(matrix)))
      pixel <- data.frame(coordinates[c("pixel_id", "sample_id", "section_id", "x", "y")], matrix, check.names = FALSE)
      names(pixel)[-(1:5)] <- state$processed$feature_metadata$column_name
      state$moran <- compute_spatially_variable_metabolites(pixel, coordinates, n_perm = input$permutations,
        alternative = "two.sided", p_adjust_method = "BH", seed = 20260807,
        neighbors = state$neighbors)
      setProgress(value = .95, detail = "Saving the ranked spatial-feature table")
      write.csv(state$moran, file.path(session_dir, paste0("moran_", input$neighbor_method, ".csv")), row.names = FALSE)
      setProgress(value = 1, detail = "Moran analysis complete")
    })
  })
  observeEvent(input$load_moran, {
    validate(need(!is.null(state$processed), "Process the current MSI before attaching a previous result."))
    directory <- normalizePath(input$moran_result_dir, mustWork = FALSE)
    path <- file.path(directory, paste0("spatially_variable_metabolites_", input$neighbor_method, ".csv"))
    validate(need(file.exists(path), "The selected adjacency result CSV was not found.")); state$moran <- read.csv(path, check.names = FALSE)
  })
  output$moran_gate <- renderUI(if (is.null(state$moran)) tags$div(class = "alert alert-warning", "No Moran result is loaded for the current MSI.") else NULL)
  moran_concordance <- reactive({
    req(state$moran, state$metabo_result)
    validate(need(state$metabo_result_type %in% c("vip", "differential"),
      "Import a feature-level VIP or differential result to compare with Moran's I."))
    validate(need("Feature" %in% names(state$metabo_result), "The statistical result has no Feature column."))
    feature_key <- function(x) {
      raw <- sub("^mz_", "", as.character(x), ignore.case = TRUE)
      number <- suppressWarnings(as.numeric(sub("__.*$", "", raw)))
      ifelse(is.finite(number), sprintf("%.6f", number), tolower(raw))
    }
    statistical <- state$metabo_result
    statistical$.feature_key <- feature_key(statistical$Feature)
    spatial <- state$moran
    spatial$.feature_key <- feature_key(spatial$feature)
    merged <- merge(statistical, spatial, by = ".feature_key", all.x = TRUE, sort = FALSE)
    merged$spatial_fdr_significant <- is.finite(merged$adj_p_value) & merged$adj_p_value <= 0.05
    merged$spatial_coherence_class <- ifelse(is.na(merged$morans_i), "not matched to MSI feature",
      ifelse(merged$spatial_fdr_significant & merged$morans_i > 0, "positive spatial coherence", "weak/non-positive spatial coherence"))
    merged
  })
  output$moran_validation_gate <- renderUI({
    if (is.null(state$moran)) return(tags$div(class = "alert alert-warning", "Compute Moran's I in Processing & QC first."))
    if (is.null(state$metabo_result)) return(tags$div(class = "alert alert-info", "Moran's I is cached. Import a feature-level statistical result in Step 6 to check concordance."))
    tags$div(class = "alert alert-success", "Cached Moran results and the imported feature-level result are ready for descriptive concordance checking.")
  })
  output$moran_concordance_ui <- renderUI({
    result <- moran_concordance()
    tagList(
      tags$div(class = "alert alert-info", sprintf(
        "%d of %d imported features show positive Moran spatial coherence at FDR ≤ 0.05. This is a coherence check, not independent validation.",
        sum(result$spatial_coherence_class == "positive spatial coherence"), nrow(result))),
      tags$div(class = "result-card", DTOutput("moran_concordance_table")),
      downloadButton("download_moran_concordance", "Download statistical × Moran concordance"))
  })
  output$moran_concordance_table <- renderDT({
    result <- moran_concordance()
    datatable(result, options = list(pageLength = 15, scrollX = TRUE))
  })
  output$download_moran_concordance <- downloadHandler(
    filename = function() "statistical_moran_spatial_coherence.csv",
    content = function(file) write.csv(moran_concordance(), file, row.names = FALSE))
  output$moran_screening_ui <- renderUI(if (!is.null(state$moran)) tagList(
    tags$div(class = "alert alert-success", sprintf("Moran's I is cached for %d m/z features; %d pass FDR 0.05.",
      nrow(state$moran), sum(is.finite(state$moran$adj_p_value) & state$moran$adj_p_value <= 0.05))),
    selectInput("moran_feature_view", "Ion image to inspect",
      choices = setNames(state$moran$feature, sprintf("%s · I=%.3f · FDR=%.3g", state$moran$feature, state$moran$morans_i, state$moran$adj_p_value))),
    fluidRow(
      column(7, tags$div(class = "result-card", plotOutput("ion_plot", height = 480))),
      column(5, tags$div(class = "result-card", h4("Top spatial features"), DTOutput("moran_table")))),
    tags$details(class = "result-card", tags$summary("Neighbor-graph diagnostics"), DTOutput("neighbor_table"))
  ) else NULL)
  output$neighbor_table <- renderDT({ req(state$neighbors); datatable(spatial_neighbor_diagnostics(state$neighbors), options = list(pageLength = 10)) })
  output$moran_table <- renderDT({ req(state$moran); datatable(head(state$moran[order(state$moran$adj_p_value, -abs(state$moran$morans_i)), ], 100), options = list(pageLength = 10)) })
  output$ion_plot <- renderPlot({ req(state$moran, state$processed); column <- input$moran_feature_view %or% state$moran$feature[1]
    d <- data.frame(state$processed$coordinates, intensity = state$analysis_matrix[, match(column, state$processed$feature_metadata$column_name)])
    ggplot(d, aes(x, y, color = intensity)) + geom_point(size = .7) + coord_equal() + scale_y_reverse() + scale_color_viridis_c() + theme_minimal() + labs(title = paste("Top Moran feature", column)) })

  output$lcms_action <- renderUI(if (isTRUE(state$valid) && !is.null(state$processed) && present_path(active_lcms_path()))
    actionButton("load_lcms", "Read precursor window", class = "btn-primary") else tags$div(class = "alert alert-secondary", "LC-MS/MS evidence is disabled. Supply and validate an mzML file."))
  observeEvent(input$load_lcms, {
    req(state$valid, state$processed); withProgress(message = "Reading target MS2 window", value = .1, {
      state$lcms <- read_mzml_fragment_spectra(active_lcms_path(), input$lcms_precursor, 10); incProgress(.9)
      if (nrow(state$lcms$precursor_scan_metadata)) {
        state$lcms$precursor_scan_metadata$msi_to_reference_ppm <- (input$msi_mz - input$lcms_precursor) / input$lcms_precursor * 1e6
        state$lcms$precursor_scan_metadata$evidence_relationship <- active_lcms_relationship()
        write.csv(state$lcms$precursor_scan_metadata, file.path(session_dir, "precursor_scan_metadata.csv"), row.names = FALSE)
        write.csv(state$lcms$fragment_peak_table, file.path(session_dir, "fragment_peak_table.csv"), row.names = FALSE)
      }
    })
  })
  output$lcms_gate <- renderUI(if (is.null(state$lcms)) tags$div(class = "alert alert-warning", "No LC-MS/MS target window has been read.") else NULL)
  output$precursor_table <- renderDT({ req(state$lcms); datatable(state$lcms$precursor_scan_metadata, options = list(scrollX = TRUE)) })
  output$fragment_plot <- renderPlot({ req(state$lcms); d <- state$lcms$fragment_peak_table
    validate(need(nrow(d) > 0, "No target-precursor MS2 scan was found in this file."))
    ggplot(d, aes(fragment_mz, fragment_intensity)) + geom_segment(aes(xend = fragment_mz, yend = 0)) + facet_wrap(~scan_id, scales = "free_y") + theme_minimal() + labs(title = "Unassigned fragment spectrum", x = "fragment m/z", y = "intensity") })

  observeEvent(input$run_feature_matching, {
    req(state$processed, input$msi_feature_csv, input$lcms_feature_csv)
    validate(need(!identical(active_lcms_relationship(), "none"),
      "Declare the LC-MS/MS evidence relationship in Data setup before matching features."))
    msi <- read.csv(input$msi_feature_csv$datapath, check.names = FALSE, stringsAsFactors = FALSE)
    lcms <- read.csv(input$lcms_feature_csv$datapath, check.names = FALSE, stringsAsFactors = FALSE)
    require_columns(msi, "mz", "MSI feature CSV")
    require_columns(lcms, "mz", "LC-MS feature CSV")
    state$feature_matches <- cross_validate_msi_lcms(msi, lcms, ppm = input$feature_match_ppm,
      assignment_method = input$assignment_method)
    state$feature_matches$evidence_relationship <- active_lcms_relationship()
    write.csv(state$feature_matches, file.path(session_dir, "msi_lcms_feature_matches.csv"), row.names = FALSE)
    summary <- attr(state$feature_matches, "summary")
    if (!is.null(summary)) write.csv(summary, file.path(session_dir, "msi_lcms_feature_match_summary.csv"), row.names = FALSE)
  })
  output$feature_match_gate <- renderUI({ req(state$feature_matches); summary <- attr(state$feature_matches, "summary")
    tags$div(class = "alert alert-success", sprintf("Matched %d of %d MSI features using %s one-to-one assignment.",
      summary$matched_features, summary$total_msi_features, summary$assignment_method)) })
  output$feature_match_table <- renderDT({ req(state$feature_matches)
    datatable(state$feature_matches, options = list(pageLength = 10, scrollX = TRUE)) })
  output$feature_match_download <- renderUI(if (!is.null(state$feature_matches))
    downloadButton("download_feature_matches", "Download feature matches") else NULL)
  output$download_feature_matches <- downloadHandler(filename = function() "msi_lcms_feature_matches.csv",
    content = function(file) { req(state$feature_matches); write.csv(state$feature_matches, file, row.names = FALSE) })

  observeEvent(input$run_ccs_validation, {
    req(state$processed, input$observed_ccs_csv, input$reference_ccs_csv)
    validate(need(!identical(active_lcms_relationship(), "none"),
      "Declare the LC-MS/MS evidence relationship in Data setup before evaluating CCS."))
    observed <- read.csv(input$observed_ccs_csv$datapath, check.names = FALSE, stringsAsFactors = FALSE)
    reference <- read.csv(input$reference_ccs_csv$datapath, check.names = FALSE, stringsAsFactors = FALSE)
    state$ccs_evidence <- validate_ccs_evidence(observed, reference, ccs_tolerance_pct = input$ccs_tolerance)
    state$ccs_evidence$evidence_relationship <- active_lcms_relationship()
    write.csv(state$ccs_evidence, file.path(session_dir, "ccs_evidence.csv"), row.names = FALSE)
  })
  output$ccs_gate <- renderUI(if (!is.null(state$ccs_evidence)) tags$div(class = "alert alert-success",
    "CCS evidence was evaluated. Measured-library and predicted references remain explicitly distinguished in the output.") else NULL)
  output$ccs_table <- renderDT({ req(state$ccs_evidence)
    datatable(state$ccs_evidence, options = list(pageLength = 10, scrollX = TRUE)) })
  output$ccs_download <- renderUI(if (!is.null(state$ccs_evidence)) downloadButton("download_ccs", "Download CCS evidence") else NULL)
  output$download_ccs <- downloadHandler(filename = function() "ccs_evidence.csv",
    content = function(file) { req(state$ccs_evidence); write.csv(state$ccs_evidence, file, row.names = FALSE) })

  output$download_gate <- renderUI(if (!isTRUE(state$valid)) tags$div(class = "alert alert-warning", "Validate an analysis before downloading results.") else NULL)
  output$download_bundle <- downloadHandler(filename = function() paste0("SpatialOmicsMSI-session-", Sys.Date(), ".zip"), content = function(file) {
    files <- list.files(session_dir, full.names = TRUE, recursive = TRUE); validate(need(length(files) > 0, "Run at least one workflow step first."))
    old <- setwd(session_dir); on.exit(setwd(old)); utils::zip(file, sub(paste0("^", session_dir, "/"), "", files))
  })
}

shinyApp(ui, server)
