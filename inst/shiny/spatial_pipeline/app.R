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
    make_tissue_mask = FALSE,
    source_note = "MassIVE MSV000090179; CC0; female mouse brain, 12 weeks, replicate 1."
  )
)

preset_spec <- function(key) {
  if (!key %in% names(example_datasets)) stop("Unknown example dataset.", call. = FALSE)
  example_datasets[[key]]
}

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
  data <- utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  require_columns(data, c("x", "y"), "Processed pixel-by-feature CSV")
  if (!nrow(data) || any(!is.finite(data$x) | !is.finite(data$y)) || anyDuplicated(data[c("x", "y")])) {
    stop("CSV x/y coordinates must be finite, non-empty and unique.", call. = FALSE)
  }
  metadata <- c("pixel_id", "sample_id", "section_id", "subject_id", "x", "y", "domain_id", "domain_label")
  candidates <- setdiff(names(data), metadata)
  numeric_feature <- vapply(data[candidates], is.numeric, logical(1))
  features <- candidates[numeric_feature]
  if (!length(features)) stop("CSV must contain at least one numeric feature column in addition to x/y.", call. = FALSE)
  normalized <- normalize_feature_names(features)
  normalized <- make.unique(normalized, sep = "__")
  values <- as.data.frame(data[features], check.names = FALSE)
  names(values) <- normalized
  pixel_id <- if ("pixel_id" %in% names(data)) data$pixel_id else seq_len(nrow(data))
  if (anyNA(pixel_id) || anyDuplicated(pixel_id)) stop("CSV pixel_id must be non-missing and unique.", call. = FALSE)
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
  if (!has_msi && !has_lcms) errors <- c(errors, "Provide MSI data, an LC-MS/MS mzML file, or both.")
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
      notes <- c(notes, "Processed pixel-by-feature CSV structure passed.")
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

ui <- navbarPage(
  id = "workflow_tabs", title = "SpatialOmicsMSI",
  header = tags$div(class = "alert alert-info", style = "margin:8px",
    strong("General-purpose spatial metabolomics workflow."),
    " Example datasets demonstrate complementary modules and are not a matched cohort."),
  tabPanel("1 Start",
    radioButtons("entry", NULL, inline = TRUE,
      choices = c("Start a new analysis" = "new", "Load an example dataset" = "example"), selected = "new"),
    conditionalPanel("input.entry == 'example'",
      selectInput("example_key", "Example datasets", choices = setNames(names(example_datasets), vapply(example_datasets, `[[`, "", "label"))),
      actionButton("load_example", "Use this example", class = "btn-info")),
    hr(),
    radioButtons("data_origin", "Input location", inline = TRUE,
      choices = c("Server files" = "server", "Browser upload" = "upload"), selected = "server"),
    selectInput("msi_input_type", "Primary MSI input", choices = c("Paired imzML + ibd" = "imzml",
      "Processed pixel × feature CSV" = "csv", "No MSI (LC-MS/MS only)" = "none")),
    conditionalPanel("input.data_origin == 'upload' && input.msi_input_type == 'imzml'",
      fileInput("upload_imzml", "imzML", accept = c(".imzML", ".imzml")), fileInput("upload_ibd", "ibd", accept = c(".ibd", ".IBD"))),
    conditionalPanel("input.data_origin == 'upload' && input.msi_input_type == 'csv'", fileInput("upload_csv", "Pixel × feature CSV", accept = ".csv")),
    conditionalPanel("input.data_origin == 'upload'",
      fileInput("upload_optical", "Optional optical/H&E/brightfield JPEG", accept = c(".jpg", ".jpeg")),
      fileInput("upload_transform", "Optional registration transform JSON", accept = ".json"),
      fileInput("upload_labels", "Optional ROI/domain CSV", accept = ".csv"),
      fileInput("upload_lcms", "Optional LC-MS/MS mzML", accept = c(".mzML", ".mzml"))),
    tags$details(tags$summary("Technical details — server paths"),
      conditionalPanel("input.data_origin == 'server' && input.msi_input_type == 'imzml'",
        textInput("server_imzml", "imzML path"), textInput("server_ibd", "ibd path")),
      conditionalPanel("input.data_origin == 'server' && input.msi_input_type == 'csv'", textInput("server_csv", "Pixel × feature CSV path")),
      conditionalPanel("input.data_origin == 'server'",
        textInput("server_optical", "Optional optical JPEG path"), textInput("server_transform", "Optional transform JSON path"),
        textInput("server_labels", "Optional ROI/domain CSV path"), textInput("server_lcms", "Optional LC-MS/MS mzML path"),
        textInput("server_attribution", "Optional attribution metadata path"))),
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
  tabPanel("2 Input & provenance", uiOutput("module_overview"), DTOutput("provenance_table"), textOutput("source_note")),
  tabPanel("3 Processing & QC", checkboxInput("tic_normalize", "TIC normalization", TRUE),
    checkboxInput("log1p_transform", "log10(x + 1) transformation", TRUE), uiOutput("processing_action"),
    uiOutput("processing_gate"), DTOutput("qc_table"), DTOutput("feature_preview"),
    plotOutput("tissue_mask_plot", height = 600), DTOutput("tissue_mask_diagnostics"), uiOutput("tissue_download")),
  tabPanel("4 Brightfield–MSI registration", uiOutput("registration_action"), uiOutput("registration_gate"),
    imageOutput("registration_plot", height = "650px"), DTOutput("registration_diagnostics")),
  tabPanel("5 Spatial domains", fileInput("domain_csv_runtime", "Add/replace ROI/domain CSV", accept = ".csv"),
    numericInput("domain_k", "Exploratory domain count k", 4, min = 2, max = 12),
    numericInput("domain_seed", "Random seed", 20260808, min = 1), numericInput("domain_pcs", "PCA components", 10, min = 2, max = 30),
    uiOutput("domain_action"), tags$div(class = "alert alert-warning",
      "Pseudoreplication warning: pixels are not biological replicates. Generated domains are descriptive and data-driven, never anatomical ROI."),
    uiOutput("domain_gate"), plotOutput("domain_plot", height = 600), DTOutput("domain_counts"), DTOutput("domain_features"), uiOutput("domain_download")),
  tabPanel("6 Cross-modal ROI validation",
    tags$div(class = "alert alert-info",
      "Compare the active spatial-domain labels with an independently generated label map. This validates label agreement; it does not turn exploratory clusters into anatomical ground truth."),
    fileInput("corroboration_csv", "Independent ROI/domain CSV (x, y and a label column)", accept = ".csv"),
    selectInput("corroboration_mapping", "Mapping rule",
      c("Mutual best match" = "mutual_best", "One-way best match" = "one_way")),
    fluidRow(
      column(6, numericInput("corroboration_fraction", "Minimum conditional overlap", 0.5, min = 0, max = 1, step = 0.05)),
      column(6, numericInput("corroboration_count", "Minimum shared pixels", 10, min = 1, step = 1))
    ),
    uiOutput("corroboration_action"), uiOutput("corroboration_gate"),
    DTOutput("corroboration_pairs"), plotOutput("corroboration_plot", height = 600),
    uiOutput("corroboration_download")),
  tabPanel("7 Moran's I", selectInput("neighbor_method", "Adjacency", c("Queen" = "queen", "Rook" = "rook")),
    numericInput("permutations", "Two-sided permutations", 499, min = 99, step = 100), uiOutput("moran_action"),
    tags$details(tags$summary("Technical details — load a previously generated result"),
      textInput("moran_result_dir", "Result directory", Sys.getenv("SPATIALOMICS_OMIX_MORAN_DIR", "")), actionButton("load_moran", "Load result")),
    uiOutput("moran_gate"), DTOutput("neighbor_table"), DTOutput("moran_table"), plotOutput("ion_plot", height = 600)),
  tabPanel("8 Multimodal evidence", numericInput("msi_mz", "MSI observed m/z", 775.55261535, step = .00000001),
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
    uiOutput("ccs_gate"), DTOutput("ccs_table"), uiOutput("ccs_download")),
  tabPanel("9 Download results", tags$p("Each session writes to a unique temporary directory; source data are read-only."),
    uiOutput("download_gate"), downloadButton("download_bundle", "Download session bundle"),
    tags$details(tags$summary("Technical details — temporary session location"), textOutput("session_path")))
)

server <- function(input, output, session) {
  state <- reactiveValues(valid = FALSE, spec = NULL, validation = NULL, processed = NULL,
    analysis_matrix = NULL, tissue_gate = NULL, tissue_mask = NULL, registration = NULL,
    domains = NULL, domain_counts = NULL, domain_features = NULL, neighbors = NULL,
    corroboration = NULL, corroborated_pixels = NULL, moran = NULL, lcms = NULL,
    feature_matches = NULL, ccs_evidence = NULL, example_note = NULL)
  session_dir <- file.path(session_root, paste0("SpatialOmicsMSI-session-", Sys.getpid(), "-", substr(session$token, 1, 8)))
  dir.create(session_dir, recursive = TRUE, showWarnings = FALSE)
  output$session_path <- renderText(session_dir)

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
      ion_mode_source = input$ion_mode_source %or% "", spectral_processing = input$spectral_processing,
      alignment_ppm = input$alignment_ppm, peak_pick_snr = input$peak_pick_snr,
      min_detection_fraction = input$min_detection_fraction, make_tissue_mask = isTRUE(input$make_tissue_mask))
  })

  reset_analysis <- function() {
    state$valid <- FALSE; state$validation <- NULL; state$processed <- NULL; state$analysis_matrix <- NULL
    state$tissue_gate <- NULL; state$tissue_mask <- NULL; state$registration <- NULL
    state$domains <- NULL; state$domain_counts <- NULL; state$domain_features <- NULL
    state$corroboration <- NULL; state$corroborated_pixels <- NULL
    state$neighbors <- NULL; state$moran <- NULL; state$lcms <- NULL
    state$feature_matches <- NULL; state$ccs_evidence <- NULL
  }

  observeEvent(input$load_example, {
    p <- preset_spec(input$example_key); reset_analysis(); state$example_note <- p$source_note
    updateRadioButtons(session, "data_origin", selected = "server")
    updateSelectInput(session, "msi_input_type", selected = p$msi_input_type)
    fields <- c(imzml = "server_imzml", ibd = "server_ibd", csv = "server_csv", optical = "server_optical",
      transform = "server_transform", labels = "server_labels", lcms = "server_lcms", attribution = "server_attribution")
    for (name in names(fields)) updateTextInput(session, fields[[name]], value = p[[paste0(name, "_path")]] %or% "")
    for (name in c("sample_id", "section_id", "subject_id", "ion_mode_source")) updateTextInput(session, name, value = p[[name]])
    updateSelectInput(session, "ion_mode", selected = p$ion_mode)
    updateSelectInput(session, "spectral_processing", selected = p$spectral_processing)
    updateNumericInput(session, "alignment_ppm", value = p$alignment_ppm)
    updateNumericInput(session, "peak_pick_snr", value = p$peak_pick_snr)
    updateNumericInput(session, "min_detection_fraction", value = p$min_detection_fraction)
    updateCheckboxInput(session, "make_tissue_mask", value = p$make_tissue_mask)
  })

  observeEvent(input$validate_input, {
    reset_analysis(); state$spec <- current_spec()
    state$validation <- withProgress(message = "Checking input metadata", value = .2, {
      out <- validate_analysis_spec(state$spec, inspect = TRUE); incProgress(.8); out
    })
    state$valid <- state$validation$valid
  })

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
  output$continue_processing <- renderUI(if (isTRUE(state$valid) && isTRUE(state$validation$capabilities$msi))
    actionButton("continue_button", "Continue to Processing", class = "btn-success") else NULL)
  observeEvent(input$continue_button, updateNavbarPage(session, "workflow_tabs", selected = "3 Processing & QC"))
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
    datatable(data.frame(field = c("sample_id", "section_id", "subject_id", "ion_mode", "ion_mode_source", "input_type"),
      value = c(s$sample_id, s$section_id, s$subject_id, s$ion_mode, s$ion_mode_source, s$msi_input_type)), options = list(dom = "t")) })

  output$processing_action <- renderUI(if (isTRUE(state$valid) && isTRUE(state$validation$capabilities$processing))
    actionButton("run_processing", "Run processing", class = "btn-primary") else tags$div(class = "alert alert-secondary", "Processing is available after validating MSI input."))
  observeEvent(input$run_processing, {
    req(state$valid); s <- state$spec
    withProgress(message = "Preparing MSI", value = 0, {
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
            progress_callback = function(value, detail) incProgress(value, detail = detail))
        }
      }
      feature_names <- state$processed$feature_metadata$column_name
      matrix <- as.matrix(state$processed$pixel_feature_matrix[, feature_names, drop = FALSE]); storage.mode(matrix) <- "double"
      tic <- rowSums(matrix, na.rm = TRUE); analysis <- matrix
      if (isTRUE(input$tic_normalize)) analysis <- analysis / pmax(tic, .Machine$double.eps)
      if (isTRUE(input$log1p_transform)) analysis <- log10(analysis + 1)
      state$analysis_matrix <- analysis; state$tic <- tic
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
      write.csv(state$processed$qc_summary, file.path(session_dir, "qc_summary.csv"), row.names = FALSE)
      write.csv(state$processed$feature_metadata, file.path(session_dir, "feature_metadata.csv"), row.names = FALSE)
      write.csv(state$processed$coordinates, file.path(session_dir, "coordinates.csv"), row.names = FALSE)
      write.csv(state$processed$provenance, file.path(session_dir, "provenance_manifest.csv"), row.names = FALSE)
      write.csv(data.frame(tic_normalization = isTRUE(input$tic_normalize),
        log_transform = if (isTRUE(input$log1p_transform)) "log10(x + 1)" else "none",
        alignment_ppm = s$alignment_ppm, peak_pick_snr = s$peak_pick_snr,
        min_detection_fraction = s$min_detection_fraction), file.path(session_dir, "processing_parameters.csv"), row.names = FALSE)
    })
  })
  output$processing_gate <- renderUI(if (is.null(state$processed)) tags$div(class = "alert alert-warning", "No processed MSI is available yet.")
    else tags$div(class = "alert alert-success", "Processing completed. The full matrix remains in this R session; the browser receives previews only."))
  output$qc_table <- renderDT({ req(state$processed); datatable(state$processed$qc_summary, options = list(dom = "t", scrollX = TRUE)) })
  output$feature_preview <- renderDT({ req(state$processed); datatable(head(state$processed$feature_metadata, 25), options = list(pageLength = 10)) })
  output$tissue_mask_plot <- renderPlot({ req(state$tissue_gate); d <- state$tissue_gate$mask
    ggplot(d, aes(x, y, color = tissue_status)) + geom_point(size = .7) + coord_equal() + scale_y_reverse() +
      scale_color_manual(values = c("tissue" = "#2166AC", "unclassified/background" = "#D9D9D9")) + theme_minimal() +
      labs(title = "Full acquisition field and reproducible tissue gate") })
  output$tissue_mask_diagnostics <- renderDT({ req(state$tissue_gate); datatable(state$tissue_gate$diagnostics, options = list(dom = "t", scrollX = TRUE)) })
  output$tissue_download <- renderUI(if (!is.null(state$tissue_gate)) downloadButton("download_tissue_mask", "Download tissue mask") else NULL)
  output$download_tissue_mask <- downloadHandler(filename = function() "tissue_mask.csv", content = function(file) {
    req(state$tissue_gate); write.csv(state$tissue_gate$mask, file, row.names = FALSE) })

  output$registration_action <- renderUI({
    if (!isTRUE(state$valid) || !isTRUE(state$validation$capabilities$registration))
      return(tags$div(class = "alert alert-secondary", "Registration is disabled. Supply MSI, an optical JPEG and a compatible 3 × 3 transform JSON."))
    if (is.null(state$processed)) return(tags$div(class = "alert alert-warning", "Process MSI before registration."))
    actionButton("run_registration", "Apply supplied transform", class = "btn-primary")
  })
  observeEvent(input$run_registration, {
    req(state$processed); s <- state$spec
    state$registration <- register_metaspace_optical(state$processed$coordinates, s$transform_path,
      s$optical_path, if (present_path(s$attribution_path)) s$attribution_path else NULL, tic = state$tic,
      representative_ions = { m <- state$analysis_matrix; ii <- head(order(apply(m, 2, var), decreasing = TRUE), 3)
        stats::setNames(lapply(ii, function(j) m[, j]), format(state$processed$feature_metadata$mz[ii], digits = 9)) },
      output_dir = file.path(session_dir, "live_registration"))
    state$registration$diagnostics$diagnostic_source <- "live run using user-supplied optical image, transform and current MSI"
    write.csv(state$registration$registered_coordinates, file.path(session_dir, "registered_coordinates.csv"), row.names = FALSE)
  })
  output$registration_gate <- renderUI(if (is.null(state$registration)) NULL else tags$div(class = "alert alert-success",
    "Transform applied exactly as supplied; coordinate origin was inferred from current MSI minima and no implicit y inversion was introduced."))
  output$registration_diagnostics <- renderDT({ req(state$registration); datatable(state$registration$diagnostics, options = list(dom = "t")) })
  output$registration_plot <- renderImage({ req(state$registration); src <- unname(state$registration$output_files["measurement_mask_overlay"])
    list(src = src, contentType = "image/png", alt = "Live registered measurement mask") }, deleteFile = FALSE)

  observeEvent(input$domain_csv_runtime, {
    req(state$processed); info <- validate_label_file(input$domain_csv_runtime$datapath)
    state$domains <- resolve_label_mapping(info, state$processed$coordinates, state$spec$sample_id, state$spec$section_id)
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
    state$domains <- resolve_label_mapping(state$validation$labels, state$processed$coordinates, state$spec$sample_id, state$spec$section_id)
    state$domain_counts <- as.data.frame(table(state$domains$domain_label)); names(state$domain_counts) <- c("domain", "pixel_count")
  })
  observeEvent(input$generate_domains, {
    req(state$processed, state$analysis_matrix); keep_tissue <- if (is.null(state$tissue_mask)) rep(TRUE, nrow(state$analysis_matrix)) else state$tissue_mask
    matrix <- state$analysis_matrix[keep_tissue, , drop = FALSE]; keep <- apply(matrix, 2, var) > 0
    validate(need(sum(keep) >= 2, "At least two non-zero-variance features are required."))
    pc <- prcomp(matrix[, keep, drop = FALSE], center = TRUE, scale. = TRUE); npc <- min(input$domain_pcs, ncol(pc$x))
    set.seed(input$domain_seed); km <- kmeans(pc$x[, seq_len(npc), drop = FALSE], centers = input$domain_k, nstart = 25)
    ids <- rep(-1L, nrow(state$analysis_matrix)); ids[keep_tissue] <- km$cluster - 1L
    state$domains <- data.frame(state$processed$coordinates, tissue_status = ifelse(keep_tissue, "tissue", "unclassified/background"),
      domain_id = ids, domain_label = ifelse(ids == -1L, "unclassified/background", paste0("data-driven exploratory metabolic domain ", ids)))
    state$domain_counts <- as.data.frame(table(state$domains$domain_label)); names(state$domain_counts) <- c("domain", "pixel_count")
    between <- vapply(seq_len(ncol(matrix)), function(j) var(tapply(matrix[, j], km$cluster, mean)), numeric(1)); order <- head(order(between, decreasing = TRUE), 25)
    state$domain_features <- data.frame(feature = state$processed$feature_metadata$column_name[order],
      mz = state$processed$feature_metadata$mz[order], between_domain_variance = between[order])
    write.csv(state$domains, file.path(session_dir, "exploratory_metabolic_domains.csv"), row.names = FALSE)
    write.csv(data.frame(k = input$domain_k, seed = input$domain_seed, n_pcs = npc,
      inference = "descriptive_only_no_pixel_replicates"), file.path(session_dir, "domain_parameters.csv"), row.names = FALSE)
  })
  output$domain_gate <- renderUI(if (is.null(state$domains)) tags$div(class = "alert alert-warning", "No labels or generated domains are currently loaded.") else NULL)
  output$domain_plot <- renderPlot({ req(state$domains); ggplot(state$domains, aes(x, y, color = domain_label)) + geom_point(size = .7) + coord_equal() + scale_y_reverse() + theme_minimal() })
  output$domain_counts <- renderDT({ req(state$domain_counts); datatable(state$domain_counts, options = list(dom = "t")) })
  output$domain_features <- renderDT({ req(state$domain_features); datatable(state$domain_features, options = list(pageLength = 10)) })
  output$domain_download <- renderUI(if (!is.null(state$domains)) downloadButton("download_domains", "Download domain CSV") else NULL)
  output$download_domains <- downloadHandler(filename = function() "spatial_domains.csv", content = function(file) { req(state$domains); write.csv(state$domains, file, row.names = FALSE) })

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
    actionButton("run_moran", "Run Moran's I", class = "btn-primary")
  })
  observeEvent(input$run_moran, {
    req(state$processed, state$analysis_matrix); keep <- if (is.null(state$tissue_mask)) rep(TRUE, nrow(state$analysis_matrix)) else state$tissue_mask
    coordinates <- state$processed$coordinates[keep, , drop = FALSE]; matrix <- state$analysis_matrix[keep, , drop = FALSE]
    threshold <- if (input$neighbor_method == "queen") sqrt(2) else 1
    state$neighbors <- build_spatial_neighbors(coordinates$x, coordinates$y, input$neighbor_method, threshold)
    pixel <- data.frame(coordinates[c("pixel_id", "sample_id", "section_id", "x", "y")], matrix, check.names = FALSE)
    names(pixel)[-(1:5)] <- state$processed$feature_metadata$column_name
    state$moran <- compute_spatially_variable_metabolites(pixel, coordinates, n_perm = input$permutations,
      alternative = "two.sided", p_adjust_method = "BH", seed = 20260807, neighbors = state$neighbors)
    write.csv(state$moran, file.path(session_dir, paste0("moran_", input$neighbor_method, ".csv")), row.names = FALSE)
  })
  observeEvent(input$load_moran, {
    validate(need(!is.null(state$processed), "Process the current MSI before attaching a previous result."))
    directory <- normalizePath(input$moran_result_dir, mustWork = FALSE)
    path <- file.path(directory, paste0("spatially_variable_metabolites_", input$neighbor_method, ".csv"))
    validate(need(file.exists(path), "The selected adjacency result CSV was not found.")); state$moran <- read.csv(path, check.names = FALSE)
  })
  output$moran_gate <- renderUI(if (is.null(state$moran)) tags$div(class = "alert alert-warning", "No Moran result is loaded for the current MSI.") else NULL)
  output$neighbor_table <- renderDT({ req(state$neighbors); datatable(spatial_neighbor_diagnostics(state$neighbors), options = list(pageLength = 10)) })
  output$moran_table <- renderDT({ req(state$moran); datatable(head(state$moran[order(state$moran$adj_p_value, -abs(state$moran$morans_i)), ], 100), options = list(pageLength = 10)) })
  output$ion_plot <- renderPlot({ req(state$moran, state$processed); i <- which.max(state$moran$morans_i); column <- state$moran$feature[i]
    d <- data.frame(state$processed$coordinates, intensity = state$analysis_matrix[, match(column, state$processed$feature_metadata$column_name)])
    ggplot(d, aes(x, y, color = intensity)) + geom_point(size = .7) + coord_equal() + scale_y_reverse() + scale_color_viridis_c() + theme_minimal() + labs(title = paste("Top Moran feature", column)) })

  output$lcms_action <- renderUI(if (isTRUE(state$valid) && isTRUE(state$validation$capabilities$lcms))
    actionButton("load_lcms", "Read precursor window", class = "btn-primary") else tags$div(class = "alert alert-secondary", "LC-MS/MS evidence is disabled. Supply and validate an mzML file."))
  observeEvent(input$load_lcms, {
    req(state$valid); withProgress(message = "Reading target MS2 window", value = .1, {
      state$lcms <- read_mzml_fragment_spectra(state$spec$lcms_path, input$lcms_precursor, 10); incProgress(.9)
      if (nrow(state$lcms$precursor_scan_metadata)) {
        state$lcms$precursor_scan_metadata$msi_to_reference_ppm <- (input$msi_mz - input$lcms_precursor) / input$lcms_precursor * 1e6
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
    req(input$msi_feature_csv, input$lcms_feature_csv)
    msi <- read.csv(input$msi_feature_csv$datapath, check.names = FALSE, stringsAsFactors = FALSE)
    lcms <- read.csv(input$lcms_feature_csv$datapath, check.names = FALSE, stringsAsFactors = FALSE)
    require_columns(msi, "mz", "MSI feature CSV")
    require_columns(lcms, "mz", "LC-MS feature CSV")
    state$feature_matches <- cross_validate_msi_lcms(msi, lcms, ppm = input$feature_match_ppm,
      assignment_method = input$assignment_method)
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
    req(input$observed_ccs_csv, input$reference_ccs_csv)
    observed <- read.csv(input$observed_ccs_csv$datapath, check.names = FALSE, stringsAsFactors = FALSE)
    reference <- read.csv(input$reference_ccs_csv$datapath, check.names = FALSE, stringsAsFactors = FALSE)
    state$ccs_evidence <- validate_ccs_evidence(observed, reference, ccs_tolerance_pct = input$ccs_tolerance)
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
