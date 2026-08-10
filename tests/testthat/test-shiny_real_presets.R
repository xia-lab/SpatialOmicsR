testthat::test_that("regression: test-shiny_real_presets", {
if (dir.exists(".lib")) .libPaths(c(normalizePath(".lib"), .libPaths()))
if (!requireNamespace("shiny", quietly = TRUE)) { cat("SHINY_GENERAL_APP_SKIPPED=TRUE\n"); quit(status = 0) }
checking_installed_package <- nzchar(Sys.getenv("_R_CHECK_PACKAGE_NAME_"))
if (!checking_installed_package && !requireNamespace("SpatialOmicsMSI", quietly = TRUE)) {
  cat("SHINY_SOURCE_APP_SKIPPED_PACKAGE_NOT_INSTALLED=TRUE\n")
  quit(status = 0)
}
app_path <- if (checking_installed_package) {
  system.file("shiny", "spatial_pipeline", "app.R", package = "SpatialOmicsMSI")
} else file.path(normalizePath(getwd()), "inst", "shiny", "spatial_pipeline", "app.R")
text <- paste(readLines(app_path, warn = FALSE), collapse = "\n")

stopifnot(all(vapply(c("Start a new analysis", "Load an example dataset", "Example datasets",
  "Paired imzML + ibd", "Processed pixel × feature CSV", "Biological subject ID",
  "Optional optical/H&E/brightfield JPEG", "Optional ROI/domain CSV", "Optional LC-MS/MS mzML",
  "Continue to Processing", "Technical details", "Pseudoreplication warning",
  "Cross-modal ROI validation", "MSI–LC-MS feature-table matching",
  "Optimal (recommended)", "CCS / ion-mobility evidence"),
  grepl, logical(1), x = text, fixed = TRUE)))
stopifnot(!grepl("Registration is available for Brain01 only", text, fixed = TRUE))
stopifnot(!grepl("Select the MSV000090179 preset first", text, fixed = TRUE))

env <- new.env(parent = globalenv())
old <- Sys.getenv("SPATIALOMICS_PROJECT_ROOT", unset = NA)
Sys.setenv(SPATIALOMICS_PROJECT_ROOT = normalizePath(getwd()))
sys.source(app_path, envir = env)

# A general processed CSV follows exactly the same validation and loading path.
csv <- tempfile(fileext = ".csv")
write.csv(data.frame(x = c(1, 2, 1, 2), y = c(1, 1, 2, 2),
  `100.123456` = c(1, 2, 3, 4), arbitrary_feature = c(4, 3, 2, 1), check.names = FALSE),
  csv, row.names = FALSE)
base <- list(msi_input_type = "csv", imzml_path = "", ibd_path = "", csv_path = csv,
  optical_path = "", transform_path = "", labels_path = "", lcms_path = "", attribution_path = "",
  sample_id = "sample_custom", section_id = "section_custom", subject_id = "subject_custom",
  ion_mode = "positive", ion_mode_source = "user_acquisition_record",
  spectral_processing = "processed_peak_lists", alignment_ppm = 7,
  peak_pick_snr = 5, min_detection_fraction = 0, make_tissue_mask = FALSE)
validation <- env$validate_analysis_spec(base, inspect = FALSE)
loaded <- env$read_pixel_feature_csv(csv, base$sample_id, base$section_id, base$subject_id)
stopifnot(validation$valid, validation$capabilities$processing, validation$capabilities$moran,
  !validation$capabilities$registration, !validation$capabilities$lcms,
  nrow(loaded$pixel_feature_matrix) == 4L, nrow(loaded$feature_metadata) == 2L,
  loaded$coordinates$subject_id[1] == "subject_custom")

# Optional modules are capability-gated from supplied files, not preset identity.
optical <- tempfile(fileext = ".jpg"); jpeg::writeJPEG(array(.8, c(12, 12, 3)), optical)
transform <- tempfile(fileext = ".json")
writeLines('{"transform":[[1,0,0],[0,1,0],[0,0,1]]}', transform)
labels <- tempfile(fileext = ".csv")
write.csv(data.frame(x = c(1, 2, 1, 2), y = c(1, 1, 2, 2), domain_id = c(-1, 0, 1, 1)), labels, row.names = FALSE)
lcms <- tempfile(fileext = ".mzML"); writeLines('<mzML><run><spectrumList count="0"/></run></mzML>', lcms)
with_modules <- base; with_modules$optical_path <- optical; with_modules$transform_path <- transform
with_modules$labels_path <- labels; with_modules$lcms_path <- lcms
module_validation <- env$validate_analysis_spec(with_modules, inspect = FALSE)
stopifnot(module_validation$valid, all(unlist(module_validation$capabilities)))

# Explicit imzML/ibd pairing must not depend on a shared basename.
pair_dir <- tempfile("nonfixed_pair_"); dir.create(pair_dir)
imz <- file.path(pair_dir, "user-selected-section.imzML")
ibd <- file.path(pair_dir, "binary-from-instrument.ibd")
writeLines("<mzML><fileDescription/></mzML>", imz)
writeBin(as.raw(seq_len(16)), ibd)
nonfixed <- base; nonfixed$msi_input_type <- "imzml"; nonfixed$csv_path <- ""
nonfixed$imzml_path <- imz; nonfixed$ibd_path <- ibd
stopifnot(env$validate_analysis_spec(nonfixed, inspect = FALSE)$valid)
readable <- env$cardinal_pair_path(validate_imzml_ibd_pair(imz, ibd))
stopifnot(readable$temporary, basename(readable$imzml_path) == "paired.imzML",
          file.exists(file.path(dirname(readable$imzml_path), "paired.ibd")))

# All examples populate the same analysis specification and validation function.
for (key in names(env$example_datasets)) {
  p <- env$preset_spec(key)
  stopifnot(!"kind" %in% names(p), p$msi_input_type %in% c("imzml", "csv", "none"),
            nzchar(p$sample_id), nzchar(p$section_id), nzchar(p$subject_id),
            p$ion_mode %in% c("positive", "negative"), nzchar(p$ion_mode_source))
  spec <- base
  for (name in names(p)) spec[[name]] <- p[[name]]
  for (name in c("imzml_path", "ibd_path", "csv_path", "optical_path", "transform_path",
                 "labels_path", "lcms_path", "attribution_path")) if (is.null(spec[[name]])) spec[[name]] <- ""
  required <- if (spec$msi_input_type == "imzml") c(spec$imzml_path, spec$ibd_path) else if (spec$msi_input_type == "csv") spec$csv_path else spec$lcms_path
  if (all(file.exists(required))) {
    result <- env$validate_analysis_spec(spec, inspect = FALSE)
    stopifnot(result$valid)
  }
}

shiny::testServer(env$server, {
  stopifnot(input$entry == "new")
  session$setInputs(data_origin = "server", msi_input_type = "csv", server_csv = csv,
    server_imzml = "", server_ibd = "", server_optical = "", server_transform = "",
    server_labels = "", server_lcms = "", server_attribution = "",
    sample_id = "sample_custom", section_id = "section_custom", subject_id = "subject_custom",
    ion_mode = "positive", ion_mode_source = "user_acquisition_record",
    spectral_processing = "processed_peak_lists", alignment_ppm = 7,
    peak_pick_snr = 5, min_detection_fraction = 0, make_tissue_mask = FALSE,
    tic_normalize = TRUE, log1p_transform = TRUE, validate_input = 1)
  session$flushReact()
  stopifnot(state$valid, state$validation$capabilities$processing)
  session$setInputs(run_processing = 1); session$flushReact()
  stopifnot(nrow(state$processed$pixel_feature_matrix) == 4L,
            ncol(state$analysis_matrix) == 2L, all(state$tissue_mask))

  msi_features <- tempfile(fileext = ".csv")
  lcms_features <- tempfile(fileext = ".csv")
  write.csv(data.frame(feature_id = c("m1", "m2"), mz = c(100, 200),
    ion_mode = "positive", log2fc = c(1, -1)), msi_features, row.names = FALSE)
  write.csv(data.frame(feature_id = c("l1", "l2"), mz = c(100.0002, 200.0003),
    ion_mode = "positive", log2fc = c(2, -2)), lcms_features, row.names = FALSE)
  session$setInputs(msi_feature_csv = list(datapath = msi_features),
    lcms_feature_csv = list(datapath = lcms_features), feature_match_ppm = 5,
    assignment_method = "optimal", run_feature_matching = 1)
  session$flushReact()
  stopifnot(nrow(state$feature_matches) == 2L,
    all(state$feature_matches$assignment_method == "optimal"),
    all(state$feature_matches$direction_agreement))

  observed_ccs <- tempfile(fileext = ".csv")
  reference_ccs <- tempfile(fileext = ".csv")
  write.csv(data.frame(candidate_id = "c1", observed_ccs = 150,
    observed_source = "lcms_empirical"), observed_ccs, row.names = FALSE)
  write.csv(data.frame(candidate_id = "c1", reference_ccs = 151,
    reference_source = "measured_library"), reference_ccs, row.names = FALSE)
  session$setInputs(observed_ccs_csv = list(datapath = observed_ccs),
    reference_ccs_csv = list(datapath = reference_ccs), ccs_tolerance = 2,
    run_ccs_validation = 1)
  session$flushReact()
  stopifnot(nrow(state$ccs_evidence) == 1L)
  session$setInputs(entry = "example", example_key = "brain01", load_example = 1)
  session$flushReact()
  stopifnot(input$entry == "example")
})
if (is.na(old)) Sys.unsetenv("SPATIALOMICS_PROJECT_ROOT") else Sys.setenv(SPATIALOMICS_PROJECT_ROOT = old)
cat("SHINY_GENERAL_APP_TEST_OK=TRUE\n")

  testthat::succeed()
})
