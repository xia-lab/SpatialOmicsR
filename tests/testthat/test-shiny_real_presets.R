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

launcher_formals <- formals(SpatialOmicsMSI::run_spatial_app)
stopifnot(all(c("launch.browser", "port", "host") %in% names(launcher_formals)),
  identical(eval(launcher_formals$host), "127.0.0.1"),
  inherits(try(SpatialOmicsMSI::run_spatial_app(host = ""), silent = TRUE), "try-error"),
  inherits(try(SpatialOmicsMSI::run_spatial_app(port = 70000), silent = TRUE), "try-error"))

stopifnot(all(vapply(c("Spatial MSI analysis", "Spatial MSI", "LC-MS/MS evidence (optional)",
  "Choose an example dataset",
  "Paired imzML + ibd", "Processed pixel × feature CSV", "Biological subject ID",
  "Optional optical/H&E/brightfield JPEG", "Optional ROI/domain CSV", "LC-MS/MS mzML",
  "Continue to MSI processing", "Technical details", "Pseudoreplication warning",
  "Optional advanced check", "second coordinate-aligned map", "MSI–LC-MS feature-table matching",
  "Optimal (recommended)", "CCS / ion-mobility evidence", "Spatial niches",
  "Hellinger (recommended)", "Niche assignment ambiguity",
  "ROI/subregion summaries", "Import statistical result and map it back",
  "Results & export", "Choose how much tissue enters domain detection",
  "Complete tissue (recommended; choose domains afterward)",
  "Working…", "The page will update automatically when this step finishes.",
  "Use selected domains as downstream ROI",
  "MetaboAnalyst lipid enrichment export", "Prepare lipid export",
  "View Brain01 QC example", "View OMIX016317 QC example",
  "MSV000090179 external-reference example", "MSIflow study-matched evidence",
  "Kasarla 2025 kidney ROI validation", "Module-only LC-MS/MS evidence",
  "Kasarla 2025 kidney MSI + LMD-LC-MS/MS study bundle", "Load selected example"),
  grepl, logical(1), x = text, fixed = TRUE)))
stopifnot(all(vapply(c("1. Input & provenance", "2. Processing & QC",
  "3. Registration (optional)", "4. Optional pre-analysis ROI", "5. Domains & niches",
  "6. ROI summaries & statistics", "7. Spatial validation",
  "8. Matched transcriptomics & H&E", "9. LC-MS/MS evidence (optional)",
  "10. Results & export",
  "Vicari 2024", "Load selected example"), grepl, logical(1), x = text, fixed = TRUE)))
stopifnot(!grepl("Start LC-MS/MS", text, fixed = TRUE),
  !grepl("Start matched analysis", text, fixed = TRUE),
  !grepl("Matched integration", text, fixed = TRUE),
  !grepl('selectInput("example_key"', text, fixed = TRUE))
stopifnot(!grepl("Registration is available for Brain01 only", text, fixed = TRUE))
stopifnot(!grepl("Select the MSV000090179 preset first", text, fixed = TRUE))
stopifnot(!grepl("manual-busy", text, fixed = TRUE))
stopifnot(!grepl("shiny-busy-overlay", text, fixed = TRUE))
stopifnot(grepl("server-busy-card", text, fixed = TRUE))

env <- new.env(parent = globalenv())
old <- Sys.getenv("SPATIALOMICS_PROJECT_ROOT", unset = NA)
test_project_root <- Sys.getenv("SPATIALOMICS_PROJECT_ROOT", unset = "")
if (!file.exists(file.path(test_project_root, "DESCRIPTION"))) {
  candidate <- normalizePath(file.path(getwd(), "../.."), mustWork = FALSE)
  test_project_root <- if (file.exists(file.path(candidate, "DESCRIPTION"))) candidate else normalizePath(getwd())
}
Sys.setenv(SPATIALOMICS_PROJECT_ROOT = test_project_root)
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
  lcms_relationship = "none",
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
with_modules$lcms_relationship <- "external_reference"
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
stopifnot(identical(env$example_workflow(env$example_datasets$brain01), "spatial"),
  identical(env$example_workflow(env$example_datasets$omix), "spatial"),
  identical(env$example_workflow(env$example_datasets$lcms), "lcms"),
  identical(env$example_workflow(env$example_datasets$msiflow), "matched"))
stopifnot(setequal(names(env$spatial_example_datasets), c("brain01", "omix", "msiflow",
  "vicari_mpd1", "vicari_mpd3", "vicari_mpd4", "kasarla_kidney")))
stopifnot(setequal(unname(unlist(env$example_dataset_choices)), names(env$example_datasets)),
  "lcms" %in% env$example_dataset_choices[["Module-only LC-MS/MS evidence"]])
for (key in c("vicari_mpd1", "vicari_mpd3", "vicari_mpd4")) {
  p <- env$preset_spec(key)
  vicari_paths <- unlist(p[c("csv_path", "st_expression_path", "st_positions_path",
    "histology_path", "st_region_path", "st_lesion_path", "st_loupe_path")], use.names = FALSE)
  stopifnot(all(nzchar(vicari_paths)))
  # External example data are present in the development workspace, but are
  # intentionally excluded from the source package and CRAN checks.
  if (dir.exists(file.path(getwd(), "data_raw", "vicari_sma"))) {
    stopifnot(all(file.exists(vicari_paths)))
  }
}
for (key in names(env$spatial_example_datasets)) {
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
kasarla_spec <- env$preset_spec("kasarla_kidney")
if (file.exists(kasarla_spec$csv_path)) {
  kasarla_loaded <- env$read_pixel_feature_csv(kasarla_spec$csv_path,
    kasarla_spec$sample_id, kasarla_spec$section_id, kasarla_spec$subject_id)
  stopifnot(nrow(kasarla_loaded$pixel_feature_matrix) == 45623L,
    nrow(kasarla_loaded$feature_metadata) == 236L,
    identical(kasarla_spec$lcms_relationship, "study_matched"))
}

shiny::testServer(env$server, {
  session$setInputs(data_origin = "server", msi_input_type = "csv", server_csv = csv,
    server_imzml = "", server_ibd = "", server_optical = "", server_transform = "",
    server_labels = "", server_lcms = "", server_attribution = "",
    sample_id = "sample_custom", section_id = "section_custom", subject_id = "subject_custom",
    ion_mode = "positive", ion_mode_source = "user_acquisition_record",
    lcms_relationship = "external_reference",
    spectral_processing = "processed_peak_lists", alignment_ppm = 7,
    peak_pick_snr = 5, min_detection_fraction = 0, make_tissue_mask = FALSE,
    tic_normalize = TRUE, log1p_transform = TRUE, validate_input = 1)
  session$flushReact()
  stopifnot(state$valid, state$validation$capabilities$processing)
  session$setInputs(run_processing = 1); session$flushReact()
  stopifnot(nrow(state$processed$pixel_feature_matrix) == 4L,
            ncol(state$analysis_matrix) == 2L, all(state$tissue_mask))
  session$setInputs(open_kasarla_validation = 1); session$flushReact()
  stopifnot(!is.null(state$kasarla_validation), nrow(state$kasarla_validation) == 57L,
    !is.null(state$kasarla_bundle), nrow(state$kasarla_bundle) == 13L,
    sum(state$kasarla_bundle$role == "LMD-LC-MS/MS") == 12L,
    identical(state$kasarla_bundle$dataset_id[1], "2023-12-19_13h02m39s"),
    all(state$kasarla_validation$n_pairs == 3L),
    all(is.finite(state$kasarla_validation$correlation)),
    sum(state$kasarla_validation$correlation > 0) >= 51L)
  session$setInputs(neighbor_method = "queen", permutations = 99, run_moran = 1)
  session$flushReact()
  stopifnot(!is.null(state$moran), nrow(state$moran) >= 1L,
    all(c("feature", "morans_i", "adj_p_value") %in% names(state$moran)))
  state$metabo_result <- data.frame(Feature = state$moran$feature[1], p.value = 0.01,
    stringsAsFactors = FALSE)
  state$metabo_result_type <- "differential"
  concordance <- moran_concordance()
  stopifnot(nrow(concordance) == 1L,
    all(c("morans_i", "spatial_coherence_class") %in% names(concordance)))
  state$metabo_result <- NULL; state$metabo_result_type <- NULL
  # Domain detection can start on the complete tissue without a prior ROI, and
  # selected data-driven domains can subsequently become the downstream ROI.
  stopifnot(is.null(state$rois))
  session$setInputs(domain_k = 2, domain_pcs = 2, domain_seed = 19,
    domain_feature_source = "all", domain_moran_top_n = 2,
    generate_domains = 1)
  session$flushReact()
  stopifnot(!is.null(state$domains),
    sum(state$domains$domain_label != "unclassified/background") == 4L)
  detected_labels <- unique(state$domains$domain_label)
  detected_labels <- detected_labels[detected_labels != "unclassified/background"]
  session$setInputs(domain_compare_labels = detected_labels[1:2]); session$flushReact()
  domain_pair <- domain_comparison()
  stopifnot(nrow(domain_pair$contact) == 1L, nrow(domain_pair$features) == 2L,
    identical(domain_pair$selected, detected_labels[1:2]))
  session$setInputs(domain_roi_labels = detected_labels[1],
    domain_roi_mode = "separate", use_domains_as_roi = 1)
  session$flushReact()
  stopifnot(identical(state$roi_source, "data_driven_domains"),
    any(!is.na(state$rois$roi_id)), !is.null(state$domains))

  roi_path <- tempfile(fileext = ".csv")
  write.csv(data.frame(
    roi_id = c("left", "right"),
    x_min = c(1, 2), x_max = c(1, 2),
    y_min = c(1, 1), y_max = c(2, 2)
  ), roi_path, row.names = FALSE)
  session$setInputs(roi_table_upload = list(datapath = roi_path, name = "two_rois.csv"))
  session$flushReact()
  stopifnot(!is.null(state$rois),
    identical(sort(unique(na.omit(state$rois$roi_id))), c("left", "right")))
  session$setInputs(domain_csv_runtime = list(datapath = labels))
  session$flushReact()
  stopifnot(nrow(state$domains) == 4L)
  session$setInputs(
    niche_neighborhood = "knn", niche_k = 1, niche_min_neighbors = 1,
    niche_count = 2, niche_transform = "hellinger", niche_include_self = TRUE,
    niche_domain_alignment = "joint", niche_seed = 31, run_niches = 1
  )
  session$flushReact()
  stopifnot(!is.null(state$niches), nrow(state$niches$matrix) == 3L,
    all(state$niches$matrix$eligible),
    length(unique(state$niches$matrix$niche_id)) == 2L,
    identical(state$niches$settings$transform, "hellinger"))
  niche_labels <- sort(unique(state$niches$matrix$niche_label))
  session$setInputs(niche_compare_labels = niche_labels[1:2]); session$flushReact()
  niche_pair <- niche_comparison()
  stopifnot(nrow(niche_pair$contact) == 1L, nrow(niche_pair$composition) == 2L,
    nrow(niche_pair$features) == 2L)

  session$setInputs(
    sampling_grid_size = 2, sampling_min_pixels = 1,
    sampling_grid_scope = "roi", run_sampling = 1
  )
  session$flushReact()
  stopifnot(!is.null(state$sampling), nrow(state$sampling$sample_matrix) >= 2L,
    nrow(state$metabo_input) == nrow(state$sampling$sample_matrix),
    all(c("Sample", "Group") %in% names(state$metabo_input)),
    !is.null(state$functional_peak_table),
    identical(state$functional_peak_table$feature[1:2], c("Sample", "Group")),
    all(is.finite(as.numeric(state$functional_peak_table$feature[-(1:2)]))),
    !is.null(state$functional_peak_list),
    identical(names(state$functional_peak_list), c("m.z", "p.value", "t.score")),
    all(vapply(state$functional_peak_list, function(x) all(is.finite(x)), logical(1))),
    !any(grepl("^grid_", names(state$functional_peak_list))))
  lipid_annotation_path <- tempfile(fileext = ".csv")
  write.csv(data.frame(
    lipid_name = "PC(16:0/18:1)",
    mz = state$functional_peak_list$m.z[1],
    lipid_id = "LMGP01010001"
  ), lipid_annotation_path, row.names = FALSE)
  session$setInputs(
    lipid_annotation_csv = list(datapath = lipid_annotation_path, name = "lipids.csv"),
    lipid_match_ppm = 5, lipid_p_cutoff = 1
  )
  session$flushReact()
  session$setInputs(prepare_lipid_export = 1)
  session$flushReact()
  stopifnot(nrow(state$lipid_export) == 1L,
    identical(state$lipid_name_list, "PC(16:0/18:1)"),
    all(state$lipid_export$selected_for_enrichment),
    all(abs(state$lipid_export$ppm_error) < sqrt(.Machine$double.eps)))
  score_path <- tempfile(fileext = ".csv")
  write.csv(data.frame(
    Sample = state$sampling$sample_matrix$sample_id,
    PC1 = seq_len(nrow(state$sampling$sample_matrix)),
    PC2 = rev(seq_len(nrow(state$sampling$sample_matrix)))
  ), score_path, row.names = FALSE)
  session$setInputs(metabo_result_csv = list(datapath = score_path, name = "pca_scores.csv"))
  session$flushReact()
  session$setInputs(metabo_score_column = "PC1")
  session$flushReact()
  stopifnot(identical(state$metabo_result_type, "pca_scores"),
    !is.null(state$mapped_scores), nrow(state$mapped_scores) == 4L,
    all(is.finite(state$mapped_scores$score)))

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
    all(state$feature_matches$direction_agreement),
    all(state$feature_matches$evidence_relationship == "external_reference"))

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
  stopifnot(nrow(state$ccs_evidence) == 1L,
    state$ccs_evidence$evidence_relationship == "external_reference")
  session$setInputs(use_brain01_qc = 1)
  session$flushReact()
  stopifnot(identical(state$qc_example$key, "brain01"), state$valid)
})

# A matched Vicari preset enters through the same workflow and retains every
# section-specific Visium/H&E companion for later modules. The real-data branch
# runs only in the development workspace because these large external files are
# intentionally absent from the source package.
vicari_mpd1 <- env$preset_spec("vicari_mpd1")
vicari_inputs <- unlist(vicari_mpd1[c("csv_path", "st_expression_path",
  "st_positions_path", "histology_path")], use.names = FALSE)
if (all(file.exists(vicari_inputs))) shiny::testServer(env$server, {
  session$setInputs(example_dataset = "vicari_mpd1", load_example_dataset = 1)
  session$flushReact()
  stopifnot(state$valid, identical(state$preset_key, "vicari_mpd1"),
    identical(state$spec$sample_id, "mPD1"),
    identical(state$spec$section_id, "V11L12-038_D1"),
    length(matched_paths()) == 6L,
    all(file.exists(matched_paths())))
})

# The Kasarla study bundle starts at Step 1 with a pixel-level MSI matrix and
# retains its matched ROI-level MALDI/LMD evidence for Step 9.
kasarla_inputs <- env$preset_spec("kasarla_kidney")$csv_path
if (file.exists(kasarla_inputs)) shiny::testServer(env$server, {
  session$setInputs(example_dataset = "kasarla_kidney", load_example_dataset = 1)
  session$flushReact()
  stopifnot(state$valid, identical(state$preset_key, "kasarla_kidney"),
    nrow(state$validation$metadata$coordinates) == 45623L,
    isTRUE(state$kasarla_pixel_msi_available),
    !is.null(state$kasarla_validation), nrow(state$kasarla_validation) == 57L,
    !is.null(state$kasarla_bundle), nrow(state$kasarla_bundle) == 13L)
  session$setInputs(tic_normalize = TRUE, log1p_transform = TRUE, run_processing = 1)
  session$flushReact()
  stopifnot(nrow(state$processed$pixel_feature_matrix) == 45623L,
    ncol(state$analysis_matrix) == 236L, all(state$tissue_mask))
})

# The same selector also exposes module-only evidence without pretending that
# an LC-MS/MS file is a pixel-level MSI input. It opens the function it can run.
msv_input <- env$preset_spec("lcms")$lcms_path
if (file.exists(msv_input)) shiny::testServer(env$server, {
  session$setInputs(example_dataset = "lcms", load_example_dataset = 1)
  session$flushReact()
  stopifnot(identical(state$preset_key, "lcms"), !state$valid,
    identical(state$lcms_example$example_key, "msv"),
    identical(state$lcms_example$relationship, "external_reference"),
    nrow(state$lcms_example$precursor_scan_metadata) > 0L)
})

# Module-local examples render directly without replacing the active analysis.
module_example_files <- c(
  env$example_datasets$msiflow$csv_path,
  env$example_datasets$msiflow$labels_path,
  env$example_datasets$msiflow$lcms_path,
  env$example_datasets$brain01$imzml_path,
  env$example_datasets$brain01$ibd_path,
  env$example_datasets$brain01$optical_path,
  env$example_datasets$brain01$transform_path
)
if (all(file.exists(module_example_files))) shiny::testServer(env$server, {
  session$setInputs(use_msiflow_evidence = 1); session$flushReact()
  stopifnot(identical(state$lcms_example$relationship, "study_matched"),
    nrow(state$lcms_example$precursor_scan_metadata) > 0L)
  session$setInputs(use_brain01_qc = 1); session$flushReact()
  session$setInputs(use_brain01_registration = 1); session$flushReact()
  session$setInputs(use_msiflow_domains = 1); session$flushReact()
  stopifnot(identical(state$qc_example$key, "brain01"),
    identical(state$registration_example$key, "brain01"),
    identical(state$domain_example$key, "msiflow"),
    !is.null(state$processed), !is.null(state$domains),
    nrow(state$domains) == nrow(state$processed$coordinates))
})

# LC-MS/MS alone cannot start the MSI workflow.
shiny::testServer(env$server, {
  session$setInputs(data_origin = "server", msi_input_type = "imzml",
    server_imzml = "", server_ibd = "", server_csv = "",
    server_optical = "", server_transform = "", server_labels = "",
    server_lcms = lcms, server_attribution = "", sample_id = "", section_id = "",
    subject_id = "", ion_mode = "positive", ion_mode_source = "MassIVE metadata",
    lcms_relationship = "external_reference",
    spectral_processing = "processed_peak_lists", alignment_ppm = 10,
    peak_pick_snr = 6, min_detection_fraction = 0, make_tissue_mask = FALSE,
    validate_input = 1)
  session$flushReact()
  stopifnot(!state$valid, !state$validation$capabilities$msi)
})

# MSI plus study-level LC-MS/MS evidence remains a spatial MSI analysis.
shiny::testServer(env$server, {
  session$setInputs(msi_input_type = "csv", data_origin = "server", server_csv = csv,
    server_imzml = "", server_ibd = "",
    server_optical = "", server_transform = "", server_labels = "",
    server_lcms = lcms, server_attribution = "", sample_id = "matched_sample",
    section_id = "matched_section", subject_id = "matched_subject", ion_mode = "positive",
    ion_mode_source = "acquisition metadata", lcms_relationship = "study_matched",
    spectral_processing = "processed_peak_lists",
    alignment_ppm = 10, peak_pick_snr = 6, min_detection_fraction = 0,
    make_tissue_mask = FALSE, validate_input = 1)
  session$flushReact()
  stopifnot(state$valid, state$validation$capabilities$msi,
    state$validation$capabilities$lcms)
})
if (is.na(old)) Sys.unsetenv("SPATIALOMICS_PROJECT_ROOT") else Sys.setenv(SPATIALOMICS_PROJECT_ROOT = old)
cat("SHINY_GENERAL_APP_TEST_OK=TRUE\n")

  testthat::succeed()
})
