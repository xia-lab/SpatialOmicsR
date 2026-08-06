checking_installed_package <- nzchar(Sys.getenv("_R_CHECK_PACKAGE_NAME_"))
if (checking_installed_package) {
  library(SpatialOmicsMSI)
  orient_shared_intensity_matrix <- getFromNamespace("orient_shared_intensity_matrix", "SpatialOmicsMSI")
  stable_mz_column_names <- getFromNamespace("stable_mz_column_names", "SpatialOmicsMSI")
  project_root <- NULL
} else {
  project_root <- normalizePath(getwd(), mustWork = TRUE)
  source(file.path(project_root, "R", "msi_pipeline.R"))
}

expect_error <- function(expression, pattern = NULL) {
  error <- tryCatch(expression, error = function(condition) condition)
  stopifnot(inherits(error, "error"))
  if (!is.null(pattern)) stopifnot(grepl(pattern, conditionMessage(error), fixed = TRUE))
  invisible(error)
}

# Both supported intensity orientations produce spectra x features output.
feature_by_spectrum <- matrix(seq_len(12), nrow = 3, ncol = 4)
oriented_a <- orient_shared_intensity_matrix(feature_by_spectrum, 3, 4)
stopifnot(identical(dim(oriented_a$values), c(4L, 3L)))
stopifnot(identical(oriented_a$source_orientation, "features_x_spectra"))
spectrum_by_feature <- t(feature_by_spectrum)
oriented_b <- orient_shared_intensity_matrix(spectrum_by_feature, 3, 4)
stopifnot(identical(dim(oriented_b$values), c(4L, 3L)))
stopifnot(identical(oriented_b$source_orientation, "spectra_x_features"))
stopifnot(identical(oriented_a$values, oriented_b$values))

# Stable names are unique while scientific values remain separate numeric metadata.
precise_mz <- c(311.26712036132812, 311.26712036132812, 912.68975830078125)
stable_names <- stable_mz_column_names(precise_mz)
stopifnot(!anyDuplicated(stable_names))
stopifnot(identical(precise_mz, as.numeric(precise_mz)))

# Missing ibd is rejected before attempting to read imzML.
missing_pair_dir <- tempfile("missing_ibd_")
dir.create(missing_pair_dir)
missing_imzml <- file.path(missing_pair_dir, "missing.imzML")
invisible(file.create(missing_imzml))
expect_error(
  load_centroided_msi_features(
    missing_imzml,
    sample_id = "sample",
    section_id = "section",
    ion_mode = "positive",
    ion_mode_source = "test"
  ),
  "companion .ibd file does not exist"
)

# Polarity and its source must be explicit, even when an imzML/ibd pair exists.
fake_pair_dir <- tempfile("explicit_polarity_")
dir.create(fake_pair_dir)
fake_imzml <- file.path(fake_pair_dir, "fake.imzML")
fake_ibd <- file.path(fake_pair_dir, "fake.ibd")
invisible(file.create(fake_imzml, fake_ibd))
expect_error(
  load_centroided_msi_features(
    fake_imzml,
    sample_id = "sample",
    section_id = "section",
    ion_mode = NULL,
    ion_mode_source = "test"
  ),
  "ion_mode must be supplied explicitly"
)
expect_error(
  load_centroided_msi_features(
    fake_imzml,
    sample_id = "sample",
    section_id = "section",
    ion_mode = "positive",
    ion_mode_source = NULL
  ),
  "ion_mode_source must be supplied explicitly"
)

# Accurate-mass annotation preserves one-to-many candidates without confirmation labels.
annotation <- data.frame(
  `Input Mass` = c(500, 500, 700),
  `Matched Mass` = c(500.0001, 500.0002, 710),
  Delta = c(0.0001, 0.0002, 10),
  Name = c("candidate_a", "candidate_b", "outside"),
  Formula = c("C1H1", "C2H2", "C3H3"),
  Ion = c("[M+H]+", "[M+Na]+", "[M+H]+"),
  `LMSD Examples` = c("a", "b", "c"),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
annotation_path <- tempfile(fileext = ".tsv")
write.table(annotation, annotation_path, sep = "\t", row.names = FALSE, quote = FALSE)
feature_metadata <- data.frame(
  feature_id = c("feature_0001", "feature_0002"),
  column_name = c("mz_500", "mz_600"),
  mz = c(500, 600),
  stringsAsFactors = FALSE
)
annotation_mapping <- prepare_lipid_annotation_mapping(annotation_path, feature_metadata, ppm = 5)
stopifnot(nrow(annotation_mapping) == 3L)
stopifnot(sum(annotation_mapping$msi_feature_id == "feature_0001", na.rm = TRUE) == 2L)
stopifnot(identical(annotation_mapping$Name[1:2], c("candidate_a", "candidate_b")))
stopifnot(all(annotation_mapping$annotation_evidence == "putative_accurate_mass_annotation"))
annotation_text <- paste(c(names(annotation_mapping), unlist(annotation_mapping, use.names = FALSE)), collapse = " ")
stopifnot(!grepl("confirmed", annotation_text, ignore.case = TRUE))

# Domain labels must match one-to-one by sample + x + y and remain data-driven labels.
pixels <- data.frame(
  pixel_id = 1:4,
  sample_id = "sample_a",
  section_id = "section_a",
  x = c(1, 2, 1, 2),
  y = c(1, 1, 2, 2),
  mz_500 = 1:4,
  stringsAsFactors = FALSE
)
domains <- data.frame(
  sample = "sample_a",
  x = c(2, 1, 2, 1),
  y = c(2, 1, 1, 2),
  label = c(2, -1, 0, 1),
  stringsAsFactors = FALSE
)
domain_path <- tempfile(fileext = ".csv")
write.csv(domains, domain_path, row.names = FALSE)
mapped_domains <- map_spatial_domain_labels(pixels, domain_path, sample_id = "sample_a")
stopifnot(identical(mapped_domains$domain_id, c("-1", "0", "1", "2")))
stopifnot(mapped_domains$domain_label[1] == "unclassified/background")
stopifnot(all(mapped_domains$domain_type[-1] == "data-driven metabolic domain"))
stopifnot(!any(grepl("anatomical", mapped_domains$domain_type, ignore.case = TRUE)))

duplicate_domains <- rbind(domains, domains[1, , drop = FALSE])
duplicate_domain_path <- tempfile(fileext = ".csv")
write.csv(duplicate_domains, duplicate_domain_path, row.names = FALSE)
expect_error(
  map_spatial_domain_labels(pixels, duplicate_domain_path, sample_id = "sample_a"),
  "Domain labels are not one-to-one"
)

# Local real-data contract test; skipped in source-package checks where data_raw is excluded.
if (!is.null(project_root)) {
  real_root <- file.path(project_root, "data_raw", "msiflow", "ly6g_heterogeneity_signatures")
  real_imzml <- file.path(real_root, "msi", "UPEC_12.imzML")
  if (file.exists(real_imzml)) {
    real <- load_centroided_msi_features(
      real_imzml,
      sample_id = "UPEC_12",
      section_id = "UPEC_12",
      ion_mode = "positive",
      ion_mode_source = "inferred_from_annotation_adducts"
    )
    stopifnot(nrow(real$pixel_feature_matrix) == 7036L)
    stopifnot(nrow(real$feature_metadata) == 252L)
    stopifnot(length(feature_columns(real$pixel_feature_matrix)) == 252L)
    stopifnot(!anyDuplicated(real$feature_metadata$column_name))
    stopifnot(identical(real$feature_metadata$mz, as.numeric(Cardinal::mz(Cardinal::readMSIData(real_imzml)))) )
    stopifnot(real$qc_summary$intensity_source_orientation == "features_x_spectra")
    stopifnot(identical(real$parameters$polarity_confirmed_by_imzml_cv, FALSE))
    stopifnot(real$provenance$value[real$provenance$key == "polarity_cv_metadata_status"] == "not_confirmed")
    real_domains <- map_spatial_domain_labels(
      real$pixel_feature_matrix,
      file.path(real_root, "umap_data.csv"),
      sample_id = "UPEC_12"
    )
    stopifnot(nrow(real_domains) == 7036L)
    stopifnot(!any(is.na(real_domains$domain_id)))
  }
}

cat("IMZML_ADAPTER_TEST_OK=TRUE\n")
