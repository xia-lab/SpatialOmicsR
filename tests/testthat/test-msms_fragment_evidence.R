testthat::test_that("regression: test-msms_fragment_evidence", {
checking_installed_package <- nzchar(Sys.getenv("_R_CHECK_PACKAGE_NAME_"))
if (checking_installed_package) {
  library(SpatialOmicsMSI)
} else {
  source(file.path(normalizePath(getwd()), "R", "msi_pipeline.R"))
}

scans <- data.frame(
  scan_number = c(1, 2, 3, 4),
  precursor_mz = c(500.0005, 500.0004, 499.9995, 500.02),
  isolation_target_mz = c(500, 502, 500, 500),
  isolation_lower_offset = rep(0.5, 4),
  isolation_upper_offset = rep(0.5, 4),
  base_peak_intensity = rep(1000, 4),
  stringsAsFactors = FALSE
)
peaks <- data.frame(
  scan_number = c(1, 3, 2, 1, 1, 3, 1, 3),
  fragment_mz = c(100.0005, 99.9997, 100.0001, 200.001, 300.001, 300.001, 400.2, 400.2),
  fragment_intensity = c(100, 120, 500, 50, 5, 9, 100, 100),
  peak_source = "profile_derived",
  extraction_method = "local maxima",
  extraction_parameters = "prominence=0",
  software_version = "test 1.0",
  stringsAsFactors = FALSE
)
diagnostics <- data.frame(
  label = c("repeat", "single", "trace", "absent"),
  target_mz = c(100, 200, 300, 400),
  evidence_role = c("role_a", "role_b", "role_c", "role_d"),
  stringsAsFactors = FALSE
)

result <- validate_msms_fragment_evidence(
  scans, peaks, diagnostics,
  precursor_target = 500,
  precursor_ppm_tolerance = 10,
  fragment_tolerance = 0.02,
  fragment_tolerance_unit = "Da",
  min_relative_intensity = 0.01,
  min_core_scans = 2,
  required_diagnostic_labels = c("repeat", "absent")
)

evidence <- result$scan_diagnostic_evidence
summary <- setNames(result$diagnostic_summary$evidence_grade,
                    result$diagnostic_summary$diagnostic_label)
stopifnot(!evidence$core_scan[evidence$scan_number == 2][1])
stopifnot(identical(unname(summary["single"]), "single_scan_support"))
stopifnot(identical(unname(summary["repeat"]), "repeated_diagnostic_support"))
stopifnot(identical(unname(summary["trace"]), "trace_match"))
stopifnot(identical(unname(summary["absent"]), "not_detected"))
stopifnot(!result$overall$support)
stopifnot(all(evidence$peak_source[evidence$mass_match] == "profile_derived"))
stopifnot(all(evidence$extraction_method[evidence$mass_match] == "local maxima"))

# A mass match in a non-core isolation scan cannot create qualifying evidence.
repeat_scan_2 <- evidence[evidence$scan_number == 2 & evidence$diagnostic_label == "repeat", ]
stopifnot(repeat_scan_2$mass_match, !repeat_scan_2$core_scan,
          !repeat_scan_2$qualifies_core_strength_rule)

# PPM tolerance is supported independently of Da tolerance.
ppm_diagnostic <- data.frame(label = "ppm_test", target_mz = 100,
                             evidence_role = "tolerance_test")
ppm_pass <- validate_msms_fragment_evidence(
  scans[1, ], peaks[1, ], ppm_diagnostic,
  precursor_target = 500, precursor_ppm_tolerance = 10,
  fragment_tolerance = 6, fragment_tolerance_unit = "ppm",
  min_relative_intensity = 0.01, min_core_scans = 2,
  required_diagnostic_labels = "ppm_test"
)
ppm_fail <- validate_msms_fragment_evidence(
  scans[1, ], peaks[1, ], ppm_diagnostic,
  precursor_target = 500, precursor_ppm_tolerance = 10,
  fragment_tolerance = 4, fragment_tolerance_unit = "ppm",
  min_relative_intensity = 0.01, min_core_scans = 2,
  required_diagnostic_labels = "ppm_test"
)
stopifnot(ppm_pass$scan_diagnostic_evidence$mass_match)
stopifnot(!ppm_fail$scan_diagnostic_evidence$mass_match)

# Input order cannot change any returned result.
reordered <- validate_msms_fragment_evidence(
  scans[4:1, ], peaks[nrow(peaks):1, ], diagnostics[4:1, ],
  precursor_target = 500, precursor_ppm_tolerance = 10,
  fragment_tolerance = 0.02, fragment_tolerance_unit = "Da",
  min_relative_intensity = 0.01, min_core_scans = 2,
  required_diagnostic_labels = c("repeat", "absent")
)
stopifnot(identical(result, reordered))

# No output field or value may use confirmation terminology.
serialized <- tolower(paste(capture.output(dput(result)), collapse = " "))
stopifnot(!grepl("confirmed", serialized, fixed = TRUE))

cat("MSMS_FRAGMENT_EVIDENCE_TEST_OK=TRUE\n")

  testthat::succeed()
})
