orig_wd <- getwd()
checking_installed_package <- nzchar(Sys.getenv("_R_CHECK_PACKAGE_NAME_"))
if (checking_installed_package) {
  library(SpatialOmicsMSI)
} else {
  source(file.path(orig_wd, "R", "msi_pipeline.R"))
  source(file.path(orig_wd, "R", "spatial_registration.R"))
  source(file.path(orig_wd, "R", "multimodal_integration.R"))
}

# Greedy consumes the only candidate available to b; optimal assignment first
# maximizes cardinality and therefore retains both feasible matches.
msi <- data.frame(
  feature_id = c("a", "b"), mz = c(100, 99.9992), ion_mode = "positive",
  stringsAsFactors = FALSE
)
lcms <- data.frame(
  id = c("x", "y"), mz = c(100.0001, 100.0008), ion_mode = "positive",
  stringsAsFactors = FALSE
)
candidates <- find_msi_lcms_candidates(msi, lcms, ppm = 10)
stopifnot(nrow(candidates) == 3L, candidates$candidate_count[candidates$msi_row == 1L] == 2L)
greedy <- cross_validate_msi_lcms(msi, lcms, ppm = 10, assignment_method = "greedy")
optimal <- cross_validate_msi_lcms(msi, lcms, ppm = 10, assignment_method = "optimal")
stopifnot(nrow(greedy) == 1L, nrow(optimal) == 2L)
stopifnot(identical(sort(optimal$msi_feature_id), c("a", "b")))
stopifnot(identical(unique(optimal$assignment_method), "optimal"))

# Without an exact-cost tie, input row order does not alter the selected stable
# feature-ID pairs.
reordered <- cross_validate_msi_lcms(msi[2:1, ], lcms[2:1, ], ppm = 10,
                                     assignment_method = "optimal")
pair_key <- function(x) sort(paste(x$msi_feature_id, x$lcms_feature_id, sep = ":"))
stopifnot(identical(pair_key(optimal), pair_key(reordered)))

observed <- data.frame(
  observation_id = c("obs1", "obs1"), candidate_id = c("lipid_a", "lipid_b"),
  observed_ccs = c(200, 200), observed_source = "lcms_empirical",
  adduct = "[M+H]+", charge = 1L, ion_mode = "positive",
  stringsAsFactors = FALSE
)
reference <- data.frame(
  candidate_id = c("lipid_a", "lipid_b"), reference_ccs = c(201, 202),
  reference_source = c("measured_library", "ml_predicted"),
  adduct = "[M+H]+", charge = 1L, ion_mode = "positive",
  stringsAsFactors = FALSE
)
ccs <- validate_ccs_evidence(
  observed, reference, observation_id_column = "observation_id",
  ccs_tolerance_pct = 2
)
stopifnot(all(ccs$within_tolerance), all(ccs$ambiguous_candidate_set))
stopifnot(!any(ccs$msi_identity_confirmed_by_ccs))
stopifnot(all(ccs$evidence_scope == "lcms_candidate_characterization_only"))

observed_msi <- observed[1L, ]
observed_msi$observed_source <- "msi_empirical"
ccs_msi <- validate_ccs_evidence(observed_msi, reference[1L, ], ccs_tolerance_pct = 2)
stopifnot(ccs_msi$msi_identity_confirmed_by_ccs)

polygon <- roi_pixels_to_polygon(c(0, 1, 0, 1), c(0, 0, 1, 1), 1)
stopifnot(nrow(polygon) == 9L, polygon$x[1L] == polygon$x[nrow(polygon)])
stopifnot(polygon$y[1L] == polygon$y[nrow(polygon)], all(polygon$ring_type == "exterior"))
separate <- roi_pixels_to_polygon(c(0, 2), c(0, 0), 1)
stopifnot(length(unique(separate$ring_id)) == 2L)

control <- data.frame(
  histology_x = c(0, 10, 0, 10), histology_y = c(0, 0, 10, 10),
  msi_x = c(100, 120, 100, 120), msi_y = c(200, 200, 230, 230)
)
registration <- fit_histology_msi_registration(control)
roi_pixels <- data.frame(roi_id = "r1", x = c(0, 1, 0, 1), y = c(0, 0, 1, 1))
lcm_target <- export_lcm_targets(roi_pixels, registration, 1)
stopifnot(all(is.finite(lcm_target$stage_x)), all(is.finite(lcm_target$stage_y)))
stopifnot(all(c("registration_rmse", "registration_max_error") %in% names(lcm_target)))

msi_roi <- data.frame(
  roi_id = paste0("r", 1:5), section_id = "s1", subject = c("a", "a", "b", "b", "c"),
  mz_100 = 1:5
)
lcm_roi <- data.frame(
  roi_id = paste0("r", 1:5), section_id = "s1", compound_a = 2 * (1:5),
  area = 1, internal_standard = 2
)
mapping <- data.frame(msi_feature = "mz_100", lcm_feature = "compound_a")
comparison <- compare_msi_lcm_quantification(
  msi_roi, lcm_roi, mapping, lcm_area_column = "area",
  lcm_internal_standard_column = "internal_standard"
)
stopifnot(comparison$n_pairs == 5L, abs(comparison$correlation - 1) < 1e-12)
stopifnot(comparison$lcm_area_normalized, comparison$lcm_internal_standard_normalized)

cat("MULTIMODAL_INTEGRATION_TEST_OK=TRUE\n")
