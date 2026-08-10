\name{multimodal_integration}
\alias{find_msi_lcms_candidates}
\alias{assign_msi_lcms_candidates}
\alias{validate_ccs_evidence}
\alias{roi_pixels_to_polygon}
\alias{export_lcm_targets}
\alias{compare_msi_lcm_quantification}
\title{Integrate MSI with LC-MS, CCS, and LCM evidence}
\description{
Functions for auditable multimodal integration: enumerate every LC-MS candidate
within an MSI mass window, optionally select a globally optimal one-to-one
assignment, evaluate explicitly sourced CCS evidence, construct grid-cell ROI
boundaries, export vendor-neutral registered LCM coordinates, and assess
ROI-level MSI--LCM quantitative concordance.
}
\details{
\code{find_msi_lcms_candidates()} preserves the complete candidate set rather
than silently treating a one-to-one assignment as chemical identity.
\code{assign_msi_lcms_candidates()} supports historical greedy selection and a
base-R rectangular Hungarian solver whose objective first maximizes feasible
match cardinality and then minimizes total ppm error.

\code{validate_ccs_evidence()} requires an explicit calibrated CCS tolerance.
An empirical CCS measured only in LC-MS/IMS characterizes LC-MS candidates; it
does not confirm which candidate generated an MSI signal when MSI has no
ion-mobility measurement.

\code{roi_pixels_to_polygon()} treats x and y as regular-grid pixel centers and
returns closed exterior/hole rings at half-pixel boundaries.
\code{export_lcm_targets()} produces vendor-neutral coordinates only. Its output
must be validated against instrument units, axes, fiducials, topology, minimum
area, and cutting constraints before physical use.

\code{compare_msi_lcm_quantification()} computes ROI-level correlation with
optional LCM area/internal-standard normalization and optional within-subject
centering. Correlation is concordance evidence, not compound identification;
within-subject centering is not a mixed-effects model.
}
\seealso{
\code{cross_validate_msi_lcms()},
\code{fit_histology_msi_registration()},
\code{registration_diagnostics()}
}
