\name{distance_gradient}
\alias{compute_reference_distance}
\alias{bin_distance_to_rings}
\alias{aggregate_distance_profiles}
\alias{fit_distance_gam}
\alias{compute_distance_variable_metabolites}
\alias{analyze_distance_gradient}
\title{Auditable Distance-to-Domain Gradient Analysis}
\description{
Computes section-aware physical distances to domain interfaces, observed-grid
perimeters, or component centroids; constructs signed distance rings; aggregates
subject profiles; and fits descriptive pixel-level or population-level GAMs.
}
\usage{
compute_reference_distance(pixel_matrix, domain_column, target_domain,
  x_resolution, y_resolution, topology_x_step, topology_y_step,
  distance_unit, reference = c("boundary", "centroid"),
  boundary_type = c("domain_interface", "observed_domain_perimeter",
    "domain_interface_and_perimeter"),
  x_col = "x", y_col = "y", section_column = NULL,
  neighbor_method = c("queen", "rook"),
  component_action = c("error", "nearest", "largest"))

bin_distance_to_rings(distance_result, method, n_rings = 5,
  ring_width = NULL, breaks = NULL, separate_sides = TRUE)

aggregate_distance_profiles(pixel_matrix, distance_result, rings,
  subject_column, section_column = NULL,
  features = feature_columns(pixel_matrix),
  summary_method = c("mean", "median"),
  min_pixels_per_subject_ring = 10, min_subjects_per_ring = 2)

fit_distance_gam(intensity, distance,
  analysis_mode = c("exploratory_pixel", "population_subject"),
  subject_id = NULL, section_id = NULL, data_are_aggregated = FALSE,
  k = 10, subject_k = NULL, min_subjects = 5,
  family = stats::gaussian())

compute_distance_variable_metabolites(pixel_matrix, distance_result,
  features = feature_columns(pixel_matrix),
  analysis_mode = c("exploratory_pixel", "population_subject"),
  aggregated_profiles = NULL, k = 10, subject_k = NULL,
  min_subjects = 5, p_adjust_method = "BH")

analyze_distance_gradient(pixel_matrix, domain_column, target_domain,
  ring_method, x_resolution, y_resolution, topology_x_step,
  topology_y_step, distance_unit, domain_source_features = NULL,
  test_features = NULL, reference = c("boundary", "centroid"),
  boundary_type = c("domain_interface", "observed_domain_perimeter",
    "domain_interface_and_perimeter"),
  x_col = "x", y_col = "y", section_column = NULL,
  neighbor_method = c("queen", "rook"),
  component_action = c("error", "nearest", "largest"),
  n_rings = 5, ring_width = NULL, ring_breaks = NULL,
  separate_sides = TRUE,
  analysis_mode = c("exploratory_pixel", "population_subject"),
  subject_column = NULL, gam_k = 10, subject_k = NULL,
  min_subjects = 5, profile_summary = c("mean", "median"),
  min_pixels_per_subject_ring = 10, min_subjects_per_ring = 2,
  allow_circular_analysis = FALSE, p_adjust_method = "BH")
}
\arguments{
  \item{pixel_matrix}{Pixel table containing coordinates, domains, and numeric features.}
  \item{domain_column, target_domain}{Domain-label column and target value.}
  \item{x_resolution, y_resolution}{Physical units represented by one coordinate unit. These have no defaults.}
  \item{topology_x_step, topology_y_step}{Coordinate increments representing adjacent grid positions. Topology is built after dividing by these steps; physical distances are computed separately.}
  \item{distance_unit}{Explicit output unit, such as \code{"um"} or \code{"pixel"}.}
  \item{reference}{Distance to a boundary or nearest retained-component centroid.}
  \item{boundary_type}{Target/non-target interface, observed acquisition-grid perimeter, or their union.}
  \item{x_col, y_col, section_column}{Coordinate and optional section columns.}
  \item{neighbor_method}{Rook or Queen grid topology.}
  \item{component_action}{Stop for disconnected targets, retain all and use the nearest reference, or retain only the largest component.}
  \item{distance_result}{Object returned by \code{compute_reference_distance()}.}
  \item{method, ring_method}{Required ring strategy: quantile, fixed width, or fixed breaks.}
  \item{n_rings, ring_width, breaks, ring_breaks}{Ring construction parameters.}
  \item{separate_sides}{For boundary distances, construct inner and outer rings separately.}
  \item{rings}{Object returned by \code{bin_distance_to_rings()}.}
  \item{subject_column, subject_id}{Biological-subject identifiers.}
  \item{features, test_features}{Numeric features to analyze.}
  \item{summary_method, profile_summary}{Subject-profile aggregation statistic.}
  \item{min_pixels_per_subject_ring, min_subjects_per_ring}{Ring eligibility thresholds.}
  \item{intensity, distance}{Response and distance vectors.}
  \item{analysis_mode}{Descriptive pixels or inferential aggregated subject profiles.}
  \item{section_id}{Optional section identifiers for aggregated profiles.}
  \item{data_are_aggregated}{Required TRUE safeguard for population mode.}
  \item{k, gam_k}{Population smooth basis dimension.}
  \item{subject_k}{Optional factor-smooth basis dimension for subject-specific curve deviations.}
  \item{min_subjects}{Required number of biological subjects in population mode.}
  \item{family}{Response family passed to \code{mgcv::gam()}.}
  \item{aggregated_profiles}{Subject-by-section-by-ring profiles.}
  \item{p_adjust_method}{Multiple-testing adjustment method.}
  \item{domain_source_features}{Features used to construct the domain; use \code{character(0)} only for an independently defined domain.}
  \item{allow_circular_analysis}{Permit explicitly descriptive reuse of domain-defining features.}
}
\details{
Topology and metric coordinates are deliberately separate. Rook/Queen graphs
are built on grid coordinates normalized by the declared topology steps,
whereas distances use coordinates multiplied by their declared physical
resolutions. Sections never share graph edges or reference points.

Boundary distance is signed: retained target pixels are non-positive and other
observed pixels non-negative. Zero denotes a selected boundary pixel.
\code{observed_domain_perimeter} means incomplete expected degree on the
observed grid. It includes acquisition edges, holes, and missing grid positions;
it must not automatically be called a biological tissue boundary.

Population GAMs reject raw pixels and require aggregated profiles. The model
contains a population distance smooth, a subject random intercept, and an
optional nested subject-section random intercept and factor-smooth deviation.
Smooth-term p-values are approximate and require residual, basis-dimension, and
subject-coverage diagnostics. Exploratory pixel GAMs return no p-value.

If the analyzed feature helped define the target domain, its distance gradient
is circular evidence. The integrated function rejects overlap unless explicitly
permitted and records the declaration in provenance.
}
\value{
The low-level functions return structured objects with matrices, settings, QC,
breaks, fitted models, and failure messages. \code{analyze_distance_gradient()}
returns distances, rings, aggregation/QC, continuous and discrete results,
annotated pixels, provenance, and interpretation.
}
\references{
Wood SN. Generalized Additive Models: An Introduction with R. Second edition.
CRC Press (2017).
}
