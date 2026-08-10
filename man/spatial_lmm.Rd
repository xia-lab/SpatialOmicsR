\name{spatial_lmm}
\alias{differential_region_analysis_spatial_lmm}
\alias{differential_region_analysis_glmm}
\title{Gaussian Spatial Mixed-Effects Region Analysis}
\description{
Fits one Gaussian linear mixed-effects model per MSI feature with nested
subject/section random intercepts and an exponential, Gaussian, or spherical
within-field spatial residual correlation. Provides marginal omnibus group
tests, fixed-effect pairwise contrasts, covariance parameters, model-selection
summaries, and auditable spatial-field QC.
}
\usage{
differential_region_analysis_spatial_lmm(sample_matrix,
  group_column = "roi_id", subject_column, section_column,
  x_col = NULL, y_col = NULL,
  x_bounds = c("x_min", "x_max"),
  y_bounds = c("y_min", "y_max"),
  x_resolution, y_resolution, distance_unit,
  correlation_structure = c("exponential", "gaussian", "spherical"),
  reference_group = NULL, p_adjust_method = "BH",
  min_subjects_per_group = 5,
  min_unique_coordinates_per_field = 4,
  field_action = c("error", "drop"), confidence_level = 0.95,
  keep_fits = FALSE)

differential_region_analysis_glmm(...)
}
\arguments{
  \item{sample_matrix}{Subregion/sample-level table with \code{mz_} features.}
  \item{group_column}{Region/group fixed-effect column.}
  \item{subject_column}{Biological-subject identifier; required.}
  \item{section_column}{Section/spatial-field identifier nested within subject; required.}
  \item{x_col, y_col}{Optional subregion-center coordinates.}
  \item{x_bounds, y_bounds}{Minimum/maximum columns used to calculate centers when explicit center columns are absent.}
  \item{x_resolution, y_resolution}{Physical units represented by one coordinate unit.}
  \item{distance_unit}{Explicit coordinate unit saved in results.}
  \item{correlation_structure}{Stationary isotropic spatial covariance model.}
  \item{reference_group}{Optional first factor/reference group.}
  \item{p_adjust_method}{Adjustment method for feature tests and pairwise contrasts.}
  \item{min_subjects_per_group}{Required biological subjects represented in every group.}
  \item{min_unique_coordinates_per_field}{Required unique centers in every retained subject-section field.}
  \item{field_action}{Stop on an invalid spatial field or explicitly drop it.}
  \item{confidence_level}{Pairwise contrast confidence level.}
  \item{keep_fits}{Retain potentially large fitted model objects.}
  \item{...}{Arguments forwarded by the compatibility alias.}
}
\details{
The fitted model is Gaussian, not a generalized linear mixed model. It is
appropriate only when the modeled feature scale has a defensible approximately
Gaussian conditional response. Binary, count, proportion, zero-inflated, and
hurdle responses require separately validated non-Gaussian models.

Random intercepts are nested as subject/section. Residual spatial correlation
is restricted to the same subject-section using
\code{~ x + y | subject/section}; different sections are never connected.
Repeated coordinates within a spatial field are rejected rather than jittered.

The spatial model is the inferential model. The independence model and its
AIC/BIC are sensitivity summaries using the same fixed and random effects.
No ordinary chi-square likelihood-ratio p-value is reported for adding the
spatial covariance because the independence null is a covariance-boundary case.

Range and nugget are extracted on their natural constrained scale. The range
parameter has covariance-function-specific meaning; it must not generically be
called the distance beyond which correlation vanishes. Pairwise estimates are
fixed-effect contrasts in the direction group B minus group A.
}
\value{
A list containing feature-level omnibus tests and covariance summaries,
pairwise contrasts, field QC, subject counts, optional fits, settings, and an
interpretation statement. The legacy \code{differential_region_analysis_glmm()}
name emits a warning and forwards to the Gaussian spatial LMM.
}
\references{
Pinheiro JC, Bates DM. Mixed-Effects Models in S and S-PLUS. Springer (2000).
}
