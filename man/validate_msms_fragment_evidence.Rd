\name{validate_msms_fragment_evidence}
\alias{validate_msms_fragment_evidence}
\title{Validate Explicit Diagnostic-Ion Evidence in MS/MS Scans}
\usage{
validate_msms_fragment_evidence(
  precursor_scans,
  fragment_peaks,
  diagnostic_ions,
  precursor_target,
  precursor_ppm_tolerance = 10,
  fragment_tolerance = 0.02,
  fragment_tolerance_unit = c("Da", "ppm"),
  min_relative_intensity = 0.01,
  min_core_scans = 2L,
  required_diagnostic_labels
)
}
\arguments{
  \item{precursor_scans}{Data frame with scan number, selected precursor m/z,
  isolation-window target and offsets, and base-peak intensity.}
  \item{fragment_peaks}{Data frame with scan number, fragment m/z, fragment
  intensity, and optional peak-extraction provenance.}
  \item{diagnostic_ions}{Data frame containing unique \code{label}, numeric
  \code{target_mz}, and explicit \code{evidence_role} columns.}
  \item{precursor_target}{Positive numeric precursor target m/z.}
  \item{precursor_ppm_tolerance}{Non-negative precursor tolerance in ppm.}
  \item{fragment_tolerance}{Non-negative fragment matching tolerance.}
  \item{fragment_tolerance_unit}{Either \code{"Da"} or \code{"ppm"}.}
  \item{min_relative_intensity}{Minimum fragment intensity divided by scan
  base-peak intensity.}
  \item{min_core_scans}{Minimum number of qualifying core scans required for
  repeated diagnostic support; must be at least two.}
  \item{required_diagnostic_labels}{Diagnostic labels that must all achieve
  repeated diagnostic support for the overall rule to pass.}
}
\value{
A list with \code{scan_diagnostic_evidence}, \code{diagnostic_summary}, and
\code{overall}. The recommended overall label is
\code{fragment_level_orthogonal_support}.
}
\description{
Evaluates user-defined diagnostic ions without inferring chemical identity. A
core scan must have a selected precursor within the requested ppm tolerance and
an isolation window that covers the target precursor. Fragment evidence is
graded from mass matching, relative intensity, and the number of qualifying
core scans.
}
