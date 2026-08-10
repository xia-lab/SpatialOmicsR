\name{real_data_adapters}
\alias{validate_imzml_ibd_pair}
\alias{load_variable_mz_msi_features}
\alias{read_metaspace_transform}
\alias{apply_metaspace_transform}
\alias{register_metaspace_optical}
\alias{read_mzml_fragment_spectra}
\alias{build_msi_tissue_mask}
\title{Adapters for variable-axis MSI, METASPACE optical registration, and mzML MS2}
\description{
Read-only adapters used to standardize processed variable-axis imzML data,
apply an externally supplied METASPACE transform, and extract MS2 precursor and
fragment arrays from open mzML files.
}
\usage{
validate_imzml_ibd_pair(imzml_path, ibd_path = NULL)

load_variable_mz_msi_features(imzml_path, ibd_path = NULL,
  sample_id, section_id, ion_mode, ion_mode_source,
  processing = c("processed_peak_lists", "profile_diff"),
  alignment_ppm = 10, peak_pick_snr = 6,
  min_detection_fraction = 0, progress_callback = NULL)

read_metaspace_transform(path)

apply_metaspace_transform(coordinates, transform,
  origin = c(min(coordinates$x), min(coordinates$y)))

register_metaspace_optical(coordinates, transform_json, optical_path,
  attribution_path = NULL, tic = NULL, representative_ions = NULL,
  output_dir = NULL, tissue_threshold = 165 / 255)

read_mzml_fragment_spectra(mzml_path, precursor_target = NULL,
  precursor_ppm_tolerance = 10)

build_msi_tissue_mask(coordinates, tic, peak_count,
  method = c("kmeans_log_tic_peak_count", "score_quantile"),
  score_quantile = 0.50, seed = 20260808)
}
\arguments{
  \item{imzml_path, ibd_path}{Paths to paired imzML and ibd files.}
  \item{sample_id, section_id}{Explicit sample and section identifiers.}
  \item{ion_mode, ion_mode_source}{Explicit polarity and provenance for that claim.}
  \item{processing}{Either alignment of stored processed peak lists or profile peak picking followed by alignment.}
  \item{alignment_ppm}{Positive alignment tolerance in ppm.}
  \item{peak_pick_snr}{Signal-to-noise threshold for the profile-difference peak picker.}
  \item{min_detection_fraction}{Minimum fraction of pixels detecting an aligned feature.}
  \item{progress_callback}{Optional function accepting progress and detail.}
  \item{path, transform_json}{Path to the unmodified platform-transform JSON.}
  \item{coordinates}{Data frame containing x and y.}
  \item{transform}{Finite 3 by 3 homogeneous transformation matrix.}
  \item{origin}{MSI coordinate origin; by default inferred from coordinate minima.}
  \item{optical_path}{Path to the optical-brightfield image.}
  \item{attribution_path}{Optional attribution/license JSON.}
  \item{tic}{Per-pixel TIC values; optional for registration overlays and required for tissue masking.}
  \item{representative_ions}{Optional named list of pixel intensity vectors.}
  \item{output_dir}{Optional directory for registration overlay images.}
  \item{tissue_threshold}{Optical grayscale threshold used only for overlap diagnostics.}
  \item{mzml_path}{Path to an open mzML file.}
  \item{precursor_target}{Optional selected-precursor target m/z.}
  \item{precursor_ppm_tolerance}{Positive precursor window in ppm.}
  \item{peak_count}{Raw per-pixel detected-array-value counts.}
  \item{method}{Transparent tissue-candidate rule.}
  \item{score_quantile}{Quantile used by the sensitivity threshold rule.}
  \item{seed}{Fixed random seed for the two-class k-means rule.}
}
\value{
Lists containing validated paths or standardized matrices, coordinates,
metadata, parameters and provenance. The mzML reader returns precursor scan
metadata and an unassigned fragment peak table; it does not infer identities.
}
\details{
Spectrum representation is reported from CV metadata and is never inferred from
the filename. The METASPACE transform direction is relative MSI coordinates to
optical pixels. No chemical identity or confirmation label is inferred from an
mzML fragment array.
}
