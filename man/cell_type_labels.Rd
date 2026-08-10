\name{cell_type_labels}
\alias{prepare_cell_type_labels}
\alias{deconvolve_spatial_transcriptomics}
\alias{segment_cells_from_histology}
\alias{transfer_cell_type_to_msi_pixels}
\alias{compute_neighborhood_composition_soft}
\title{Prepare Multimodal Labels for Niche Analysis}
\description{
Validates external labels, estimates light-weight spatial-transcriptomics
composition scores, performs exploratory watershed nuclear morphology
clustering, transfers labels into an MSI coordinate system, and aggregates
soft cell-type proportions over fixed-size spatial windows.
}
\usage{
prepare_cell_type_labels(labels, proportions = NULL,
  required_fields = c("object_id", "x", "y", "cell_type"),
  entity_type = NULL, label_mode = NULL, coordinate_system = NULL)

deconvolve_spatial_transcriptomics(count_matrix, spot_coordinates,
  reference_signature, method = c("nnls", "cosine_label_transfer"),
  min_confidence = 0.1, min_shared_genes = 10,
  region_column = NULL,
  coordinate_system = "spatial_transcriptomics")

segment_cells_from_histology(image_path, tissue_mask = NULL,
  min_nucleus_size = 20, max_nucleus_size = 500,
  n_morphology_classes = 3, opening_size = 3,
  watershed_tolerance = 1, watershed_ext = 1, seed = 1,
  region_id = "region_1", coordinate_system = "histology_pixel")

transfer_cell_type_to_msi_pixels(pixel_matrix, cell_type_labels,
  method = c("shared_coordinates", "nearest_neighbor"),
  max_distance = NULL, pixel_region_column = "section_id",
  registration = NULL)

compute_neighborhood_composition_soft(proportions, x, y, k = 10,
  pixel_id = seq_len(nrow(proportions)),
  region_id = rep("region_1", nrow(proportions)), include_self = TRUE)
}
\arguments{
  \item{labels}{Data frame containing spatial object labels.}
  \item{proportions}{Optional position-by-cell-type proportion matrix.}
  \item{required_fields}{Required label-table fields.}
  \item{entity_type, label_mode, coordinate_system}{Explicit provenance fields.}
  \item{count_matrix}{Non-negative, linear-scale spot-by-gene matrix.}
  \item{spot_coordinates}{Spot identifiers and coordinates.}
  \item{reference_signature}{Non-negative gene-by-cell-type reference matrix.}
  \item{method}{Deconvolution/transfer method.}
  \item{min_confidence}{Dominant-score threshold below which a label is unresolved.}
  \item{min_shared_genes}{Minimum shared target/reference genes.}
  \item{region_column, pixel_region_column, region_id}{Section identifiers used to prohibit cross-section matching.}
  \item{image_path}{Path to one two-dimensional histology image.}
  \item{tissue_mask}{Optional logical mask with image spatial dimensions.}
  \item{min_nucleus_size, max_nucleus_size}{Accepted nuclear areas in image pixels.}
  \item{n_morphology_classes}{Number of exploratory morphology clusters.}
  \item{opening_size}{Disc-brush size for morphological opening.}
  \item{watershed_tolerance, watershed_ext}{Arguments passed to \code{EBImage::watershed()}.}
  \item{seed}{Optional k-means seed.}
  \item{pixel_matrix}{MSI pixel table.}
  \item{cell_type_labels}{A validated label object or compatible data frame.}
  \item{max_distance}{Required finite matching radius for nearest-neighbor transfer.}
  \item{registration}{Optional result from \code{fit_histology_msi_registration()}.}
  \item{x, y, k, pixel_id, include_self}{Spatial-window inputs.}
}
\details{
NNLS is a transparent baseline, not a replacement for a validated probabilistic
deconvolution workflow. Cosine outputs are normalized similarity scores and
must not be described as cell fractions. For multi-cell spots, preserve the
soft NNLS proportions and use \code{compute_neighborhood_composition_soft()}.

The histology path applies global luminance Otsu thresholding, opening, hole
filling, a distance transform, watershed segmentation, and k-means on area,
circularity, mean intensity, and intensity variation. Its labels are nuclear
morphology clusters rather than biological cell types. Their confidence is
relative separation from the two closest k-means centers, not classification
probability. Histology coordinates use pixels with an upper-left origin and
normally require registration before MSI transfer.
}
\value{
Label preparation and source functions return a \code{cell_type_labels} list.
Transfer returns the augmented pixel matrix, transferred proportions, matching
table, QC, and settings. Soft composition returns a
\code{neighborhood_composition} object accepted by \code{define_niches()}.
}
\references{
Pau G, Fuchs F, Sklyar O, Boutros M, Huber W. EBImage: an R package for image
processing with applications to cellular phenotypes. Bioinformatics 26,
979--981 (2010).
}
