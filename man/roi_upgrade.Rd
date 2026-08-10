\name{roi_upgrade}
\alias{corroborate_cluster_labels}
\alias{generate_roi_candidates_from_domains}
\alias{point_in_polygon_with_holes}
\alias{import_qupath_geojson}
\alias{validate_lmd_shapes}
\title{Construct and validate multimodal and externally annotated ROIs}
\description{
Functions supporting cross-modal label agreement, connected-domain automatic
ROI candidates, polygons with holes, QuPath GeoJSON annotations, and read-only
physical-geometry QC for laser microdissection.
}
\details{
\code{corroborate_cluster_labels()} retains only pixels whose actual pair of
labels belongs to an accepted one-way or mutual-best overlap. Its result is a
cross-modal agreement summary, not independent validation or biological
identity assignment.

\code{generate_roi_candidates_from_domains()} constructs one candidate per
connected component of a classification column. Rook connectivity is the
default; topology steps must describe the coordinate grid. The
\code{cluster_column} used by this and the other ROI functions may contain
histology clusters, chemical spatial domains, or explicitly corroborated
cross-modal labels. The provenance must be reported with downstream results.

Manual polygon ROI tables use the hierarchy
\code{roi_id -> polygon_part_id -> ring_id}. Every part must have exactly one
\code{outer} ring and may have zero or more \code{hole} rings. Multiple parts
with the same ROI identity are combined only after holes have been subtracted
from their owning part.

\code{import_qupath_geojson()} accepts Polygon and MultiPolygon Feature data.
Coordinates are recorded as full-resolution image pixels with a top-left
origin and downward y-axis, following QuPath export semantics. Registration is
required before transferring them to MSI or physical stage coordinates.

\code{validate_lmd_shapes()} reports area, spacing, ring-role and basic
self-intersection failures. It deliberately does not split, buffer, simplify,
or export instrument commands. Manufacturer constraints and calibrated stage
coordinates must be independently verified before cutting.
}
\seealso{
\code{select_rois()}, \code{optimize_histology_rois()},
\code{fit_histology_msi_registration()}, \code{export_lcm_targets()}
}
