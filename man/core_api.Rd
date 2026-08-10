\name{SpatialOmicsMSI-core-api}
\alias{apply_matched_region_labels}
\alias{apply_roi_labels}
\alias{backmap_sample_scores}
\alias{build_spatial_neighbors}
\alias{cluster_diagnostics}
\alias{cluster_histology_superpixels}
\alias{cluster_pixels}
\alias{compare_contamination_methods}
\alias{compute_morans_i_grid}
\alias{compute_spatially_variable_metabolites}
\alias{cross_validate_annotations}
\alias{cross_validate_msi_lcms}
\alias{define_rois}
\alias{detect_metaboanalyst_result}
\alias{detect_ubiquitous_features}
\alias{differential_region_analysis}
\alias{feature_columns}
\alias{feature_stats}
\alias{filter_background_cardinal_spatial_kmeans}
\alias{fit_histology_msi_registration}
\alias{generate_roi_candidates}
\alias{import_peakpicked_msi_csv}
\alias{infer_feature_mapping}
\alias{link_metaboanalyst_pathway_features}
\alias{load_centroided_msi_features}
\alias{load_msi_target_features}
\alias{load_peakpicked_msi_series}
\alias{load_serial_msi_target_features}
\alias{make_feature_mapping}
\alias{make_metaboanalyst_data}
\alias{make_pipeline_manifest}
\alias{map_spatial_domain_labels}
\alias{nearest_intensity_ppm}
\alias{normalize_metaboanalyst_result}
\alias{optimize_histology_rois}
\alias{otsu_threshold}
\alias{plot_feature_region_boxplot}
\alias{plot_feature_thumbnails}
\alias{plot_ion_image}
\alias{plot_metabo_feature_rank}
\alias{plot_metabo_feature_view}
\alias{plot_metabo_score_plot}
\alias{plot_pathway_feature_images}
\alias{plot_pathway_score_map}
\alias{plot_sample_score_map}
\alias{prepare_lipid_annotation_mapping}
\alias{preprocess_matrix}
\alias{preprocess_sections}
\alias{preprocess_select_features}
\alias{reduce_pixel_matrix}
\alias{registration_diagnostics}
\alias{run_spatial_app}
\alias{run_spatial_metabolomics_pipeline}
\alias{sample_coordinate_rois}
\alias{sample_matched_regions}
\alias{sample_subregions}
\alias{score_roi_indices}
\alias{select_features}
\alias{select_manual_rois}
\alias{select_rois}
\alias{select_shared_features}
\alias{spatial_neighbor_diagnostics}
\alias{transform_histology_coordinates}
\alias{write_pipeline_outputs}
\alias{write_validation_plots}
\title{Core SpatialOmicsMSI data-processing and analysis API}
\description{
Public functions for MSI import, quality control, preprocessing, feature
selection, spatial graphs and autocorrelation, clustering, ROI construction,
registration, differential analysis, annotation integration, pathway
back-mapping, visualization, provenance, and export.
}
\details{
The package is organized as a function-first workflow. Import functions create
pixel, coordinate, feature-metadata, QC, parameter, and provenance tables.
Preprocessing functions provide explicit background filtering, TIC
normalization, log transformation, and feature selection. Spatial functions
construct reusable neighbor graphs for autocorrelation and clustering. ROI and
registration functions preserve coordinate and section provenance. Statistical
functions distinguish exploratory pixel-level summaries from analyses that use
biological subjects as the inference unit.

Individual function signatures are discoverable with \code{args()} and are
demonstrated in the package README and workflow vignette. Specialized methods
with additional scientific assumptions have dedicated help pages.
}
\seealso{
\code{help(package = "SpatialOmicsMSI")}
}
