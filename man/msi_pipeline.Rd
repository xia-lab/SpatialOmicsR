\name{msi_pipeline}
\alias{msi_pipeline}
\alias{apply_matched_region_labels}
\alias{feature_columns}
\alias{infer_feature_mapping}
\alias{import_peakpicked_msi_csv}
\alias{load_peakpicked_msi_series}
\alias{load_msi_target_features}
\alias{load_serial_msi_target_features}
\alias{make_feature_mapping}
\alias{nearest_intensity_ppm}
\alias{feature_stats}
\alias{select_features}
\alias{select_shared_features}
\alias{reduce_pixel_matrix}
\alias{preprocess_select_features}
\alias{preprocess_matrix}
\alias{cluster_pixels}
\alias{cluster_diagnostics}
\alias{sample_subregions}
\alias{sample_matched_regions}
\alias{make_metaboanalyst_data}
\alias{detect_metaboanalyst_result}
\alias{normalize_metaboanalyst_result}
\alias{backmap_sample_scores}
\alias{plot_ion_image}
\alias{plot_feature_thumbnails}
\alias{plot_feature_region_boxplot}
\alias{plot_metabo_feature_rank}
\alias{plot_metabo_feature_view}
\alias{plot_metabo_score_plot}
\alias{plot_sample_score_map}
\alias{preprocess_sections}
\alias{run_spatial_app}
\alias{write_pipeline_outputs}
\alias{write_validation_plots}
\title{Spatial MSI analysis pipeline}
\description{
Function-first tools for targeted MSI feature extraction, feature selection,
preprocessing, spatial clustering, sub-region sampling, MetaboAnalyst export,
MetaboAnalyst result import normalization, and spatial back-mapping.
}
\usage{
feature_columns(data)
apply_matched_region_labels(pixel_matrix, region_labels,
  region_column = "matched_region_label")
infer_feature_mapping(pixel_matrix)
import_peakpicked_msi_csv(path, section_id = NULL, section_order = 1,
  sample_id = NULL)
load_peakpicked_msi_series(paths, section_ids = NULL,
  section_order = NULL, sample_ids = NULL)
load_msi_target_features(imzml_path, mzmine_features, ppm = 10,
  section_id = "Section01", section_order = 1)
load_serial_msi_target_features(imzml_paths, mzmine_features, ppm = 10,
  section_ids = NULL, section_order = NULL)
make_feature_mapping(mzmine_features)
nearest_intensity_ppm(spectrum_mz, spectrum_intensity, targets, ppm = 10)
feature_stats(pixel_matrix)
select_features(pixel_matrix, feature_mapping = NULL, cv_top_percent = 70,
  mean_min = 0, nonzero_min = 0.3, manual_columns = character(),
  combine_mode = c("union", "intersection"))
select_shared_features(pixel_matrix, feature_mapping = NULL,
  section_column = "section_id", min_section_fraction = 1,
  cv_top_percent = 70, mean_min = 0, nonzero_min = 0.3,
  manual_columns = character(), combine_mode = c("union", "intersection"))
reduce_pixel_matrix(pixel_matrix, selected_features)
preprocess_select_features(pixel_matrix, feature_mapping = NULL,
  serial = NULL, section_column = "section_id", min_section_fraction = 1,
  cv_top_percent = 70, mean_min = 0, nonzero_min = 0.3,
  manual_columns = character(), combine_mode = c("union", "intersection"),
  do_background = TRUE,
  background_method = c("log_tic_kmeans", "tic_percentile", "median_fraction"),
  background_percent = 1, background_percentile = 10,
  min_nonzero_features = 1, do_tic = TRUE, do_log = TRUE)
preprocess_matrix(pixel_matrix, do_background = TRUE,
  background_method = c("log_tic_kmeans", "tic_percentile", "median_fraction"),
  background_percent = 1, background_percentile = 10,
  min_nonzero_features = 1, do_tic = TRUE, do_log = TRUE, do_scale = FALSE)
preprocess_sections(pixel_matrix, section_column = "section_id",
  do_background = TRUE,
  background_method = c("log_tic_kmeans", "tic_percentile", "median_fraction"),
  background_percent = 1, background_percentile = 10,
  min_nonzero_features = 1, do_tic = TRUE, do_log = TRUE,
  do_scale = FALSE, scale_scope = c("global", "section"))
cluster_pixels(pixel_matrix, k = 3, nstart = 25, pca_components = NULL)
cluster_diagnostics(pixel_matrix, max_k = 10, pca_components = NULL)
sample_subregions(clustered_matrix, grid_size = 5, min_pixels = 30)
sample_matched_regions(pixel_matrix, region_column = "matched_region_label",
  section_column = "section_id", section_order_column = "section_order",
  min_pixels = 30)
make_metaboanalyst_data(sample_matrix, group_column = NULL,
  group_prefix = NULL)
detect_metaboanalyst_result(result_data)
normalize_metaboanalyst_result(result_data, source_name = "")
backmap_sample_scores(score_data, sample_mapping, score_column)
plot_ion_image(pixel_matrix, feature, transform = c("identity", "log10"),
  title = NULL)
plot_feature_thumbnails(pixel_matrix, page = 1, per_page = 16,
  features = NULL, transform = c("identity", "log10"))
plot_feature_region_boxplot(sample_matrix, feature, region_column = "cluster",
  region_prefix = "Region_", title = NULL)
plot_metabo_feature_rank(result_data, top_n = 20, metric = NULL,
  decreasing = NULL, source_name = "")
plot_metabo_feature_view(result_data, pixel_matrix, feature = NULL, rank = 1,
  metric = NULL, top_n = 20, source_name = "",
  transform = c("identity", "log10"))
plot_metabo_score_plot(score_data, sample_matrix = NULL, x = NULL, y = NULL,
  source_name = "")
plot_sample_score_map(score_data, sample_mapping, pixel_matrix,
  score_column = NULL, source_name = "")
run_spatial_app(launch.browser = TRUE, port = NULL)
write_pipeline_outputs(output_dir, pixel_matrix = NULL, feature_mapping = NULL,
  section_mapping = NULL, background_stats = NULL, selected_features = NULL,
  reduced_matrix = NULL, preprocessed_matrix = NULL,
  clustered_matrix = NULL, sample_matrix = NULL, sample_mapping = NULL,
  metaboanalyst_data = NULL)
write_validation_plots(output_dir, pixel_matrix = NULL, sample_matrix = NULL,
  sample_mapping = NULL, metaboanalyst_data = NULL, feature = NULL,
  width = 8, height = 6, dpi = 150)
}
\arguments{
\item{data, pixel_matrix, clustered_matrix, sample_matrix}{Data frames produced by earlier pipeline steps. Pixel-level matrices must include \code{pixel_id}, \code{x}, \code{y}, and feature columns named \code{mz_*}. For MetaboAnalyst feature views, pass the same normalized/log-transformed pixel matrix used to construct the MetaboAnalyst input, usually the reduced matrix from \code{preprocess_select_features()}.}
\item{region_labels}{Matched-region label table containing \code{matched_region_label} plus \code{pixel_id}, \code{section_id + local_pixel_id}, or \code{section_id + x + y}.}
\item{path, paths}{Peak-picked MSI CSV path(s) with \code{x}, \code{y}, and numeric m/z columns.}
\item{sample_id, sample_ids}{Source sample labels stored with imported peak-picked MSI pixels.}
\item{imzml_path}{Path to an imzML file. The matching ibd file should be next to it.}
\item{imzml_paths}{Character vector of serial-section imzML file paths.}
\item{mzmine_features}{MZmine feature list with at least \code{id} and \code{mz} columns.}
\item{ppm}{Parts-per-million matching tolerance for target m/z extraction. Defaults to \code{10}, matching the pipeline outline.}
\item{section_id, section_ids, section_order}{Section metadata stored with extracted pixels.}
\item{spectrum_mz, spectrum_intensity}{Numeric vectors for one spectrum.}
\item{targets}{Numeric vector of target m/z values.}
\item{feature_mapping}{Feature mapping table from \code{make_feature_mapping()} or \code{infer_feature_mapping()}.}
\item{cv_top_percent, mean_min, nonzero_min}{Rule-based feature selection thresholds.}
\item{section_column, min_section_fraction}{Serial-section feature selection settings. Shared features must be selected in at least this fraction of sections.}
\item{serial}{Whether to use serial-section preprocessing and shared feature selection. If \code{NULL}, this is inferred from whether \code{section_column} is present and contains more than one unique section.}
\item{manual_columns}{Feature column names manually selected by the user.}
\item{combine_mode}{Whether manual and rule selections are combined by union or intersection.}
\item{selected_features}{Output from \code{select_features()}.}
\item{do_background, do_tic, do_log, do_scale}{Logical switches for preprocessing steps. Auto-scaling defaults to \code{FALSE} so exported sub-region samples can be scaled at the MetaboAnalystR sample level.}
\item{background_method}{Background removal method.}
\item{background_percent, background_percentile, min_nonzero_features}{Background removal thresholds.}
\item{k, nstart, max_k}{K-means clustering and diagnostic settings.}
\item{pca_components}{Optional number of leading principal components used for k-means and clustering diagnostics. If \code{NULL}, k-means uses all feature columns directly. PCA is mean-centered but not auto-scaled.}
\item{scale_scope}{Serial-section auto-scaling scope. \code{"global"} scales after section-wise preprocessing and merging; \code{"section"} scales each section independently.}
\item{grid_size, min_pixels}{Sub-region grid and filtering settings. Grid outputs include observed x/y spans for reporting effective grid dimensions.}
\item{section_order_column}{Column used to preserve serial-section order during matched-region sampling.}
\item{group_column, group_prefix}{Column and optional prefix used to create MetaboAnalyst group labels.}
\item{result_data}{A MetaboAnalyst result data frame.}
\item{source_name}{Optional uploaded file name used to help infer result format.}
\item{score_data}{PCA or PLS-DA score table with a \code{Sample} column.}
\item{sample_mapping}{Sub-region to pixel mapping table.}
\item{score_column}{Score column to map back to pixels, such as \code{PC1} or \code{Comp1}.}
\item{feature}{A feature column name such as \code{mz_834.5291}, or a numeric m/z value.}
\item{transform}{Plot intensity transform.}
\item{title}{Optional plot title.}
\item{page, per_page}{Thumbnail pagination settings.}
\item{features}{Optional feature vector for thumbnail plotting.}
\item{region_column}{Column used as a matched-region label in \code{sample_matched_regions()} or as the region/group label in \code{plot_feature_region_boxplot()}.}
\item{region_prefix}{Prefix added to cluster labels when \code{region_column = "cluster"}.}
\item{top_n}{Number of ranked MetaboAnalyst features to show.}
\item{rank}{Ranked feature index to show when \code{feature} is omitted.}
\item{metric}{Metric column to rank by, such as \code{VIP}, \code{FDR}, or \code{p.value}.}
\item{decreasing}{Whether larger metric values should rank first. If omitted, VIP ranks high-to-low and p-value/FDR ranks low-to-high.}
\item{x, y}{Score columns for the score plot, such as \code{PC1} and \code{PC2}.}
\item{output_dir}{Directory where CSV outputs should be written.}
\item{width, height, dpi}{PNG validation plot size and resolution.}
\item{section_mapping}{Optional serial-section metadata table.}
\item{background_stats}{Optional background-removal QC table returned by preprocessing functions.}
\item{launch.browser, port}{Arguments passed to \code{shiny::runApp()}.}
\item{reduced_matrix, preprocessed_matrix}{Optional intermediate matrices to write.}
\item{metaboanalyst_data}{MetaboAnalyst-ready data table.}
}
\value{
Depending on the function, returns a data frame, a list of data frames, a
preprocessing result list, a clustering result list, or a detected result type.
}
\details{
The package is intentionally scriptable. A typical workflow is:
load targeted m/z intensities from imzML, preprocess pixels without auto-scaling,
select spatially informative features on the normalized/log-transformed matrix,
reduce the matrix, cluster pixels with optional PCA scores, summarize cluster-by-grid
sub-regions, export to MetaboAnalyst, then map imported MetaboAnalyst results
back to features or pixels.

\code{preprocess_select_features()} is the recommended helper for the
preprocessing-to-selection boundary. It performs background removal, TIC
normalization, and log transform with \code{do_scale = FALSE}, then runs
\code{select_features()} or \code{select_shared_features()} and returns both the
normalized full matrix and the reduced selected-feature matrix. This keeps
CV-based feature selection from being computed on raw loading differences or on
auto-scaled features.

\code{sample_subregions()} and \code{sample_matched_regions()} summarize each
region with column means of the current feature matrix. In the default workflow
that matrix is already TIC-normalized and log-transformed, so the exported
sample intensities are arithmetic means on the log scale, equivalent to
geometric means of normalized intensities on the original scale.
\code{sample_matched_regions()} escapes \code{section_id} and
\code{matched_region_label} before combining them into \code{sample_id}, so
underscores or other separators in user labels cannot silently collide. The
original section and region labels remain available as separate columns.

\code{cluster_pixels(pca_components = ...)} and
\code{cluster_diagnostics(pca_components = ...)} run k-means on mean-centered
PCA scores without feature auto-scaling. Use the same \code{pca_components}
setting for diagnostics and final clustering when using diagnostics to choose
\code{k}. Diagnostics fit PCA on the full matrix before any pixel subsampling,
so the PC space matches final clustering. If Cardinal returns a square
intensity matrix where both dimensions equal the pixel count, the extractor
warns and assumes a features-by-pixels orientation. In serial-section feature
selection, sections that yield no selected features are skipped with a warning
before shared-feature fractions are computed. Region-label joins by
\code{section_id + x + y} drop duplicate label keys after warning, preserving
one row per pixel.

\code{make_feature_mapping()} writes m/z feature columns with high precision and
stops if duplicate column names would be generated. \code{sample_subregions()}
coerces \code{grid_size} to at least 2. \code{backmap_sample_scores()} stops
when no MetaboAnalyst samples match the sample mapping, instead of returning an
empty spatial map. \code{plot_sample_score_map()} uses the median mapped score
as the midpoint of its diverging color scale, avoiding misleading zero-centered
maps when imported PCA/PLS-DA scores are offset from zero.

For serial-section studies, use \code{load_serial_msi_target_features()} to
combine section-wise extractions, \code{select_shared_features()} to keep a
comparable feature set across sections, avoid pre-scaling by default so
MetaboAnalystR can scale exported samples, add a matched region label from
manual ROI annotation or registration, and use \code{sample_matched_regions()}
to construct section-by-region samples for MetaboAnalyst. Preprocessing returns
\code{background_stats} for QC of background removal.
}
\keyword{utilities}
