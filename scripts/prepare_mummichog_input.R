# Prepare all-feature mummichog input from spatial subregion samples.
#
# This script reuses the existing all-cluster sample mapping, but averages all
# 137 preprocessed mz features instead of only the feature-selected subset.
#
# R console usage:
# setwd("/Users/ly/Documents/Spatial Omics")
# source("scripts/prepare_mummichog_input.R")

if (file.exists("R/msi_pipeline.R")) {
  source("R/msi_pipeline.R")
} else {
  library(SpatialOmicsMSI)
}

if (!exists("data_dir", inherits = TRUE)) {
  data_dir <- "/Users/ly/Desktop/Jeff Xia/rat_brain_data"
}
if (!exists("test_out_dir", inherits = TRUE)) {
  test_out_dir <- file.path(data_dir, "spatial_test_outputs")
}
if (!exists("analysis_out_dir", inherits = TRUE)) {
  analysis_out_dir <- file.path(data_dir, "pathway_inputs")
}
dir.create(analysis_out_dir, recursive = TRUE, showWarnings = FALSE)

preprocessed_csv <- file.path(test_out_dir, "preprocessed_matrix.csv")
sample_mapping_csv <- file.path(test_out_dir, "metaboanalyst_all_clusters_sample_mapping.csv")
feature_mapping_csv <- file.path(test_out_dir, "feature_mapping.csv")
annotation_validation_csv <- file.path(data_dir, "annotation_cross_validation", "annotation_cross_validation.csv")

for (path in c(preprocessed_csv, sample_mapping_csv, feature_mapping_csv)) {
  if (!file.exists(path)) stop("Missing required file: ", path, call. = FALSE)
}

preprocessed <- read.csv(preprocessed_csv, check.names = FALSE, stringsAsFactors = FALSE)
sample_mapping <- read.csv(sample_mapping_csv, check.names = FALSE, stringsAsFactors = FALSE)
feature_mapping <- read.csv(feature_mapping_csv, check.names = FALSE, stringsAsFactors = FALSE)
annotation_validation <- if (file.exists(annotation_validation_csv)) {
  read.csv(annotation_validation_csv, check.names = FALSE, stringsAsFactors = FALSE)
} else {
  NULL
}

fcols <- feature_columns(preprocessed)
if (length(fcols) == 0) stop("No mz feature columns found in preprocessed matrix.", call. = FALSE)
required_columns(preprocessed, "pixel_id", "Preprocessed matrix")
required_columns(sample_mapping, c("sample_id", "roi_id", "pixel_ids"), "Sample mapping")

pixel_lookup <- seq_len(nrow(preprocessed))
names(pixel_lookup) <- as.character(preprocessed$pixel_id)

message("1. Averaging all ", length(fcols), " preprocessed mz features into spatial samples")
sample_rows <- lapply(seq_len(nrow(sample_mapping)), function(i) {
  ids <- strsplit(as.character(sample_mapping$pixel_ids[i]), ",", fixed = TRUE)[[1]]
  rows <- unname(pixel_lookup[ids])
  rows <- rows[!is.na(rows)]
  if (length(rows) == 0) return(NULL)

  means <- colMeans(preprocessed[rows, fcols, drop = FALSE], na.rm = TRUE)
  data.frame(
    sample_id = sample_mapping$sample_id[i],
    Group = sample_mapping$roi_id[i],
    n_pixels = length(rows),
    t(means),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
})
sample_rows <- sample_rows[!vapply(sample_rows, is.null, logical(1))]
all_feature_sample_matrix <- do.call(rbind, sample_rows)
rownames(all_feature_sample_matrix) <- NULL

message("2. Computing one-way ANOVA and group-level effect sizes for all features")
groups <- factor(all_feature_sample_matrix$Group)
if (length(levels(groups)) < 2) stop("Need at least two groups for ANOVA.", call. = FALSE)

stats_rows <- lapply(fcols, function(feature) {
  y <- as.numeric(all_feature_sample_matrix[[feature]])
  fit <- stats::aov(y ~ groups)
  aov_table <- summary(fit)[[1]]
  f_stat <- as.numeric(aov_table[["F value"]][1])
  p_value <- as.numeric(aov_table[["Pr(>F)"]][1])

  group_means <- tapply(y, groups, mean, na.rm = TRUE)
  max_group <- names(group_means)[which.max(group_means)]
  min_group <- names(group_means)[which.min(group_means)]
  max_mean <- max(group_means, na.rm = TRUE)
  min_mean <- min(group_means, na.rm = TRUE)

  # The preprocessed features are log10-transformed; convert max-min log10
  # difference to an approximate signed log2 fold-change for ranking.
  signed_log2fc <- (max_mean - min_mean) * log2(10)

  data.frame(
    feature = feature,
    mz = suppressWarnings(as.numeric(sub("^mz_", "", feature))),
    F.stat = f_stat,
    p.value = p_value,
    FDR = NA_real_,
    max_group = max_group,
    min_group = min_group,
    max_group_mean = max_mean,
    min_group_mean = min_mean,
    signed_log2fc = signed_log2fc,
    statistic = signed_log2fc,
    stringsAsFactors = FALSE
  )
})
all_feature_stats <- do.call(rbind, stats_rows)
all_feature_stats$FDR <- stats::p.adjust(all_feature_stats$p.value, method = "BH")

feature_meta <- feature_mapping[, intersect(c("column_name", "mzmine_id", "feature_id", "mz"), names(feature_mapping)), drop = FALSE]
names(feature_meta)[names(feature_meta) == "mz"] <- "feature_mapping_mz"
all_feature_stats <- merge(
  all_feature_stats,
  feature_meta,
  by.x = "feature",
  by.y = "column_name",
  all.x = TRUE,
  sort = FALSE
)

if (!is.null(annotation_validation)) {
  all_feature_stats <- merge(
    all_feature_stats,
    annotation_validation,
    by.x = "mzmine_id",
    by.y = "feature_id",
    all.x = TRUE,
    sort = FALSE
  )
}

mummichog_input <- all_feature_stats[, c(
  "mz", "p.value", "statistic", "signed_log2fc", "F.stat", "FDR",
  "feature", "mzmine_id", "consistency_class", "swiss_neutral_formulas",
  "swiss_names", "metaspace_formulas"
), drop = FALSE]
names(mummichog_input)[names(mummichog_input) == "mz"] <- "m.z"

all_feature_sample_path <- file.path(analysis_out_dir, "all_feature_sample_matrix.csv")
all_feature_stats_path <- file.path(analysis_out_dir, "all_feature_anova_effects.csv")
mummichog_csv_path <- file.path(analysis_out_dir, "mummichog_input_all_features.csv")
mummichog_tsv_path <- file.path(analysis_out_dir, "mummichog_input_all_features.tsv")

write.csv(all_feature_sample_matrix, all_feature_sample_path, row.names = FALSE)
write.csv(all_feature_stats, all_feature_stats_path, row.names = FALSE)
write.csv(mummichog_input, mummichog_csv_path, row.names = FALSE)
utils::write.table(mummichog_input, mummichog_tsv_path, sep = "\t", quote = FALSE, row.names = FALSE)

cat("\n=== Mummichog input preparation complete ===\n")
cat("Samples:", nrow(all_feature_sample_matrix), "\n")
cat("Features:", length(fcols), "\n")
cat("Groups:\n")
print(table(all_feature_sample_matrix$Group))
cat("Significant features at FDR < 0.05:", sum(all_feature_stats$FDR < 0.05, na.rm = TRUE), "\n")
cat("All-feature sample matrix:", normalizePath(all_feature_sample_path), "\n")
cat("ANOVA/effect table:", normalizePath(all_feature_stats_path), "\n")
cat("Mummichog CSV:", normalizePath(mummichog_csv_path), "\n")
cat("Mummichog TSV:", normalizePath(mummichog_tsv_path), "\n")
