# Cross-validate MZmine features, SwissLipids candidates, and METASPACE annotations.
#
# R console usage:
# setwd("/Users/ly/Documents/Spatial Omics")
# feature_csv <- "/Users/ly/Desktop/Jeff Xia/rat_brain_data/Filtered_7E4.csv"
# swisslipids_csv <- "/path/to/swisslipids_species_filtered.csv"
# metaspace_csv <- "/Users/ly/Desktop/Jeff Xia/rat_brain_data/annotations.csv"
# source("scripts/cross_validate_annotations.R")

if (file.exists("R/msi_pipeline.R")) {
  source("R/msi_pipeline.R")
} else {
  library(SpatialOmicsMSI)
}

if (!exists("feature_csv", inherits = TRUE)) {
  feature_csv <- "/Users/ly/Desktop/Jeff Xia/rat_brain_data/Filtered_7E4.csv"
}
if (!exists("metaspace_csv", inherits = TRUE)) {
  metaspace_csv <- "/Users/ly/Desktop/Jeff Xia/rat_brain_data/annotations.csv"
}
if (!exists("swisslipids_csv", inherits = TRUE)) {
  swisslipids_csv <- "/Users/ly/Downloads/swisslipids_species_filtered.csv"
}
if (!exists("ppm", inherits = TRUE)) ppm <- 5

if (!file.exists(feature_csv)) stop("Missing MZmine feature CSV: ", feature_csv, call. = FALSE)
if (!file.exists(swisslipids_csv)) {
  stop(
    "Missing SwissLipids CSV: ", swisslipids_csv, "\n",
    "Set `swisslipids_csv` before sourcing this script.",
    call. = FALSE
  )
}
if (!file.exists(metaspace_csv)) stop("Missing METASPACE annotation CSV: ", metaspace_csv, call. = FALSE)

out_dir <- file.path(dirname(feature_csv), "annotation_cross_validation")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

validation <- cross_validate_annotations(
  feature_csv = feature_csv,
  swisslipids_csv = swisslipids_csv,
  metaspace_csv = metaspace_csv,
  ppm = ppm
)

summary_table <- attr(validation, "summary")
raw_summary_table <- as.data.frame(table(validation$raw_consistency_class), stringsAsFactors = FALSE)
names(raw_summary_table) <- c("raw_consistency_class", "feature_count")
class_transition_table <- as.data.frame(
  table(raw_consistency_class = validation$raw_consistency_class, consistency_class = validation$consistency_class),
  stringsAsFactors = FALSE
)
class_transition_table <- class_transition_table[class_transition_table$Freq > 0, , drop = FALSE]

write.csv(validation, file.path(out_dir, "annotation_cross_validation.csv"), row.names = FALSE)
write.csv(summary_table, file.path(out_dir, "annotation_cross_validation_summary.csv"), row.names = FALSE)
write.csv(raw_summary_table, file.path(out_dir, "annotation_cross_validation_raw_summary.csv"), row.names = FALSE)
write.csv(class_transition_table, file.path(out_dir, "annotation_cross_validation_class_transition.csv"), row.names = FALSE)

cat("\n=== Annotation cross-validation complete ===\n")
cat("MZmine features:", normalizePath(feature_csv), "\n")
cat("SwissLipids candidates:", normalizePath(swisslipids_csv), "\n")
cat("METASPACE annotations:", normalizePath(metaspace_csv), "\n")
cat("PPM tolerance:", ppm, "\n")
cat("Classification summary:\n")
print(summary_table)
cat("Raw classification summary before neutral formula correction:\n")
print(raw_summary_table)
cat("Class transition after neutral formula correction:\n")
print(class_transition_table)
cat("Output directory:", normalizePath(out_dir), "\n")
