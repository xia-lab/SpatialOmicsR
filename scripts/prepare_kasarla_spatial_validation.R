args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("Usage: prepare_kasarla_spatial_validation.R table_s12.csv output_dir", call. = FALSE)
source <- read.csv(args[1], check.names = FALSE, stringsAsFactors = FALSE)
output_dir <- args[2]
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

metabolite <- trimws(source[[1]])
feature <- make.unique(paste0("metabolite_", gsub("[^a-z0-9]+", "_", tolower(metabolite))))
lcm_columns <- grep("LMD-ROI", names(source), fixed = TRUE)
msi_columns <- grep("MALDI-ROI", names(source), fixed = TRUE)
region_of <- function(columns) {
  label <- tolower(names(source)[columns])
  ifelse(grepl("cortex", label), "cortex",
    ifelse(grepl("medulla", label), "medulla",
      ifelse(grepl("pelvis", label), "renal_pelvis", NA_character_)))
}
mean_by_region <- function(columns) {
  regions <- region_of(columns)
  values <- lapply(c("cortex", "medulla", "renal_pelvis"), function(region) {
    selected <- columns[regions == region]
    matrix <- sapply(source[selected], function(x) suppressWarnings(as.numeric(x)))
    if (is.null(dim(matrix))) matrix <- matrix(matrix, ncol = 1L)
    rowMeans(matrix, na.rm = TRUE)
  })
  result <- data.frame(roi_id = c("cortex", "medulla", "renal_pelvis"),
    section_id = "kasarla2025_kidney_serial_sections", stringsAsFactors = FALSE)
  for (i in seq_along(feature)) result[[feature[i]]] <- vapply(values, `[`, numeric(1), i)
  result
}

lcm <- mean_by_region(lcm_columns)
msi <- mean_by_region(msi_columns)
mapping <- data.frame(
  metabolite_name = metabolite,
  msi_feature = feature,
  lcm_feature = feature,
  evidence = "Kasarla 2025 Table S12; region means across reported technical/ROI replicates",
  stringsAsFactors = FALSE)

write.csv(msi, file.path(output_dir, "kasarla_kidney_maldi_roi_matrix.csv"), row.names = FALSE)
write.csv(lcm, file.path(output_dir, "kasarla_kidney_lmd_lcms_roi_matrix.csv"), row.names = FALSE)
write.csv(mapping, file.path(output_dir, "kasarla_kidney_feature_mapping.csv"), row.names = FALSE)
cat(sprintf("Prepared %d matched metabolites across 3 kidney regions.\n", nrow(mapping)))
