#!/usr/bin/env Rscript

# Create the compact, deterministic real-data fixture used by the executable
# package vignette. The source matrix is reconstructed from public METASPACE
# annotation ion images for Kasarla et al. (2025); see inst/extdata/kasarla2025/SOURCE.md.

root <- normalizePath(Sys.getenv("SPATIALOMICS_PROJECT_ROOT", unset = getwd()))
source_rds <- file.path(root, "data_raw/kasarla2025/processed/kasarla_kidney_pixel_feature_matrix.rds")
source_csv <- file.path(root, "data_raw/kasarla2025/processed/kasarla_kidney_pixel_feature_matrix.csv.gz")
output_path <- file.path(root, "inst/extdata/kasarla2025/kasarla_kidney_vignette_pixels.csv.gz")
metadata_path <- file.path(root, "inst/extdata/kasarla2025/kasarla_kidney_vignette_features.csv")

if (file.exists(source_rds)) {
  data <- readRDS(source_rds)
} else if (file.exists(source_csv)) {
  data <- utils::read.csv(source_csv, check.names = FALSE)
} else {
  stop("Missing source matrix. Expected ", source_rds, " or ", source_csv, call. = FALSE)
}

# A fixed 50 x 50 measured field gives a connected grid while keeping package
# size and vignette build time modest.
data <- data[data$x >= 76 & data$x <= 125 & data$y >= 102 & data$y <= 151, , drop = FALSE]
stopifnot(nrow(data) == 2500L)

feature_columns <- grep("^mz_", names(data), value = TRUE)
raw <- as.matrix(data[feature_columns])
storage.mode(raw) <- "double"
tic <- rowSums(raw, na.rm = TRUE)
analysis <- log10(raw / pmax(tic, .Machine$double.eps) + 1)
detection <- colMeans(raw > 0, na.rm = TRUE)
variance <- apply(analysis, 2, stats::var, na.rm = TRUE)
eligible <- detection > 0.05 & detection < 0.95 & is.finite(variance) & variance > 0
ranked <- names(sort(variance[eligible], decreasing = TRUE))
if (length(ranked) < 24L) {
  stop("Source matrix has fewer than 24 eligible real features.", call. = FALSE)
}
selected <- ranked[seq_len(24L)]

fixture <- data.frame(
  pixel_id = seq_len(nrow(data)),
  x = data$x - min(data$x) + 1,
  y = data$y - min(data$y) + 1,
  data[selected],
  check.names = FALSE
)
dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
connection <- gzfile(output_path, open = "wt")
on.exit(close(connection), add = TRUE)
utils::write.csv(fixture, connection, row.names = FALSE)
close(connection)
on.exit(NULL, add = FALSE)

source_metadata <- utils::read.csv(
  file.path(root, "data_raw/kasarla2025/processed/kasarla_kidney_metaspace_feature_metadata.csv"),
  check.names = FALSE, stringsAsFactors = FALSE
)
feature_metadata <- source_metadata[match(selected, source_metadata$column_name), , drop = FALSE]
if (anyNA(feature_metadata$column_name)) {
  stop("Selected fixture features are missing from source metadata.", call. = FALSE)
}
feature_metadata$detection_fraction_in_fixture <- unname(detection[selected])
feature_metadata$variance_after_tic_log10 <- unname(variance[selected])
utils::write.csv(feature_metadata, metadata_path, row.names = FALSE)

message("Wrote ", nrow(fixture), " pixels x ", length(selected), " features to ", output_path)
