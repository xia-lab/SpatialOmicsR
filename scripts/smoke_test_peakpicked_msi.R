# ============================================================
# SpatialOmicsMSI full workflow smoke test
# ============================================================
#
# Terminal usage:
# Rscript scripts/smoke_test_peakpicked_msi.R path/to/your_data.csv
#
# R console usage:
# input_csv <- "path/to/your_data.csv"
# source("scripts/smoke_test_peakpicked_msi.R")

source("scripts/_bootstrap.R")
load_spatialomics_code()

args <- commandArgs(trailingOnly = TRUE)
if (!exists("input_csv", inherits = TRUE) || !nzchar(input_csv)) {
  input_csv <- if (length(args) >= 1) args[[1]] else ""
}

if (!nzchar(input_csv)) {
  stop(
    "Please provide an input CSV path.\n",
    "Terminal usage: Rscript scripts/smoke_test_peakpicked_msi.R path/to/your_data.csv\n",
    "R console usage: input_csv <- 'path/to/your_data.csv'; source('scripts/smoke_test_peakpicked_msi.R')",
    call. = FALSE
  )
}
if (!file.exists(input_csv)) {
  stop("Input CSV does not exist: ", input_csv, call. = FALSE)
}

# ------ Step 1: Import ------
pixel_matrix <- import_peakpicked_msi_csv(input_csv)

cat("=== Step 1: Import ===\n")
cat("Input CSV:", normalizePath(input_csv), "\n")
cat("Dimensions:", dim(pixel_matrix), "\n")
cat("Feature columns:", length(feature_columns(pixel_matrix)), "\n")
cat("Pixel count:", nrow(pixel_matrix), "\n")
cat("x range:", range(pixel_matrix$x), "\n")
cat("y range:", range(pixel_matrix$y), "\n\n")

# ------ Step 2: Feature mapping ------
feature_mapping2 <- infer_feature_mapping(pixel_matrix)

cat("=== Step 2: Feature Mapping ===\n")
cat("Features in mapping:", nrow(feature_mapping2), "\n")
print(head(feature_mapping2))
cat("\n")

# ------ Step 3: Preprocessing + feature selection ------
warnings_step3 <- character(0)

withCallingHandlers(
  {
    feature_selection2 <- preprocess_select_features(
      pixel_matrix,
      feature_mapping = feature_mapping2,
      background_method = "gap_otsu_log_tic",
      ubiquitous_method = "gap",
      do_background = TRUE,
      do_tic = TRUE,
      do_log = TRUE,
      cv_top_percent = 70,
      mean_min = 0,
      nonzero_min = 0.2
    )
  },
  warning = function(w) {
    warnings_step3 <<- c(warnings_step3, conditionMessage(w))
    invokeRestart("muffleWarning")
  }
)

cat("=== Step 3: Preprocess + Feature Selection ===\n")
if (length(warnings_step3) == 0) {
  cat("warnings(): CLEAN - no Step 3 warnings\n")
} else {
  cat("WARNINGS detected in Step 3:\n")
  for (msg in warnings_step3) cat(" -", msg, "\n")
}

cat("\nBackground stats:\n")
stat_cols <- c("threshold", "pixels_before", "pixels_after", "removed_fraction")
stat_cols <- intersect(stat_cols, names(feature_selection2$background_stats))
print(feature_selection2$background_stats[, stat_cols, drop = FALSE])

cat("\nSelected features:", nrow(feature_selection2$selected_features), "\n")
print(feature_selection2$selected_features$column_name)

cat("\nreduced_matrix dimensions:", dim(feature_selection2$reduced_matrix), "\n")
cat("reduced_matrix columns:", names(feature_selection2$reduced_matrix), "\n\n")

# ------ Step 4: Clustering ------
stopifnot(nrow(feature_selection2$selected_features) > 0)
stopifnot(length(feature_columns(feature_selection2$reduced_matrix)) > 0)
stopifnot(nrow(feature_selection2$reduced_matrix) >= 3)

cluster_result <- cluster_pixels(
  feature_selection2$reduced_matrix,
  k = 3,
  seed = 1
)

cat("=== Step 4: Clustering (k=3) ===\n")
cat("Cluster table:\n")
print(table(cluster_result$matrix$cluster))
cat("\n")

# ------ Step 5: Quick visualization check ------
cat("=== Step 5: Cluster Map ===\n")
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop("Package 'ggplot2' is required for plotting.", call. = FALSE)
}

plot_data <- cluster_result$matrix
plot_data$cluster <- as.factor(plot_data$cluster)

p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = x, y = -y, color = cluster)) +
  ggplot2::geom_point(size = 0.3) +
  ggplot2::scale_color_manual(values = c("1" = "#E41A1C", "2" = "#377EB8", "3" = "#4DAF4A")) +
  ggplot2::coord_equal() +
  ggplot2::theme_minimal() +
  ggplot2::labs(title = "Pixel Clustering (k=3)", color = "Cluster")

print(p)

dir.create("outputs", recursive = TRUE, showWarnings = FALSE)
plot_path <- file.path("outputs", "test_cluster_map.png")
ggplot2::ggsave(plot_path, p, width = 6, height = 6, dpi = 300)
cat("Cluster map saved to", normalizePath(plot_path), "\n")
