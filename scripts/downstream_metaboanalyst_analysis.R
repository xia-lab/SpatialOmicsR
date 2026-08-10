# Downstream MetaboAnalystR statistical analysis
#
# Run from the repository root after setting SPATIALOMICS_RAT_BRAIN_DIR.
#
# Default: build an all-cluster MetaboAnalystR input from spatial_test_outputs/clustered_matrix.csv
# source("scripts/downstream_metaboanalyst_analysis.R")
#
# Or provide your own MetaboAnalyst-format CSV:
# input_csv <- "/path/to/metaboanalyst_data.csv"
# output_dir <- "/path/to/metaboanalyst_results"
# source("scripts/downstream_metaboanalyst_analysis.R")

source("scripts/_bootstrap.R")

required_pkgs <- c("MetaboAnalystR", "RSclient", "factoextra", "ggplot2", "pls")
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop(
    "Missing required package(s): ",
    paste(missing_pkgs, collapse = ", "),
    ". Install them before running this script.",
    call. = FALSE
  )
}

library(MetaboAnalystR)

run_local_pca <- function(mSet) {
  pca <- stats::prcomp(mSet$dataSet$norm, center = TRUE, scale. = FALSE)
  sum_pca <- summary(pca)
  imp_pca <- sum_pca$importance
  std_pca <- imp_pca[1, ]
  var_pca <- imp_pca[2, ]
  cum_pca <- imp_pca[3, ]
  contrib <- pca$rotation^2
  contrib <- sweep(contrib, 2, colSums(contrib), "/", check.margin = FALSE) * 100
  if (ncol(contrib) > 10) {
    contrib <- contrib[, 1:10, drop = FALSE]
  }
  mSet$analSet$pca <- append(pca, list(
    std = std_pca,
    variance = var_pca,
    cum.var = cum_pca,
    contrib = contrib,
    loading.type = "all"
  ))
  mSet$custom.cmpds <- c()
  utils::write.csv(signif(mSet$analSet$pca$x, 5), file = "pca_score.csv")
  utils::write.csv(signif(mSet$analSet$pca$rotation, 5), file = "pca_loadings.csv")
  mSet
}

plot_local_score <- function(scores, groups, x, y, output_file, title, x_label = x, y_label = y) {
  if (length(groups) != nrow(scores) && exists("input_data", inherits = TRUE)) {
    groups <- input_data$Group[match(rownames(scores), input_data$Sample)]
  }
  if (length(groups) != nrow(scores)) {
    groups <- rep("Sample", nrow(scores))
  }
  plot_data <- data.frame(
    Sample = rownames(scores),
    Group = groups,
    x = scores[, x],
    y = scores[, y],
    stringsAsFactors = FALSE
  )
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = .data[["x"]], y = .data[["y"]], color = .data[["Group"]])) +
    ggplot2::geom_point(size = 2.7, alpha = 0.9) +
    ggplot2::theme_minimal() +
    ggplot2::labs(title = title, x = x_label, y = y_label, color = "Group")
  ggplot2::ggsave(output_file, p, width = 7, height = 6, dpi = 300)
  invisible(output_file)
}

plot_local_heatmap <- function(norm_data, groups, output_file, max_features = 50) {
  feature_var <- apply(norm_data, 2, stats::var, na.rm = TRUE)
  keep_features <- names(sort(feature_var, decreasing = TRUE))[seq_len(min(max_features, length(feature_var)))]
  plot_matrix <- norm_data[, keep_features, drop = FALSE]
  scaled_matrix <- scale(plot_matrix)
  sample_order <- order(groups, rownames(norm_data))
  ordered_samples <- rownames(norm_data)[sample_order]
  heatmap_data <- as.data.frame(as.table(as.matrix(scaled_matrix)))
  names(heatmap_data) <- c("Sample", "Feature", "ScaledIntensity")
  heatmap_data$Group <- groups[match(heatmap_data$Sample, rownames(norm_data))]
  heatmap_data$Sample <- factor(heatmap_data$Sample, levels = ordered_samples)
  heatmap_data$Feature <- factor(heatmap_data$Feature, levels = rev(keep_features))
  group_sizes <- as.integer(table(factor(groups[sample_order], levels = unique(groups[sample_order]))))
  boundary_positions <- cumsum(group_sizes)[-length(group_sizes)] + 0.5
  label_positions <- cumsum(group_sizes) - group_sizes / 2 + 0.5
  label_data <- data.frame(
    x = label_positions,
    y = length(keep_features) + 1,
    Group = unique(groups[sample_order]),
    stringsAsFactors = FALSE
  )
  p <- ggplot2::ggplot(heatmap_data, ggplot2::aes(x = .data[["Sample"]], y = .data[["Feature"]], fill = .data[["ScaledIntensity"]])) +
    ggplot2::geom_tile() +
    ggplot2::geom_vline(xintercept = boundary_positions, color = "#333333", linewidth = 0.3) +
    ggplot2::geom_text(
      data = label_data,
      ggplot2::aes(x = .data[["x"]], y = .data[["y"]], label = .data[["Group"]]),
      inherit.aes = FALSE,
      size = 3.2
    ) +
    ggplot2::scale_fill_gradient2(low = "#2563eb", mid = "white", high = "#dc2626", midpoint = 0) +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(t = 20, r = 10, b = 10, l = 10)
    ) +
    ggplot2::labs(title = "Normalized feature heatmap", x = "Samples", y = "Features", fill = "Z-score")
  ggplot2::ggsave(output_file, p, width = 10, height = 8, dpi = 300)
  invisible(output_file)
}

plot_local_norm_summary <- function(raw_data, norm_data, output_file, max_features = 30) {
  feature_var <- apply(raw_data, 2, stats::var, na.rm = TRUE)
  keep_features <- names(sort(feature_var, decreasing = TRUE))[seq_len(min(max_features, length(feature_var)))]
  raw_long <- as.data.frame(as.table(as.matrix(raw_data[, keep_features, drop = FALSE])))
  norm_long <- as.data.frame(as.table(as.matrix(norm_data[, keep_features, drop = FALSE])))
  names(raw_long) <- c("Sample", "Feature", "Value")
  names(norm_long) <- c("Sample", "Feature", "Value")
  raw_long$Stage <- "Before normalization"
  norm_long$Stage <- "After normalization"
  plot_data <- rbind(raw_long, norm_long)
  plot_data$Feature <- factor(plot_data$Feature, levels = rev(keep_features))

  p_density <- ggplot2::ggplot(plot_data, ggplot2::aes(x = .data[["Value"]], color = .data[["Stage"]])) +
    ggplot2::geom_density(linewidth = 0.8, na.rm = TRUE) +
    ggplot2::facet_wrap(~Stage, scales = "free") +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "none") +
    ggplot2::labs(title = "Feature intensity distribution", x = "Value", y = "Density")

  p_box <- ggplot2::ggplot(plot_data, ggplot2::aes(x = .data[["Value"]], y = .data[["Feature"]])) +
    ggplot2::geom_boxplot(fill = "#9ae68f", outlier.size = 0.6) +
    ggplot2::facet_wrap(~Stage, scales = "free_x") +
    ggplot2::theme_minimal() +
    ggplot2::labs(x = "Value", y = "Features")

  grDevices::png(output_file, width = 900, height = 1000, res = 120)
  grid::grid.newpage()
  grid::pushViewport(grid::viewport(layout = grid::grid.layout(2, 1, heights = grid::unit(c(0.32, 0.68), "npc"))))
  print(p_density, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 1))
  print(p_box, vp = grid::viewport(layout.pos.row = 2, layout.pos.col = 1))
  grDevices::dev.off()
  invisible(output_file)
}

plot_local_sample_norm_summary <- function(raw_data, norm_data, groups, output_file) {
  sample_order <- order(groups, rownames(norm_data))
  ordered_samples <- rownames(norm_data)[sample_order]
  raw_means <- rowMeans(raw_data[ordered_samples, , drop = FALSE], na.rm = TRUE)
  norm_means <- rowMeans(norm_data[ordered_samples, , drop = FALSE], na.rm = TRUE)
  plot_data <- rbind(
    data.frame(Sample = ordered_samples, Group = groups[sample_order], Stage = "Before normalization", Mean = raw_means),
    data.frame(Sample = ordered_samples, Group = groups[sample_order], Stage = "After normalization", Mean = norm_means)
  )
  plot_data$Sample <- factor(plot_data$Sample, levels = ordered_samples)
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = .data[["Sample"]], y = .data[["Mean"]], fill = .data[["Group"]])) +
    ggplot2::geom_col(width = 0.85) +
    ggplot2::facet_wrap(~Stage, scales = "free_y", ncol = 1) +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.text.x = ggplot2::element_blank(), axis.ticks.x = ggplot2::element_blank()) +
    ggplot2::labs(title = "Sample mean intensity summary", x = "Samples", y = "Mean intensity", fill = "Group")
  ggplot2::ggsave(output_file, p, width = 10, height = 7, dpi = 300)
  invisible(output_file)
}

load_spatialomics_code()

data_dir <- spatialomics_data_dir(
  "SPATIALOMICS_RAT_BRAIN_DIR", "data_raw/rat_brain_data", "Rat-brain data"
)
test_out_dir <- file.path(data_dir, "spatial_test_outputs")

if (!exists("output_dir", inherits = TRUE)) {
  output_dir <- file.path(test_out_dir, "metaboanalyst_results")
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

if (!exists("input_csv", inherits = TRUE)) {
  clustered_path <- file.path(test_out_dir, "clustered_matrix.csv")
  if (!file.exists(clustered_path)) {
    stop(
      "Set input_csv or run scripts/test_rat_brain_workflow.R first to create clustered_matrix.csv.",
      call. = FALSE
    )
  }
  clustered <- read.csv(clustered_path, check.names = FALSE)
  all_cluster_rois <- define_rois(clustered, mode = "cluster")
  all_cluster_samples <- sample_subregions(all_cluster_rois, grid_size = 5, min_pixels = 30)
  metabo_input <- make_metaboanalyst_data(all_cluster_samples$sample_matrix)
  input_csv <- file.path(output_dir, "metaboanalyst_all_clusters_data.csv")
  write.csv(metabo_input, input_csv, row.names = FALSE)
  write.csv(
    all_cluster_samples$sample_mapping,
    file.path(output_dir, "metaboanalyst_all_clusters_sample_mapping.csv"),
    row.names = FALSE
  )
}

# Inputs produced by this package's preprocessing workflow have already been
# transformed with log10(x + 1). Do not apply MetaboAnalystR LogNorm to those
# values a second time. Set input_already_log_transformed <- FALSE before
# sourcing this script only when supplying an untransformed external input_csv.
if (!exists("input_already_log_transformed", inherits = TRUE)) {
  input_already_log_transformed <- TRUE
}
if (length(input_already_log_transformed) != 1L || is.na(input_already_log_transformed)) {
  stop("input_already_log_transformed must be TRUE or FALSE.", call. = FALSE)
}
metabo_transform <- if (isTRUE(input_already_log_transformed)) "NULL" else "LogNorm"

if (!file.exists(input_csv)) stop("Missing input file: ", input_csv, call. = FALSE)

input_data <- read.csv(input_csv, check.names = FALSE)
if (!"Group" %in% names(input_data)) {
  stop("Input CSV must contain a 'Group' column from make_metaboanalyst_data().", call. = FALSE)
}
groups <- unique(input_data$Group)
if (length(groups) < 2) {
  stop(
    "Only 1 group found in Group column (",
    groups,
    "). PLS-DA and ANOVA require >= 2 groups; use an all-cluster or multi-ROI export.",
    call. = FALSE
  )
}

old_wd <- getwd()
on.exit(setwd(old_wd), add = TRUE)
setwd(output_dir)

if (!exists("include_extra_metabo_plots", inherits = TRUE)) {
  include_extra_metabo_plots <- FALSE
}

unlink(c(
  "anova_dpi72.png",
  "heatmap_dpi72.png",
  "heatmap_local.png",
  "norm_summary_dpi72.png",
  "sample_norm_summary_dpi72.png",
  "pca_biplot_dpi72.png",
  "pca_loading_dpi72.png",
  "pca_pair_dpi72.png",
  "pca_score2d_dpi72.png",
  "pca_score2d_local.png",
  "pca_scree_dpi72.png",
  "pls_cv_dpi72.png",
  "pls_pair_dpi72.png",
  "pls_score2d_dpi72.png",
  "pls_score2d_local.png",
  "pls_vip_dpi72.png"
))

# Clear stale MetaboAnalystR state files. Normalization picks the newest
# preprocessing file, so an old broken data_proc.qs can shadow preproc.qs.
unlink(c(
  "data.edit.qs", "data.edit.qs2",
  "data.filt.qs", "data.filt.qs2",
  "data_proc.qs", "data_proc.qs2",
  "preproc.qs", "preproc.qs2",
  "preproc.orig.qs", "preproc.orig.qs2",
  "prenorm.qs", "prenorm.qs2",
  "row_norm.qs", "row_norm.qs2",
  "complete_norm.qs", "complete_norm.qs2",
  "data_orig.qs", "data_orig.qs2",
  "data_orig_0.qs", "data_orig_0.qs2"
))

# MetaboAnalystR message helpers expect this global in some execution paths.
current.msg <- ""

mSet <- InitDataObjects("conc", "stat", FALSE, default.dpi = 72)
mSet <- Read.TextData(mSet, normalizePath(input_csv), "rowu", "disc")
mSet <- SanityCheckData(mSet)

# MetaboAnalystR 4.3.0 normalization summary plots read data_proc.qs,
# while this import path writes preproc.qs/preproc.orig.qs. Mirror the
# checked preprocessing matrix so PlotNormSummary can run.
if (!file.exists("data_proc.qs") && file.exists("preproc.qs")) {
  MetaboAnalystR:::ov_qs_save(MetaboAnalystR:::ov_qs_read("preproc.qs"), "data_proc.qs")
}

mSet <- Normalization(
  mSet,
  rowNorm = "NULL",
  transNorm = metabo_transform,
  scaleNorm = "AutoNorm",
  ratio = FALSE,
  ratioNum = 20
)

write.csv(
  data.frame(
    input_already_log_transformed = isTRUE(input_already_log_transformed),
    MetaboAnalystR_transNorm = metabo_transform,
    scaleNorm = "AutoNorm",
    stringsAsFactors = FALSE
  ),
  file.path(output_dir, "metaboanalyst_normalization_parameters.csv"),
  row.names = FALSE
)

# MetaboAnalystR 4.3.0 can drop class labels after Normalization() for this
# row-wise CSV path. Re-align labels by sample name before PCA/PLS/ANOVA.
norm_samples <- rownames(mSet$dataSet$norm)
norm_groups <- input_data$Group[match(norm_samples, input_data$Sample)]
if (any(is.na(norm_groups))) {
  stop("Could not match normalized samples back to input Group labels.", call. = FALSE)
}
mSet$dataSet$cls <- factor(norm_groups)
mSet$dataSet$orig.cls <- mSet$dataSet$cls
mSet$dataSet$proc.cls <- mSet$dataSet$cls
mSet$dataSet$prenorm.cls <- mSet$dataSet$cls
mSet$dataSet$cls.num <- length(levels(mSet$dataSet$cls))
mSet$dataSet$min.grp.size <- min(table(mSet$dataSet$cls))
mSet$dataSet$meta.info <- data.frame(Class = mSet$dataSet$cls, row.names = norm_samples)

raw_proc_data <- MetaboAnalystR:::ov_qs_read("data_proc.qs")
plot_local_norm_summary(raw_proc_data, mSet$dataSet$norm, file.path(output_dir, "norm_summary_dpi72.png"))
plot_local_sample_norm_summary(
  raw_proc_data,
  mSet$dataSet$norm,
  mSet$dataSet$cls,
  file.path(output_dir, "sample_norm_summary_dpi72.png")
)

mSet <- run_local_pca(mSet)
plot_local_score(
  mSet$analSet$pca$x,
  mSet$dataSet$cls,
  "PC1",
  "PC2",
  file.path(output_dir, "pca_score2d_local.png"),
  "PCA score plot",
  x_label = paste0("PC1 (", round(100 * mSet$analSet$pca$variance[1], 1), "%)"),
  y_label = paste0("PC2 (", round(100 * mSet$analSet$pca$variance[2], 1), "%)")
)
mSet <- PlotPCAScree(mSet, "pca_scree_", "png", 72, width = NA, scree.num = 5)
mSet <- PlotPCALoading(mSet, "pca_loading_", "png", 72, width = NA, inx1 = 1, inx2 = 2)
if (include_extra_metabo_plots) {
  mSet <- tryCatch(
    PlotPCAPairSummary(mSet, "pca_pair_", "png", 72, width = NA, pc.num = 5),
    error = function(e) {
      warning("Skipping PCA pair summary plot: ", conditionMessage(e), call. = FALSE)
      mSet
    }
  )
  mSet <- tryCatch(
    PlotPCA2DScore(mSet, "pca_score2d_", "png", 72, width = NA, pcx = 1, pcy = 2, reg = 0.95, show = 1, grey.scale = 0),
    error = function(e) {
      warning("Skipping MetaboAnalystR PCA score plot: ", conditionMessage(e), call. = FALSE)
      mSet
    }
  )
  mSet <- tryCatch(
    PlotPCABiplot(mSet, "pca_biplot_", "png", 72, width = NA, inx1 = 1, inx2 = 2, topnum = 10),
    error = function(e) {
      warning("Skipping PCA biplot: ", conditionMessage(e), call. = FALSE)
      mSet
    }
  )
}

mSet <- PLSR.Anal(mSet, reg = TRUE)
if (include_extra_metabo_plots) {
  mSet <- tryCatch(
    PlotPLSPairSummary(mSet, "pls_pair_", "png", 72, width = NA, pc.num = 5),
    error = function(e) {
      warning("Skipping PLS pair summary plot: ", conditionMessage(e), call. = FALSE)
      mSet
    }
  )
}
plot_local_score(
  mSet$analSet$plsr$scores,
  mSet$dataSet$cls,
  "Comp 1",
  "Comp 2",
  file.path(output_dir, "pls_score2d_dpi72.png"),
  "PLS-DA score plot",
  x_label = "Component 1",
  y_label = "Component 2"
)
mSet <- PLSDA.CV(mSet, cvOpt = "loo", compNum = min(5, length(groups) + 2), choice = "Q2")
mSet <- PlotPLS.Classification(mSet, "pls_cv_", "png", 72, width = NA)
mSet <- PlotPLS.Imp(mSet, "pls_vip_", "png", 72, width = NA, type = "vip", feat.nm = "Comp. 1", feat.num = 15, color.BW = FALSE)

if (length(groups) == 2) {
  mSet <- Ttests.Anal(mSet, nonpar = FALSE, threshp = 0.05, paired = FALSE, equal.var = TRUE)
  mSet <- PlotTT(mSet, "tt_", "png", 72, width = NA)
}

mSet <- ANOVA.Anal(mSet, nonpar = FALSE, thresh = 0.05, all_results = TRUE)
mSet <- PlotANOVA(mSet, "anova_", "png", 72, width = NA)

mSet <- tryCatch(
  PlotHeatMap(
    mSet,
    "heatmap_",
    "png",
    72,
    width = NA,
    dataOpt = "norm",
    scaleOpt = "row",
    smplDist = "euclidean",
    clstDist = "ward.D",
    palette = "bwm",
    fzCol = 8,
    fzRow = 8,
    fzAnno = 8,
    annoPer = 1,
    unitCol = 0.25,
    unitRow = 0.25,
    rowV = TRUE,
    colV = TRUE,
    var.inx = NULL,
    border = TRUE,
    grp.ave = FALSE
  ),
  error = function(e) {
    warning("Using local heatmap fallback: ", conditionMessage(e), call. = FALSE)
    plot_local_heatmap(
      mSet$dataSet$norm,
      mSet$dataSet$cls,
      file.path(output_dir, "heatmap_local.png")
    )
    mSet
  }
)

saveRDS(mSet, file.path(output_dir, "metaboanalyst_mSet.rds"))

cat("\n=== Downstream MetaboAnalystR analysis complete ===\n")
cat("MetaboAnalystR version:", as.character(utils::packageVersion("MetaboAnalystR")), "\n")
cat("Input file:", normalizePath(input_csv), "\n")
cat("Groups:\n")
print(table(input_data$Group))
cat("Output directory:", normalizePath(output_dir), "\n")
