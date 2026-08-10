# MetaboAnalystR PCA integration test
#
# Run from the repository root after setting SPATIALOMICS_RAT_BRAIN_DIR.
source("scripts/_bootstrap.R")

required_pkgs <- c("MetaboAnalystR", "RSclient", "factoextra")
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop(
    "Missing required package(s): ",
    paste(missing_pkgs, collapse = ", "),
    ". Install them before running this test.",
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

load_spatialomics_code()

data_dir <- spatialomics_data_dir(
  "SPATIALOMICS_RAT_BRAIN_DIR", "data_raw/rat_brain_data", "Rat-brain data"
)
out_dir <- file.path(data_dir, "spatial_test_outputs")
plot_dir <- file.path(out_dir, "plots")
ma_dir <- file.path(out_dir, "metaboanalystR_test")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(ma_dir, recursive = TRUE, showWarnings = FALSE)

clustered_path <- file.path(out_dir, "clustered_matrix.csv")
if (!file.exists(clustered_path)) {
  stop(
    "Missing clustered_matrix.csv. Run scripts/test_rat_brain_workflow.R first.",
    call. = FALSE
  )
}

clustered <- read.csv(clustered_path, check.names = FALSE)
all_cluster_rois <- define_rois(clustered, mode = "cluster")
sampled_all <- sample_subregions(all_cluster_rois, grid_size = 5, min_pixels = 30)
metabo_all <- make_metaboanalyst_data(sampled_all$sample_matrix)

metabo_input_path <- file.path(ma_dir, "metaboanalyst_all_clusters_data.csv")
sample_mapping_path <- file.path(ma_dir, "metaboanalyst_all_clusters_sample_mapping.csv")
write.csv(metabo_all, metabo_input_path, row.names = FALSE)
write.csv(sampled_all$sample_mapping, sample_mapping_path, row.names = FALSE)

old_wd <- getwd()
on.exit(setwd(old_wd), add = TRUE)
setwd(ma_dir)

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

# MetaboAnalystR uses these globals internally in a few message helpers.
current.msg <- ""

mSet <- InitDataObjects("conc", "stat", FALSE, default.dpi = 72)
mSet <- Read.TextData(mSet, metabo_input_path, "rowu", "disc")
mSet <- SanityCheckData(mSet)
mSet <- Normalization(mSet, "NULL", "NULL", "NULL", ratio = FALSE)
mSet <- run_local_pca(mSet)

pca_scores <- data.frame(
  Sample = rownames(mSet$analSet$pca$x),
  PC1 = mSet$analSet$pca$x[, 1],
  PC2 = if (ncol(mSet$analSet$pca$x) >= 2) mSet$analSet$pca$x[, 2] else 0,
  stringsAsFactors = FALSE
)
write.csv(pca_scores, file.path(ma_dir, "metaboanalystR_pca_scores.csv"), row.names = FALSE)

pca_backmap <- backmap_sample_scores(pca_scores, sampled_all$sample_mapping, "PC1")
write.csv(pca_backmap, file.path(ma_dir, "metaboanalystR_pca_pc1_backmap.csv"), row.names = FALSE)

p_score_map <- plot_sample_score_map(
  pca_scores,
  sampled_all$sample_mapping,
  clustered,
  score_column = "PC1"
)
ggplot2::ggsave(
  file.path(plot_dir, "metaboanalystR_pca_pc1_score_map.png"),
  p_score_map,
  width = 7,
  height = 9,
  dpi = 300
)

cat("\n=== MetaboAnalystR PCA test complete ===\n")
cat("MetaboAnalystR version:", as.character(utils::packageVersion("MetaboAnalystR")), "\n")
cat("Input samples:", nrow(metabo_all), "\n")
cat("Groups:\n")
print(table(metabo_all$Group))
cat("PCA scores:", nrow(pca_scores), "\n")
cat("Output directory:", normalizePath(ma_dir), "\n")
cat("Score map:", normalizePath(file.path(plot_dir, "metaboanalystR_pca_pc1_score_map.png")), "\n")
