# Run Path B: MetaboAnalystR mummichog/MS Peaks-to-Pathways.
#
# Run from the repository root after setting SPATIALOMICS_RAT_BRAIN_DIR.
source("scripts/_bootstrap.R")

library(MetaboAnalystR)

if (!exists("data_dir", inherits = TRUE)) {
  data_dir <- spatialomics_data_dir(
    "SPATIALOMICS_RAT_BRAIN_DIR", "data_raw/rat_brain_data", "Rat-brain data"
  )
}
if (!exists("ppm", inherits = TRUE)) ppm <- 5.0
if (!exists("ion_mode", inherits = TRUE)) ion_mode <- "negative"
if (!exists("mummichog_library", inherits = TRUE)) mummichog_library <- "rno_kegg"
if (!exists("mummichog_library_version", inherits = TRUE)) mummichog_library_version <- "current"
if (!exists("mummichog_p_cutoff", inherits = TRUE)) mummichog_p_cutoff <- 0.05
if (!exists("mummichog_min_lib", inherits = TRUE)) mummichog_min_lib <- 3
if (!exists("mummichog_perm_num", inherits = TRUE)) mummichog_perm_num <- 100

pathway_dir <- file.path(data_dir, "pathway_inputs")
input_csv <- file.path(pathway_dir, "path_b_mummichog_input.csv")
metabo_input_csv <- file.path(pathway_dir, "path_b_mummichog_metaboanalyst_input.csv")
out_dir <- file.path(pathway_dir, "path_b_metaboanalyst_mummichog_results_update_instrument")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(input_csv)) stop("Missing Path B mummichog input: ", input_csv, call. = FALSE)

path_b <- read.csv(input_csv, check.names = FALSE, stringsAsFactors = FALSE)
metabo_input <- data.frame(
  m.z = path_b[["m.z"]],
  p.value = path_b[["p.value"]],
  t.score = path_b[["statistic"]],
  check.names = FALSE
)
write.csv(metabo_input, metabo_input_csv, row.names = FALSE)

old_wd <- getwd()
setwd(out_dir)
on.exit(setwd(old_wd), add = TRUE)

current.msg <<- list()
err.vec <<- character()

mSet <- InitDataObjects("mass_all", "mummichog", FALSE, default.dpi = 72)
mSet <- SetPeakFormat(mSet, "mpt")

# Official MetaboAnalystR path for configuring ppm tolerance and negative mode.
# Directly setting mSet$dataSet$mode is not sufficient for adduct matching.
mSet <- UpdateInstrumentParameters(mSet, ppm, ion_mode)

mSet <- SetRTincluded(mSet, "no")
mSet <- Read.PeakListData(mSet, metabo_input_csv)
mSet <- SanityCheckMummichogData(mSet)
mSet <- SetPeakEnrichMethod(mSet, "mum", version = "v2")
mSet <- SetMummichogPval(mSet, mummichog_p_cutoff)
mSet <- PerformPSEA(
  mSet,
  mummichog_library,
  mummichog_library_version,
  minLib = mummichog_min_lib,
  permNum = mummichog_perm_num
)

saveRDS(mSet, "path_b_mummichog_mSet.rds")
if (is.list(mSet) && !is.null(mSet$mummi.resmat)) {
  pathway_results <- as.data.frame(mSet$mummi.resmat, stringsAsFactors = FALSE)
  pathway_results$Pathway <- rownames(pathway_results)
  pathway_results <- pathway_results[, c("Pathway", setdiff(names(pathway_results), "Pathway"))]
  write.csv(pathway_results, "mummichog_pathway_results_named.csv", row.names = FALSE)

  plot_ok <- tryCatch(
    {
      mSet <- PlotPeaks2Paths(
        mSet,
        "path_b_peaks_to_paths",
        "png",
        300,
        width = 9,
        labels = "default",
        num_annot = 5,
        interactive = FALSE
      )
      TRUE
    },
    error = function(e) {
      warning("Skipping MetaboAnalystR PlotPeaks2Paths: ", conditionMessage(e), call. = FALSE)
      FALSE
    }
  )
  if (isTRUE(plot_ok)) saveRDS(mSet, "path_b_mummichog_mSet.rds")
}
if (is.list(mSet) && !is.null(mSet$mummi)) {
  capture.output(str(mSet$mummi), file = "mummichog_mummi_structure.txt")
}

cat("\n=== Path B mummichog complete ===\n")
cat("Input:", normalizePath(metabo_input_csv), "\n")
cat("Library:", mummichog_library, "\n")
cat("Ion mode:", ion_mode, "\n")
cat("PPM:", ppm, "\n")
if (is.list(mSet) && !is.null(mSet$mummi.resmat)) {
  cat("Pathways:", nrow(mSet$mummi.resmat), "\n")
  print(head(mSet$mummi.resmat, 20))
}
cat("Output directory:", normalizePath(out_dir), "\n")
