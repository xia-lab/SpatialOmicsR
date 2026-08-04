# Path A contamination sensitivity analysis.
#
# R console usage:
# setwd("/Users/ly/Documents/Spatial Omics")
# source("scripts/run_path_a_contamination_sensitivity.R")

source("R/msi_pipeline.R")
library(MetaboAnalystR)

if (!exists("data_dir", inherits = TRUE)) {
  data_dir <- "/Users/ly/Desktop/Jeff Xia/rat_brain_data"
}
if (!exists("kegg_organism", inherits = TRUE)) kegg_organism <- "rno"
if (!exists("kegg_library_version", inherits = TRUE)) kegg_library_version <- "current"
if (!exists("node_importance", inherits = TRUE)) node_importance <- "rbc"
if (!exists("ora_method", inherits = TRUE)) ora_method <- "hyperg"

test_out_dir <- file.path(data_dir, "spatial_test_outputs")
pathway_dir <- file.path(data_dir, "pathway_inputs")
out_dir <- file.path(pathway_dir, "path_a_contamination_sensitivity")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

pixel_matrix_csv <- file.path(test_out_dir, "preprocessed_matrix.csv")
path_a_feature_csv <- file.path(pathway_dir, "path_a_feature_chebi_stats.csv")

for (path in c(pixel_matrix_csv, path_a_feature_csv)) {
  if (!file.exists(path)) stop("Missing required file: ", path, call. = FALSE)
}

fa_names <- data.frame(
  swiss_name = c(
    "hexadecenoate",
    "hexadecanoate",
    "octadecadienoate",
    "octadecenoate",
    "octadecanoate",
    "eicosatetraenoate",
    "docosahexaenoate",
    "docosatetraenoate"
  ),
  common_name = c(
    "Palmitoleic acid",
    "Palmitic acid",
    "Linoleic acid",
    "Oleic acid",
    "Stearic acid",
    "Arachidonic acid",
    "Docosahexaenoic acid",
    "Adrenic acid"
  ),
  stringsAsFactors = FALSE
)

run_name_ora <- function(input_names, output_prefix) {
  old_wd <- getwd()
  setwd(out_dir)
  on.exit(setwd(old_wd), add = TRUE)

  current.msg <<- list()
  err.vec <<- character()

  mSet <- InitDataObjects("conc", "pathora", FALSE, default.dpi = 72)
  mSet <- Setup.MapData(mSet, input_names)
  mSet <- CrossReferencing(
    mSet,
    "name",
    hmdb = TRUE,
    pubchem = TRUE,
    chebi = TRUE,
    kegg = TRUE,
    metlin = FALSE
  )

  mapped <- !is.na(mSet$name.map$hit.inx)
  mapping_table <- as.data.frame(mSet$name.map, stringsAsFactors = FALSE)
  write.csv(mapping_table, paste0(output_prefix, "_name_mapping.csv"), row.names = FALSE)

  if (sum(mapped) < 3) {
    warning(output_prefix, ": fewer than 3 mapped compounds; skipping ORA.", call. = FALSE)
    return(list(mSet = mSet, mapped = sum(mapped), ora = data.frame()))
  }

  mSet <- SetKEGG.PathLib(mSet, kegg_organism, kegg_library_version)
  mSet$api$filter <- FALSE
  mSet <- CalculateOraScore(mSet, node_importance, ora_method)
  saveRDS(mSet, paste0(output_prefix, "_pathora_mSet.rds"))

  ora <- as.data.frame(mSet$analSet$ora.mat, stringsAsFactors = FALSE)
  ora$pathway_id <- rownames(ora)
  ora <- ora[, c("pathway_id", setdiff(names(ora), "pathway_id"))]
  write.csv(ora, paste0(output_prefix, "_kegg_ora_results.csv"), row.names = FALSE)

  current.kegglib <<- get(".get.my.lib", asNamespace("MetaboAnalystR"))(
    paste0(kegg_organism, ".qs"),
    "kegg/metpa"
  )
  mSet <- PlotPathSummary(
    mSet,
    show.grid = FALSE,
    imgName = paste0(output_prefix, "_metaboanalyst_path_summary"),
    format = "png",
    dpi = 300,
    width = 8
  )
  saveRDS(mSet, paste0(output_prefix, "_pathora_mSet.rds"))

  list(mSet = mSet, mapped = sum(mapped), ora = ora)
}

pixel_matrix <- read.csv(pixel_matrix_csv, check.names = FALSE, stringsAsFactors = FALSE)
path_a_features <- read.csv(path_a_feature_csv, check.names = FALSE, stringsAsFactors = FALSE)

contam_diag <- compare_contamination_methods(
  pixel_matrix,
  include_moran = FALSE,
  agreement_methods = c("gap", "mad", "enrichment")
)
contam_summary <- contam_diag$agreement
write.csv(contam_summary, file.path(out_dir, "contamination_diagnostic_agreement.csv"), row.names = FALSE)

known_contaminants <- c("mz_255.23298", "mz_283.26430", "mz_554.26209")
strict_contaminants <- contam_summary$column_name[contam_summary$contaminant_agreement >= 3]

path_a_feature_names <- unique(path_a_features[, c("feature", "swiss_name_master"), drop = FALSE])
path_a_feature_names <- merge(path_a_feature_names, fa_names, by.x = "swiss_name_master", by.y = "swiss_name", all.x = TRUE, sort = FALSE)

scenario_table <- rbind(
  data.frame(scenario = "original", contaminant_rule = "none", path_a_feature_names, excluded = FALSE, stringsAsFactors = FALSE),
  data.frame(
    scenario = "known_contaminants_removed",
    contaminant_rule = paste(known_contaminants, collapse = ";"),
    path_a_feature_names,
    excluded = path_a_feature_names$feature %in% known_contaminants,
    stringsAsFactors = FALSE
  ),
  data.frame(
    scenario = "strict_agreement3_removed",
    contaminant_rule = paste(strict_contaminants, collapse = ";"),
    path_a_feature_names,
    excluded = path_a_feature_names$feature %in% strict_contaminants,
    stringsAsFactors = FALSE
  )
)
write.csv(scenario_table, file.path(out_dir, "path_a_feature_contamination_scenarios.csv"), row.names = FALSE)

scenarios <- split(scenario_table, scenario_table$scenario)
results <- lapply(names(scenarios), function(scenario) {
  rows <- scenarios[[scenario]]
  kept <- rows[!rows$excluded & !is.na(rows$common_name), , drop = FALSE]
  prefix <- paste0("path_a_", scenario)
  write.csv(kept, file.path(out_dir, paste0(prefix, "_common_name_input.csv")), row.names = FALSE)
  ora_result <- run_name_ora(unique(kept$common_name), prefix)
  rno01040 <- ora_result$ora[ora_result$ora$pathway_id == "rno01040", , drop = FALSE]
  data.frame(
    scenario = scenario,
    input_compounds = length(unique(kept$common_name)),
    mapped_compounds = ora_result$mapped,
    rno01040_hits = if (nrow(rno01040) > 0) rno01040$Hits else NA,
    rno01040_raw_p = if (nrow(rno01040) > 0) rno01040[["Raw p"]] else NA,
    rno01040_fdr = if (nrow(rno01040) > 0) rno01040$FDR else NA,
    stringsAsFactors = FALSE
  )
})
sensitivity_summary <- do.call(rbind, results)
write.csv(sensitivity_summary, file.path(out_dir, "path_a_contamination_sensitivity_summary.csv"), row.names = FALSE)

cat("\n=== Path A contamination sensitivity complete ===\n")
print(sensitivity_summary)
cat("Known contaminants:", paste(known_contaminants, collapse = ", "), "\n")
cat("Strict agreement>=3 contaminants:", paste(strict_contaminants, collapse = ", "), "\n")
cat("Output directory:", normalizePath(out_dir), "\n")
