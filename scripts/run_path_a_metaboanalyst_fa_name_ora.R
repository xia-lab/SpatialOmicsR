# Run Path A: standard MetaboAnalystR KEGG ORA using FA common names.
#
# The SwissLipids CHEBI IDs for these lipid species are not present in the
# MetaboAnalystR compound_db used by CrossReferencing("chebi"). The official,
# well-tested mapping path is name-based matching, so the eight free fatty acid
# species are mapped to common names before KEGG ORA.
#
# Run from the repository root after setting SPATIALOMICS_RAT_BRAIN_DIR.
source("scripts/_bootstrap.R")

library(MetaboAnalystR)

if (!exists("data_dir", inherits = TRUE)) {
  data_dir <- spatialomics_data_dir(
    "SPATIALOMICS_RAT_BRAIN_DIR", "data_raw/rat_brain_data", "Rat-brain data"
  )
}
if (!exists("kegg_organism", inherits = TRUE)) kegg_organism <- "rno"
if (!exists("kegg_library_version", inherits = TRUE)) kegg_library_version <- "current"
if (!exists("node_importance", inherits = TRUE)) node_importance <- "rbc"
if (!exists("ora_method", inherits = TRUE)) ora_method <- "hyperg"

pathway_dir <- file.path(data_dir, "pathway_inputs")
out_dir <- file.path(pathway_dir, "path_a_metaboanalyst_fa_name_results")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

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

old_wd <- getwd()
setwd(out_dir)
on.exit(setwd(old_wd), add = TRUE)

write.csv(fa_names, "path_a_fa_common_name_input.csv", row.names = FALSE)

current.msg <<- list()
err.vec <<- character()

mSet <- InitDataObjects("conc", "pathora", FALSE, default.dpi = 72)
mSet <- Setup.MapData(mSet, fa_names$common_name)
mSet <- CrossReferencing(
  mSet,
  "name",
  hmdb = TRUE,
  pubchem = TRUE,
  chebi = TRUE,
  kegg = TRUE,
  metlin = FALSE
)

mapping_table <- as.data.frame(mSet$name.map, stringsAsFactors = FALSE)
write.csv(mapping_table, "path_a_name_mapping.csv", row.names = FALSE)

mapped <- !is.na(mSet$name.map$hit.inx)
if (sum(mapped) < 3) {
  stop("Fewer than 3 compounds mapped; KEGG ORA would not be meaningful.", call. = FALSE)
}

mSet <- SetKEGG.PathLib(mSet, kegg_organism, kegg_library_version)

# MetaboAnalystR 4.3.0 local KEGG API branch expects this flag to exist.
mSet$api$filter <- FALSE

mSet <- CalculateOraScore(mSet, node_importance, ora_method)
saveRDS(mSet, "path_a_fa_name_pathora_mSet.rds")

if (is.list(mSet) && !is.null(mSet$analSet$ora.mat)) {
  ora_results <- as.data.frame(mSet$analSet$ora.mat, stringsAsFactors = FALSE)
  ora_results$pathway_id <- rownames(ora_results)
  ora_results <- ora_results[, c("pathway_id", setdiff(names(ora_results), "pathway_id"))]
  write.csv(ora_results, "path_a_fa_name_kegg_ora_results.csv", row.names = FALSE)

  # MetaboAnalystR 4.3.0 does not load current.kegglib in local API mode,
  # but PlotPathSummary() still needs it to translate KEGG IDs to names.
  current.kegglib <<- get(".get.my.lib", asNamespace("MetaboAnalystR"))(
    paste0(kegg_organism, ".qs"),
    "kegg/metpa"
  )
  mSet <- PlotPathSummary(
    mSet,
    show.grid = FALSE,
    imgName = "path_a_metaboanalyst_path_summary",
    format = "png",
    dpi = 300,
    width = 8
  )
  saveRDS(mSet, "path_a_fa_name_pathora_mSet.rds")
}

cat("\n=== Path A FA-name KEGG ORA complete ===\n")
cat("Mapped compounds:", sum(mapped), "of", nrow(fa_names), "\n")
cat("KEGG organism:", kegg_organism, "\n")
if (is.list(mSet) && !is.null(mSet$analSet$ora.mat)) {
  cat("Pathways:", nrow(mSet$analSet$ora.mat), "\n")
  print(head(mSet$analSet$ora.mat, 20))
}
cat("Output directory:", normalizePath(out_dir), "\n")
