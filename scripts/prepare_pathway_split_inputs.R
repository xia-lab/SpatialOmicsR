# Prepare mutually exclusive Path A / Path B pathway-analysis inputs.
#
# Path A: high-confidence SwissLipids annotations with CHEBI IDs.
# Path B: all remaining features for mummichog/PSEA using m/z + statistics.
#
# Run from the repository root after setting the data environment variables.
source("scripts/_bootstrap.R")

split_ids <- function(x) {
  x <- paste(na.omit(x), collapse = "; ")
  ids <- unique(unlist(strsplit(x, "; ", fixed = TRUE)))
  ids[ids != ""]
}

best_signed_stat <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0) return(NA_real_)
  x[which.max(abs(x))]
}

if (!exists("data_dir", inherits = TRUE)) {
  data_dir <- spatialomics_data_dir(
    "SPATIALOMICS_RAT_BRAIN_DIR", "data_raw/rat_brain_data", "Rat-brain data"
  )
}
if (!exists("swisslipids_master_tsv", inherits = TRUE)) {
  swisslipids_master_tsv <- Sys.getenv(
    "SPATIALOMICS_SWISSLIPIDS_TSV",
    unset = file.path(data_dir, "swisslipids_master.tsv")
  )
}
if (!exists("pathway_out_dir", inherits = TRUE)) {
  pathway_out_dir <- file.path(data_dir, "pathway_inputs")
}
dir.create(pathway_out_dir, recursive = TRUE, showWarnings = FALSE)

annotation_csv <- file.path(data_dir, "annotation_cross_validation", "annotation_cross_validation.csv")
anova_csv <- file.path(pathway_out_dir, "all_feature_anova_effects.csv")
mummichog_all_csv <- file.path(pathway_out_dir, "mummichog_input_all_features.csv")

for (path in c(annotation_csv, anova_csv, mummichog_all_csv, swisslipids_master_tsv)) {
  if (!file.exists(path)) stop("Missing required file: ", path, call. = FALSE)
}

annotation <- read.csv(annotation_csv, check.names = FALSE, stringsAsFactors = FALSE)
anova <- read.csv(anova_csv, check.names = FALSE, stringsAsFactors = FALSE)
mummichog_all <- read.csv(mummichog_all_csv, check.names = FALSE, stringsAsFactors = FALSE)
lipids <- read.delim(swisslipids_master_tsv, check.names = FALSE, stringsAsFactors = FALSE)

high_conf <- annotation[annotation$consistency_class %in% c("exact_match", "candidate_overlap"), , drop = FALSE]
high_conf_ids <- split_ids(high_conf$swiss_ids)

lipid_xref <- lipids[lipids[["Lipid ID"]] %in% high_conf_ids, c(
  "Lipid ID", "Name", "Abbreviation*", "Level", "CHEBI", "LIPID MAPS", "HMDB"
), drop = FALSE]
names(lipid_xref) <- c("swiss_id", "swiss_name_master", "abbreviation", "level", "CHEBI", "LIPID_MAPS", "HMDB")
lipid_xref$has_CHEBI <- !is.na(lipid_xref$CHEBI) & lipid_xref$CHEBI != ""
lipid_xref$has_LIPID_MAPS <- !is.na(lipid_xref$LIPID_MAPS) & lipid_xref$LIPID_MAPS != ""
lipid_xref$has_HMDB <- !is.na(lipid_xref$HMDB) & lipid_xref$HMDB != ""

chebi_ids <- lipid_xref$swiss_id[lipid_xref$has_CHEBI]

feature_chebi_rows <- lapply(seq_len(nrow(high_conf)), function(i) {
  ids <- split_ids(high_conf$swiss_ids[i])
  ids <- intersect(ids, chebi_ids)
  if (length(ids) == 0) return(NULL)
  rows <- lipid_xref[match(ids, lipid_xref$swiss_id), , drop = FALSE]
  data.frame(
    mzmine_id = high_conf$feature_id[i],
    measured_mz = high_conf$measured_mz[i],
    consistency_class = high_conf$consistency_class[i],
    rows,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
})
feature_chebi_rows <- feature_chebi_rows[!vapply(feature_chebi_rows, is.null, logical(1))]
feature_chebi_map <- if (length(feature_chebi_rows) > 0) do.call(rbind, feature_chebi_rows) else data.frame()

path_a_feature_ids <- unique(feature_chebi_map$mzmine_id)
path_a_features <- anova[anova$mzmine_id %in% path_a_feature_ids, , drop = FALSE]
path_b_features <- anova[!anova$mzmine_id %in% path_a_feature_ids, , drop = FALSE]
path_b_mummichog <- mummichog_all[!mummichog_all$mzmine_id %in% path_a_feature_ids, , drop = FALSE]

path_a_feature_chebi_stats <- merge(
  feature_chebi_map,
  anova[, c("mzmine_id", "feature", "p.value", "FDR", "signed_log2fc", "F.stat", "max_group", "min_group"), drop = FALSE],
  by = "mzmine_id",
  all.x = TRUE,
  sort = FALSE
)

path_a_compound <- aggregate(
  cbind(p.value, FDR) ~ CHEBI,
  data = path_a_feature_chebi_stats,
  FUN = function(x) min(x, na.rm = TRUE)
)
path_a_effect <- aggregate(
  signed_log2fc ~ CHEBI,
  data = path_a_feature_chebi_stats,
  FUN = best_signed_stat
)
path_a_names <- aggregate(
  cbind(swiss_id, swiss_name_master, abbreviation) ~ CHEBI,
  data = path_a_feature_chebi_stats,
  FUN = function(x) paste(unique(x[!is.na(x) & x != ""]), collapse = "; ")
)
path_a_compound <- merge(path_a_compound, path_a_effect, by = "CHEBI", all.x = TRUE, sort = FALSE)
path_a_compound <- merge(path_a_compound, path_a_names, by = "CHEBI", all.x = TRUE, sort = FALSE)
path_a_compound$metaboanalyst_id <- paste0("CHEBI:", path_a_compound$CHEBI)

path_a_compound <- path_a_compound[, c(
  "metaboanalyst_id", "CHEBI", "swiss_id", "swiss_name_master",
  "abbreviation", "p.value", "FDR", "signed_log2fc"
), drop = FALSE]

write.csv(lipid_xref, file.path(pathway_out_dir, "pathway_high_conf_swisslipids_xrefs.csv"), row.names = FALSE)
write.csv(feature_chebi_map, file.path(pathway_out_dir, "path_a_feature_to_chebi_map.csv"), row.names = FALSE)
write.csv(path_a_feature_chebi_stats, file.path(pathway_out_dir, "path_a_feature_chebi_stats.csv"), row.names = FALSE)
write.csv(path_a_compound, file.path(pathway_out_dir, "path_a_metaboanalyst_chebi_input.csv"), row.names = FALSE)
write.csv(path_a_features, file.path(pathway_out_dir, "path_a_features_excluded_from_mummichog.csv"), row.names = FALSE)
write.csv(path_b_features, file.path(pathway_out_dir, "path_b_feature_stats.csv"), row.names = FALSE)
write.csv(path_b_mummichog, file.path(pathway_out_dir, "path_b_mummichog_input.csv"), row.names = FALSE)
utils::write.table(
  path_b_mummichog,
  file.path(pathway_out_dir, "path_b_mummichog_input.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

summary_table <- data.frame(
  metric = c(
    "total_features",
    "high_confidence_features",
    "high_confidence_swisslipids_ids",
    "high_confidence_ids_with_chebi",
    "high_confidence_ids_with_lipid_maps",
    "high_confidence_ids_with_hmdb",
    "path_a_features",
    "path_a_unique_chebi",
    "path_b_features"
  ),
  value = c(
    nrow(anova),
    nrow(high_conf),
    length(high_conf_ids),
    sum(lipid_xref$has_CHEBI),
    sum(lipid_xref$has_LIPID_MAPS),
    sum(lipid_xref$has_HMDB),
    length(path_a_feature_ids),
    length(unique(path_a_compound$CHEBI)),
    nrow(path_b_features)
  )
)
write.csv(summary_table, file.path(pathway_out_dir, "pathway_split_summary.csv"), row.names = FALSE)

cat("\n=== Pathway split input preparation complete ===\n")
print(summary_table)
cat("Path A CHEBI input:", normalizePath(file.path(pathway_out_dir, "path_a_metaboanalyst_chebi_input.csv")), "\n")
cat("Path B mummichog TSV:", normalizePath(file.path(pathway_out_dir, "path_b_mummichog_input.tsv")), "\n")
