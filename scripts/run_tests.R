#!/usr/bin/env Rscript

arguments <- commandArgs(trailingOnly = TRUE)
project_root <- normalizePath(if (length(arguments)) arguments[1] else getwd(), mustWork = TRUE)
if (!file.exists(file.path(project_root, "DESCRIPTION"))) {
  stop("Run this script from the repository root or pass the root as its first argument.", call. = FALSE)
}

test_files <- sort(list.files(file.path(project_root, "tests"), pattern = "\\.R$", full.names = TRUE))
if (!length(test_files)) stop("No tests were found.", call. = FALSE)

failures <- character()
for (test_file in test_files) {
  cat("\n=== ", basename(test_file), " ===\n", sep = "")
  status <- system2(
    file.path(R.home("bin"), "Rscript"),
    shQuote(test_file),
    stdout = "",
    stderr = "",
    env = paste0("SPATIALOMICS_PROJECT_ROOT=", shQuote(project_root))
  )
  if (!identical(status, 0L)) failures <- c(failures, basename(test_file))
}

if (length(failures)) {
  stop("Test failure(s): ", paste(failures, collapse = ", "), call. = FALSE)
}
cat("\nALL_TESTS_OK=TRUE\n")
