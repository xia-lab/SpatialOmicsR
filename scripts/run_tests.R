#!/usr/bin/env Rscript

arguments <- commandArgs(trailingOnly = TRUE)
project_root <- normalizePath(if (length(arguments)) arguments[1] else getwd(), mustWork = TRUE)
if (!file.exists(file.path(project_root, "DESCRIPTION"))) {
  stop("Run this script from the repository root or pass the root as its first argument.", call. = FALSE)
}

if (!requireNamespace("testthat", quietly = TRUE)) {
  stop("Package 'testthat' (>= 3.0.0) is required. Install it before running tests.", call. = FALSE)
}
Sys.setenv(SPATIALOMICS_PROJECT_ROOT = project_root)
temporary_library <- tempfile("spatialomics-test-library-")
dir.create(temporary_library, recursive = TRUE)
on.exit(unlink(temporary_library, recursive = TRUE, force = TRUE), add = TRUE)
install_status <- system2(
  file.path(R.home("bin"), "R"),
  c("CMD", "INSTALL", "--no-multiarch", "--with-keep.source",
    paste0("--library=", shQuote(temporary_library)), shQuote(project_root)),
  stdout = "", stderr = ""
)
if (!identical(install_status, 0L)) {
  stop("Temporary installation failed; tests were not run.", call. = FALSE)
}
.libPaths(c(temporary_library, .libPaths()))
Sys.setenv(`_R_CHECK_PACKAGE_NAME_` = "SpatialOmicsMSI")
testthat::test_dir(file.path(project_root, "tests", "testthat"), reporter = "summary")
cat("\nALL_TESTS_OK=TRUE\n")
