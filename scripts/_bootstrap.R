# Shared setup for repository analysis scripts.

find_spatialomics_root <- function(start = getwd()) {
  configured <- Sys.getenv("SPATIALOMICS_PROJECT_ROOT", unset = "")
  candidates <- unique(c(configured, start, dirname(start)))
  candidates <- candidates[nzchar(candidates)]
  for (candidate in candidates) {
    candidate <- normalizePath(candidate, mustWork = FALSE)
    if (file.exists(file.path(candidate, "DESCRIPTION")) &&
        dir.exists(file.path(candidate, "R"))) {
      return(candidate)
    }
  }
  stop(
    "Cannot locate the SpatialOmicsMSI repository. Run from the repository root ",
    "or set SPATIALOMICS_PROJECT_ROOT.",
    call. = FALSE
  )
}

spatialomics_project_root <- find_spatialomics_root()

spatialomics_path <- function(...) {
  file.path(spatialomics_project_root, ...)
}

spatialomics_data_dir <- function(environment_variable,
                                  relative_default,
                                  label = environment_variable,
                                  must_exist = TRUE) {
  configured <- Sys.getenv(environment_variable, unset = "")
  path <- if (nzchar(configured)) configured else spatialomics_path(relative_default)
  path <- normalizePath(path, mustWork = FALSE)
  if (isTRUE(must_exist) && !dir.exists(path)) {
    stop(
      label, " directory does not exist: ", path, "\n",
      "Set ", environment_variable, " to the local data directory. Raw data are not stored in Git.",
      call. = FALSE
    )
  }
  path
}

load_spatialomics_code <- function(prefer_installed = FALSE) {
  if (isTRUE(prefer_installed) && requireNamespace("SpatialOmicsMSI", quietly = TRUE)) {
    suppressPackageStartupMessages(library(SpatialOmicsMSI))
  } else {
    files <- list.files(spatialomics_path("R"), pattern = "\\.R$", full.names = TRUE)
    package_file <- files[basename(files) == "package.R"]
    files <- c(setdiff(files, package_file), package_file)
    invisible(lapply(files, sys.source, envir = .GlobalEnv))
  }
  invisible(TRUE)
}
