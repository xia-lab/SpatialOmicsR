# Run every package regression test independently and retain one log per
# functional test area. This runner deliberately keeps the repository root as
# the working directory so opt-in real-data fixtures can be discovered.

project_root <- normalizePath(getwd(), mustWork = TRUE)
output_dir <- file.path(project_root, "outputs", "all_feature_validation")
log_dir <- file.path(output_dir, "logs")
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

test_files <- sort(list.files(
  file.path(project_root, "tests", "testthat"),
  pattern = "^test-.*\\.R$", full.names = TRUE
))
local_library <- file.path(project_root, ".lib")
rscript <- file.path(R.home("bin"), "Rscript")

run_one <- function(test_file) {
  test_name <- tools::file_path_sans_ext(basename(test_file))
  log_file <- file.path(log_dir, paste0(test_name, ".log"))
  expression <- sprintf(
    "setwd(%s); sys.source(%s, envir = new.env(parent = globalenv()))",
    encodeString(project_root, quote = '"'), encodeString(test_file, quote = '"')
  )
  started <- Sys.time()
  status <- system2(
    rscript, c("-e", shQuote(expression)),
    stdout = log_file, stderr = log_file,
    env = c(
      paste0("R_LIBS_USER=", local_library),
      "SPATIALOMICS_RUN_REAL_DATA_TESTS=true",
      paste0("SPATIALOMICS_PROJECT_ROOT=", project_root)
    )
  )
  data.frame(
    test_file = basename(test_file),
    status = if (identical(status, 0L)) "passed" else "failed",
    exit_code = as.integer(status),
    elapsed_seconds = round(as.numeric(difftime(Sys.time(), started, units = "secs")), 3),
    log_file = file.path("logs", basename(log_file)),
    stringsAsFactors = FALSE
  )
}

test_manifest <- do.call(rbind, lapply(test_files, run_one))
utils::write.csv(test_manifest, file.path(output_dir, "test_manifest.csv"), row.names = FALSE)
saveRDS(test_manifest, file.path(output_dir, "test_manifest.rds"), compress = "xz")

namespace_lines <- readLines(file.path(project_root, "NAMESPACE"), warn = FALSE)
exports <- sub("^export\\(([^)]+)\\).*$", "\\1", grep("^export\\(", namespace_lines, value = TRUE))
test_text <- lapply(test_files, function(path) paste(readLines(path, warn = FALSE), collapse = "\n"))
coverage <- do.call(rbind, lapply(exports, function(fun) {
  referenced <- vapply(test_text, function(text) grepl(fun, text, fixed = TRUE), logical(1))
  data.frame(
    exported_function = fun,
    directly_referenced = any(referenced),
    test_files = paste(basename(test_files[referenced]), collapse = ";"),
    stringsAsFactors = FALSE
  )
}))
utils::write.csv(coverage, file.path(output_dir, "exported_function_coverage.csv"), row.names = FALSE)

summary <- data.frame(
  metric = c("test_files", "passed", "failed", "exported_functions", "directly_referenced_exports"),
  value = c(
    nrow(test_manifest), sum(test_manifest$status == "passed"), sum(test_manifest$status == "failed"),
    nrow(coverage), sum(coverage$directly_referenced)
  ),
  stringsAsFactors = FALSE
)
utils::write.csv(summary, file.path(output_dir, "validation_summary.csv"), row.names = FALSE)
print(test_manifest)
print(summary)
if (any(test_manifest$status == "failed")) quit(status = 1L)
