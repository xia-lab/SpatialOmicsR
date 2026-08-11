run_spatial_app <- function(launch.browser = TRUE, port = NULL, host = "127.0.0.1") {
  if (!requireNamespace("shiny", quietly = TRUE)) {
    stop("Package 'shiny' is required to run the app.", call. = FALSE)
  }
  if (length(host) != 1L || is.na(host) || !nzchar(trimws(host))) {
    stop("host must be one non-empty address, such as '127.0.0.1' or '0.0.0.0'.", call. = FALSE)
  }
  if (!is.null(port) && (length(port) != 1L || !is.finite(port) || port < 1 || port > 65535 || port != as.integer(port))) {
    stop("port must be NULL or one integer from 1 to 65535.", call. = FALSE)
  }

  app_dir <- system.file("shiny", "spatial_pipeline", package = "SpatialOmicsMSI")
  if (!nzchar(app_dir)) {
    stop("Cannot find the bundled Shiny app. Reinstall SpatialOmicsMSI.", call. = FALSE)
  }

  shiny::runApp(app_dir, launch.browser = launch.browser, port = port, host = host)
}
