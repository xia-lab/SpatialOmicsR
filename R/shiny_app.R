run_spatial_app <- function(launch.browser = TRUE, port = NULL) {
  if (!requireNamespace("shiny", quietly = TRUE)) {
    stop("Package 'shiny' is required to run the app.", call. = FALSE)
  }

  app_dir <- system.file("shiny", "spatial_pipeline", package = "SpatialOmicsMSI")
  if (!nzchar(app_dir)) {
    stop("Cannot find the bundled Shiny app. Reinstall SpatialOmicsMSI.", call. = FALSE)
  }

  shiny::runApp(app_dir, launch.browser = launch.browser, port = port)
}
