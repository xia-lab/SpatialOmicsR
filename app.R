host <- Sys.getenv("SPATIALOMICS_SHINY_HOST", unset = "127.0.0.1")
port <- suppressWarnings(as.integer(Sys.getenv("SPATIALOMICS_SHINY_PORT", unset = "3838")))
if (is.na(port) || port < 1L || port > 65535L) {
  stop("SPATIALOMICS_SHINY_PORT must be an integer from 1 to 65535.", call. = FALSE)
}
shiny::runApp("inst/shiny/spatial_pipeline", host = host, port = port, launch.browser = FALSE)
