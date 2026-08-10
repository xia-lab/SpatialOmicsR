source("scripts/_bootstrap.R")

shiny::runApp(
  spatialomics_path("inst", "shiny", "spatial_pipeline"),
  host = "127.0.0.1",
  port = 3838,
  launch.browser = TRUE
)
