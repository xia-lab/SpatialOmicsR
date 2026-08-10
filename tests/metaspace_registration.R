checking_installed_package <- nzchar(Sys.getenv("_R_CHECK_PACKAGE_NAME_"))
if (dir.exists(".lib")) .libPaths(c(normalizePath(".lib"), .libPaths()))
if (checking_installed_package) {
  library(SpatialOmicsMSI)
  registration_mask_diagnostics <- getFromNamespace("registration_mask_diagnostics", "SpatialOmicsMSI")
} else {
  source(file.path(normalizePath(getwd()), "R", "msi_pipeline.R"))
  source(file.path(normalizePath(getwd()), "R", "real_data_adapters.R"))
}

coordinates <- data.frame(pixel_id = 1:3, x = c(82, 83, 82), y = c(340, 340, 341))
matrix <- rbind(c(4, 0, 100), c(0, 4, 200), c(0, 0, 1))
mapped <- apply_metaspace_transform(coordinates, matrix)
stopifnot(identical(mapped$relative_x, c(0, 1, 0)),
          identical(mapped$relative_y, c(0, 0, 1)),
          identical(mapped$optical_x, c(100, 104, 100)),
          identical(mapped$optical_y, c(200, 200, 204)))

if (requireNamespace("jsonlite", quietly = TRUE)) {
  json <- tempfile(fileext = ".json")
  writeLines('{"data":{"rawOpticalImage":{"url":"https://example.invalid/a.jpg","transform":[[4,0,100],[0,4,200],[0,0,1]]}}}', json)
  parsed <- read_metaspace_transform(json)
  stopifnot(identical(parsed$transform, matrix),
            identical(parsed$direction, "relative_MSI_coordinates_to_optical_pixels"))
}

measurement <- matrix(FALSE, 8, 8); measurement[3:6,3:6] <- TRUE
tissue <- measurement
mask_stats <- registration_mask_diagnostics(measurement, tissue)
stopifnot(mask_stats[["overlap"]] == 1, mask_stats[["dice"]] == 1,
          mask_stats[["boundary_distance_median_px"]] == 0)
tissue[3,3] <- FALSE
changed <- registration_mask_diagnostics(measurement, tissue)
stopifnot(changed[["overlap"]] < 1, changed[["dice"]] < 1)

run_real_data_tests <- identical(
  tolower(Sys.getenv("SPATIALOMICS_RUN_REAL_DATA_TESTS", unset = "false")), "true"
)
if (!checking_installed_package && run_real_data_tests) {
  root<-file.path(getwd(),"data_raw","mouse_brain_he_msi","metaspace_brain01")
  rds<-"/tmp/brain01_final_adapter.rds"
  if(file.exists(rds)&&file.exists(file.path(root,"optical_brightfield.jpg"))){
    real<-readRDS(rds);reg<-register_metaspace_optical(real$coordinates,
      file.path(root,"optical_transform_api.json"),file.path(root,"optical_brightfield.jpg"),
      file.path(root,"attribution_license.json"))
    key_row <- which(real$coordinates$x == 154 & real$coordinates$y == 340)[1]
    stopifnot(!is.na(key_row))
    stopifnot(all.equal(reg$coordinate_origin,c(82,340)),
      abs(reg$registered_coordinates$optical_x[key_row]-668.639294041797)<1e-8,
      abs(reg$registered_coordinates$optical_y[key_row]-673.515633162620)<1e-8,
      abs(reg$diagnostics$overlap-0.9966402223)<1e-5,
      abs(reg$diagnostics$dice-0.9770005573)<2e-5,
      abs(reg$diagnostics$boundary_distance_median_px-4.123105626)<1e-5,
      abs(reg$diagnostics$boundary_distance_p95_px-15.524174696)<1e-5)
  }
}
cat("METASPACE_REGISTRATION_TEST_OK=TRUE\n")
