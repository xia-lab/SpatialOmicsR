checking_installed_package <- nzchar(Sys.getenv("_R_CHECK_PACKAGE_NAME_"))
if (dir.exists(".lib")) .libPaths(c(normalizePath(".lib"), .libPaths()))
if (checking_installed_package) library(SpatialOmicsMSI) else {
  source(file.path(normalizePath(getwd()), "R", "msi_pipeline.R"))
  source(file.path(normalizePath(getwd()), "R", "real_data_adapters.R"))
}

if (requireNamespace("base64enc", quietly = TRUE)) {
  # A minimal uncompressed, 64-bit mzML MS2 spectrum exercises metadata and arrays.
  encode <- function(x) base64enc::base64encode(writeBin(as.double(x), raw(), size = 8, endian = "little"))
  mz <- c(100, 200, 775.550137928655); intensity <- c(10, 100, 5)
  xml <- paste0('<mzML><run><spectrumList><spectrum id="scan=1">',
    '<cvParam accession="MS:1000511" value="2"/><cvParam accession="MS:1000127"/>',
    '<cvParam accession="MS:1000016" value="5.5"/><precursor>',
    '<cvParam accession="MS:1000744" value="775.550137928655"/>',
    '<cvParam accession="MS:1000827" value="775.55"/><cvParam accession="MS:1000828" value="0.5"/>',
    '<cvParam accession="MS:1000829" value="0.5"/><cvParam accession="MS:1000045" value="30"/></precursor>',
    '<binaryDataArray><cvParam accession="MS:1000514"/><cvParam accession="MS:1000523"/><binary>', encode(mz), '</binary></binaryDataArray>',
    '<binaryDataArray><cvParam accession="MS:1000515"/><cvParam accession="MS:1000523"/><binary>', encode(intensity), '</binary></binaryDataArray>',
    '</spectrum></spectrumList></run></mzML>')
  path <- tempfile(fileext = ".mzML"); writeLines(xml, path)
  out <- read_mzml_fragment_spectra(path, 775.550137928655, 10)
  stopifnot(nrow(out$precursor_scan_metadata) == 1L,
            nrow(out$fragment_peak_table) == 3L,
            out$precursor_scan_metadata$collision_energy == 30,
            out$precursor_scan_metadata$retention_time == 5.5,
            out$precursor_scan_metadata$spectrum_representation == "centroided",
            identical(out$parameters$chemical_identity_inferred, FALSE))
  cat("MZML_ADAPTER_TEST_OK=TRUE\n")
} else {
  cat("MZML_ADAPTER_TEST_SKIPPED_BASE64ENC=TRUE\n")
}
