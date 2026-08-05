# SpatialOmicsR
R functions to convert common outputs from spatial omics tools into 2D omics data and metadata tables, ready for statistical and functional analysis

## Shiny spatial pipeline

The package includes a Shiny demo app at `inst/shiny/spatial_pipeline/app.R`.

- Use the app to load MSI data, select features, preprocess, cluster, define ROIs, and export MetaboAnalyst input.
- Manual polygon ROI input supports both MSI coordinate CSV and histology (H&E) coordinate CSV.
- For histology coordinate mode, upload an H&E control points CSV with columns:
  - `histology_x`, `histology_y`, `msi_x`, `msi_y`
  - optional `section_id` for serial section workflows.
- Polygon vertex CSV should contain:
  - `roi_id`, `vertex_order`, `x`, `y`
  - optional `section_id` for serial sections.

Run the app with:

```r
shiny::runApp("inst/shiny/spatial_pipeline")
```
