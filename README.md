# SpatialOmicsMSI

`SpatialOmicsMSI` is an R package for reproducible spatial mass spectrometry
imaging (MSI) analysis. It standardizes processed CSV and imzML inputs, performs
QC and spatial analysis, connects MSI with optical images and cell-type labels,
and exports auditable results for downstream statistics and pathway analysis.

The package deliberately distinguishes exploratory pixel-level patterns from
population-level biological inference. Spatial clusters are not automatically
called anatomical regions, pixels are not treated as biological replicates, and
accurate-mass matches are not reported as confirmed molecular structures.

## Main capabilities

- Peak-picked CSV, shared-axis imzML, variable-axis imzML and mzML MS2 adapters.
- Background/tissue filtering, TIC normalization, `log10(x + 1)` and feature QC.
- Rook, Queen and distance-threshold spatial graphs.
- Moran's I, Geary's C and binSpect-inspired spatial variability analysis.
- PCA/k-means, BANKSY-inspired H0 clustering and ICM/Potts HMRF clustering.
- Fixed-window and connected-domain automatic ROI selection, target-enrichment
  scoring, multipart/hole-aware manual polygons and QuPath GeoJSON import.
- Histology-to-MSI affine registration and METASPACE optical transforms.
- Cell-label preparation from external labels, NNLS/cosine transfer or
  exploratory watershed morphology classes.
- Schürch-style neighborhood composition niches and permutation-based
  cell-type proximity enrichment.
- Expression-only and spatial ligand-receptor communication-potential analyses.
- Signed distance-to-domain rings and descriptive or subject-level GAMs.
- Subject-level t-test/Wilcoxon analysis and Gaussian spatial mixed models.
- Auditable all-candidate and global-optimal MSI/LC-MS matching, MS/MS and
  explicitly scoped CCS evidence integration.
- Vendor-neutral LCM ROI boundaries, registered target coordinates and
  ROI-level MSI--LCM quantitative concordance.
- MetaboAnalyst export/import, ORA/mummichog scripts and spatial back-mapping.
- A bundled Shiny application and provenance-aware output manifests.

Scientific interpretation boundaries are documented in the function help
pages and module-specific output metadata.

## Installation

Install the package from a local clone:

```r
install.packages(c("remotes", "BiocManager"))
remotes::install_local(".", dependencies = TRUE)
```

Core imports are available from CRAN. Optional modules use packages listed in
`Suggests`, including `Cardinal`, `EBImage`, `mgcv`, `nlme` and `nnls`.
Bioconductor dependencies can be installed explicitly when needed:

```r
BiocManager::install(c("Cardinal", "EBImage"))
```

For a future GitHub repository, installation will be:

```r
remotes::install_github("xia-lab/SpatialOmicsR", dependencies = TRUE)
```

The GitHub repository is named `SpatialOmicsR`; the installed R package and
library name remain `SpatialOmicsMSI` for compatibility.

## Minimal CSV workflow

```r
library(SpatialOmicsMSI)

pipeline <- run_spatial_metabolomics_pipeline(
  "path/to/pixel_feature_matrix.csv",
  tic_normalize = TRUE,
  do_log = TRUE
)

neighbors <- build_spatial_neighbors(
  pipeline$coordinates$x,
  pipeline$coordinates$y,
  method = "queen"
)

svg <- compute_spatially_variable_metabolites(
  pipeline$pixel_feature_matrix,
  coordinates = pipeline$coordinates,
  neighbors = neighbors,
  n_perm = 499,
  alternative = "two.sided",
  seed = 1
)
```

Input feature columns should be named `mz_*` or be numeric m/z column names.
Coordinates must be finite and `pixel_id` values must be unique.

## Executable vignette

The package contains one formal, executable vignette for the complete MSI
workflow. It runs real bundled pixels from import through Moran screening,
domains, niches, ROI summaries, exploratory statistics and study-matched
MALDI--LMD-LC-MS/MS validation:

```r
vignette("spatial_metabolomics_workflow", package = "SpatialOmicsMSI")
```

The vignette is the code-first reproducible analysis intended for review. The
Shiny application has a separate [GUI user guide](inst/shiny/spatial_pipeline/USER_GUIDE.md)
and is not treated as a second package vignette.

## Shiny application

After installation:

```r
SpatialOmicsMSI::run_spatial_app()
```

On a remote server, either listen locally and use an SSH tunnel, or explicitly
listen on all interfaces:

```r
SpatialOmicsMSI::run_spatial_app(
  host = "0.0.0.0", port = 3838, launch.browser = FALSE
)
```

From a source checkout:

```r
shiny::runApp("inst/shiny/spatial_pipeline")
```

Follow the [Shiny GUI user guide](inst/shiny/spatial_pipeline/USER_GUIDE.md) for
the ordered workflow, the two ROI paths, example-dataset coverage, remote-server
access, progress indicators and exports.

The root `app.R` also accepts `SPATIALOMICS_SHINY_HOST` and
`SPATIALOMICS_SHINY_PORT` environment variables.

Raw example datasets are intentionally excluded from Git. Their expected local
layout and provenance are documented in
[`inst/INSTALL_REAL_DATA.md`](inst/INSTALL_REAL_DATA.md).

## Reproducible scripts

Run scripts from the repository root. Local data paths are supplied through
environment variables rather than source-code edits. Copy `.env.example` as a
reference and set variables in the shell or a private `.Renviron` file.

```bash
export SPATIALOMICS_PROJECT_ROOT="$PWD"
export SPATIALOMICS_RAT_BRAIN_DIR="/path/to/rat_brain_data"
Rscript scripts/test_rat_brain_workflow.R
```

Never commit `.Renviron`, raw imzML/ibd files, downloaded databases, or generated
analysis outputs.

## Testing

Run the dependency-light regression suite:

```bash
Rscript scripts/run_tests.R
```

Tests follow the standard `testthat` layout under `tests/testthat/`. Small
synthetic fixtures may be committed there; tests using local raw datasets remain
explicitly opt-in.

Run a package check when all declared dependencies are installed:

```bash
R CMD build .
R CMD check --as-cran SpatialOmicsMSI_*.tar.gz
```

Tests that require optional packages or real datasets skip those paths when the
dependency/data is unavailable. GitHub Actions runs the same package check on
pushes and pull requests.

Slow tests against locally installed raw datasets are opt-in:

```bash
SPATIALOMICS_RUN_REAL_DATA_TESTS=true Rscript scripts/run_tests.R
```

## Project layout

```text
R/          Package functions
man/        User-facing function documentation
tests/      Regression and scientific-boundary tests
scripts/    Reproducible analysis entry points
inst/       Shiny application and data-installation instructions
vignettes/  One executable, code-first package workflow
```

## Citation and license

The software is released under the MIT license. Method-specific references are
provided in the corresponding help pages. If this package is used
in a publication, cite both the package version and the original methods that
the selected modules implement or approximate.

## Scope limitations

- BANKSY and HMRF modules are transparent approximations, not wrappers around
  the complete upstream implementations.
- Ligand-receptor results indicate expression/proximity-supported communication
  potential, not receptor activation or causal signaling.
- Metabolite identities require evidence beyond accurate mass and pathway hits.
- LC-MS-only CCS cannot confirm an MSI identity when MSI lacks an ion-mobility
  measurement. Vendor-neutral LCM coordinates require instrument-side
  calibration and validation before physical cutting.
- LMD geometry QC is read-only; the package does not silently split, buffer or
  simplify biological ROIs into manufacturer-specific cutting commands.
- Population claims require independent biological subjects.
