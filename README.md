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
- Automatic, coordinate, polygon and histology-guided ROI workflows.
- Histology-to-MSI affine registration and METASPACE optical transforms.
- Cell-label preparation from external labels, NNLS/cosine transfer or
  exploratory watershed morphology classes.
- Schürch-style neighborhood composition niches and permutation-based
  cell-type proximity enrichment.
- Expression-only and spatial ligand-receptor communication-potential analyses.
- Signed distance-to-domain rings and descriptive or subject-level GAMs.
- Subject-level t-test/Wilcoxon analysis and Gaussian spatial mixed models.
- Annotation, LC-MS and MS/MS evidence integration.
- MetaboAnalyst export/import, ORA/mummichog scripts and spatial back-mapping.
- A bundled Shiny application and provenance-aware output manifests.

See the [workflow vignette](vignettes/spatial_metabolomics_workflow.Rmd) for the
scientific workflow and interpretation boundaries.

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
remotes::install_github("OWNER/SpatialOmicsMSI", dependencies = TRUE)
```

Replace `OWNER` after the repository has been created; no Git remote is
currently configured in this working copy.

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

## Shiny application

After installation:

```r
SpatialOmicsMSI::run_spatial_app()
```

From a source checkout:

```r
shiny::runApp("inst/shiny/spatial_pipeline")
```

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
vignettes/  End-to-end workflow and interpretation guidance
inst/       Shiny application and data-installation instructions
```

## Citation and license

The software is released under the MIT license. Method-specific references are
provided in the corresponding help pages and vignette. If this package is used
in a publication, cite both the package version and the original methods that
the selected modules implement or approximate.

## Scope limitations

- BANKSY and HMRF modules are transparent approximations, not wrappers around
  the complete upstream implementations.
- Ligand-receptor results indicate expression/proximity-supported communication
  potential, not receptor activation or causal signaling.
- Metabolite identities require evidence beyond accurate mass and pathway hits.
- Population claims require independent biological subjects.
