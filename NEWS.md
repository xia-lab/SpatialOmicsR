# SpatialOmicsMSI 0.3.1

- Reorganized Shiny into an ordered MSI workflow from data validation through
  ROI selection, domains, niches, statistics, spatial checks, matched omics,
  optional LC-MS/MS evidence, and auditable export.
- Added reusable domain/niche comparisons, responsive result panels, bounded
  progress feedback, workload warnings, and remote host/port configuration.
- Added Kasarla et al. kidney MSI--LMD-LC-MS/MS example data with source and
  licensing records, reconstruction scripts, and ROI-concordance validation.
- Added one executable package vignette using a compact real MSI fixture and a
  separate Shiny GUI user guide.
- Added regression coverage for Shiny presets, niche workflows, example-data
  routing, MetaboAnalyst exports, and the package Shiny launcher.
- Standardized CSV, shared-axis imzML, variable-axis imzML and mzML MS2 input.
- Added auditable background filtering, tissue masks and preprocessing QC.
- Added reusable spatial graphs and Moran, Geary and binSpect-inspired analyses.
- Added BANKSY-inspired H0 and ICM/Potts HMRF spatial clustering.
- Added multimodal cell-label preparation, niche and proximity analysis.
- Added ligand-receptor communication-potential sensitivity analyses.
- Added signed distance-to-domain rings, GAMs and Gaussian spatial LMMs.
- Added subject-level Wilcoxon analysis with fixed contrast-level pairing.
- Added real-data Shiny presets, provenance manifests and workflow documentation.
- Parameterized repository scripts through `SPATIALOMICS_*` environment variables.
- Added a dependency-light test runner and GitHub Actions R package checks.
