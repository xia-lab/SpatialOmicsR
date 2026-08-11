# Kasarla et al. 2025 kidney spatial-validation example

These compact CSV files are derived from Table S12 of Kasarla et al.,
*Analytical Chemistry* (2025), DOI: 10.1021/acs.analchem.5c00620.

The ROI tables compare MALDI-MSI and LMD-LC-MS/MS measurements across cortex,
medulla, and renal pelvis in serial mouse-kidney sections. They support the
ROI-level spatial-validation example in Shiny and the executable package
vignette. The compact pixel fixture is a fixed measured field reconstructed
from the matched METASPACE kidney acquisition for executable MSI examples. None
of these files represents the full five-organ raw-data deposit.

`kasarla_kidney_vignette_pixels.csv.gz` contains 2,500 measured pixels and 24
fixed real m/z features. `kasarla_kidney_vignette_features.csv` records their
source annotations, detection fractions, and post-TIC/log10 variances. Rebuild
both from the local public-data reconstruction with
`scripts/prepare_vignette_fixture.R`.

The reconstruction commands require the public `metaspace` Python client and
NumPy, followed by base R:

```text
python3 scripts/prepare_kasarla_metaspace_pixel_matrix.py \
  data_raw/kasarla2025/processed/kasarla_kidney_pixel_feature_matrix.csv.gz \
  --features data_raw/kasarla2025/processed/kasarla_kidney_metaspace_feature_metadata.csv
Rscript scripts/prepare_vignette_fixture.R
```

Source records:

- Article: https://doi.org/10.1021/acs.analchem.5c00620
- Raw-data accession: MSV000096852
- METASPACE project: https://metaspace2020.eu/project/swapna-2024
- Supporting spreadsheet: ACS Figshare file 54503595

The ACS Figshare supporting material is distributed under CC BY-NC 4.0. Keep
this attribution and license notice with redistributed derived tables.

## Selected matched kidney bundle

`kasarla_kidney_bundle_manifest.csv` records the minimal study-matched bundle.
The MALDI-MSI input is METASPACE dataset `2023-12-19_13h02m39s`, named
`231218_SK_labeledkidney_HEXANE_NH4OHwash_NEDC_MALDI_40um`. The orthogonal
measurements are the 12 open-format LMD-LC-MS/MS files from MSV000096852:
cortex, medulla, and renal pelvis, each with four replicates.

The selection is based on four convergent identifiers: Table S12 labels the
experiment `DgluKidney_Washed`; the article specifies negative-ion NEDC after
basic hexane/NH4OH washing; the METASPACE acquisition is the implicitly named
`Rep1` member of the dated Rep1/Rep2/Rep3 reproducibility series; and its kidney
outline agrees with the Rep1/Figure 4 ion-image geometry. Files from other wash
conditions, the positive-ion DHB acquisition, and the explicit Rep2/Rep3
reproducibility sections are not substituted for this spatial-validation input.

The repository does not supply a raw registered H&E image for this METASPACE
dataset through its public optical-image API. The H&E panel in the article is
therefore evidence for ROI interpretation, not a registration-ready image file.
