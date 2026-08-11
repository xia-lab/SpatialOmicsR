# SpatialOmicsMSI Shiny user guide

This guide is for users who prefer the graphical interface. The reproducible,
code-first scientific analysis is the package vignette
`spatial_metabolomics_workflow`; this document is deliberately not a second
formal vignette.

## Start the app

After installing the package:

```r
SpatialOmicsMSI::run_spatial_app()
```

From a source checkout:

```r
shiny::runApp("inst/shiny/spatial_pipeline")
```

On a remote computer, start the app on the remote host and use an SSH tunnel.
For example, run the app remotely on loopback:

```r
SpatialOmicsMSI::run_spatial_app(
  host = "127.0.0.1", port = 3838, launch.browser = FALSE
)
```

Then run this on the local computer, replacing `user@server`:

```bash
ssh -L 3838:127.0.0.1:3838 user@server
```

Open `http://127.0.0.1:3838` in the local browser. Keep both the R process and
SSH connection open.

## Ordered workflow

The navigation is ordered to match the analysis rather than the package source
files.

1. **Data setup** — select a complete example or supply MSI, declare sample and
   ion-mode provenance, and validate the input.
2. **Input & provenance** — confirm the files and metadata that define the
   current session.
3. **Processing & QC** — create the pixel-feature matrix, inspect QC, and
   optionally compute Moran's I once for later screening and validation.
4. **Registration (optional)** — use a supplied METASPACE transform or a
   separately fitted histology-to-MSI transformation. Do not visually compare
   unmatched images as if they were registered.
5. **Optional pre-analysis ROI** — normally retain the complete tissue. Choose a
   prior ROI only if anatomy or acquisition coordinates defined it independently
   before MSI clustering.
6. **Domains & niches** — detect molecular-profile domains, review them, and
   optionally calculate local domain-composition niches.
7. **ROI summaries & statistics** — retain one or more reviewed domains as ROI,
   aggregate subregions, compare selected groups, and export statistical tables.
8. **Spatial validation** — reuse the cached Moran result to check spatial
   coherence and, when available, compare domains with an independently created
   label map.
9. **Matched transcriptomics & H&E** — inspect same-section companion files only
   after their coordinate relationship is declared and registration is valid.
10. **LC-MS/MS evidence (optional)** — evaluate chemical or ROI-concordance
    evidence for MSI findings; this is not a standalone LC-MS workflow.
11. **Results & export** — download auditable session outputs.

The blue **Working…** card appears only while the R server is executing an
action. The action button is disabled during that operation and restored after
success or error. A browser tab spinner without a Working card usually means the
browser is waiting for a blocked or disconnected R session; inspect the remote
terminal for the error rather than repeatedly clicking the action.

## ROI paths and niches

The default path is:

```text
complete tissue -> domain detection -> review map -> retain one or more domains as ROI
```

Two prior-ROI alternatives are available only when justified:

- **Registered histology ROI** transfers an independently drawn H&E polygon into
  MSI coordinates. It requires a validated transformation.
- **Direct MSI-coordinate ROI** selects a known rectangle in MSI x/y space. It
  does not register or interpret histology.

A niche is not another ROI algorithm. A domain describes the molecular profile
of one pixel; a niche describes the proportions of domains around that pixel.
Use the comparison selectors to compare two domains or two niches repeatedly.
Selections remain editable; running a comparison does not consume or replace
the labels.

## Example-dataset coverage

An example appears beside a function only when it can verify that function. A
complete-study preset can start at Data setup and keeps its companion files for
later steps; a module-only reference does not pretend to be pixel-level MSI.

| Example | Core MSI/QC | Moran/domains/niches/ROI | Registration or matched H&E/ST | LC-MS/MS role |
|---|---:|---:|---|---|
| Brain01 | Yes | Core MSI possible | Supplied METASPACE optical transform; the image is brightfield, not H&E | None |
| OMIX016317 | Yes | Yes; full-tissue MSI example | None | None |
| MSIflow UPEC_12 | Yes | Yes | None in the preset | Same-study chemical evidence |
| Vicari 2024 mPD1/mPD3/mPD4 | Yes | Yes after MSI processing | Matched Visium/H&E files are retained, but spatial integration requires validated registration | None in this preset |
| Kasarla 2025 kidney | Yes | Yes on the bundled pixel-level MSI | No registration-ready H&E is bundled | Study-matched serial-section LMD-LC-MS/MS ROI concordance |
| MSV000090179 | No | No | No | External-reference LC-MS/MS evidence only |

The Vicari preset retaining Visium files does not by itself prove successful
deconvolution or ligand--receptor analysis. Those functions require compatible
expression matrices, reference profiles and validated coordinate handling. The
UI must report missing prerequisites rather than marking the modules as tested.

## Moran's I in two places

Moran's I is calculated per m/z feature, not on domain labels.

- Before domains, it can screen spatially structured features for clustering.
- After ROI statistics, the same cached table can flag whether differential
  features have coherent spatial patterns.

The second use is supporting evidence, not independent validation when Moran
features were already used to construct the domains. Runtime grows roughly with
pixels × graph edges × features × permutations. The workload warning is
dataset-independent: it does not silently change requested permutations or
special-case one example.

## MetaboAnalyst and lipid pathway files

The ROI summary export is a sample-by-feature table for MetaboAnalyst statistical
analysis. Do not upload coordinate columns or internal grid fields as features.

The one-column lipid-name file is for lipid enrichment/pathway modules. It can
only be generated after reviewed annotations have mapped MSI m/z features to
lipid names or supported identifiers. MetaboAnalyst can map accepted lipid names,
but it cannot turn arbitrary m/z values into confirmed lipid identities without
annotation evidence. Download the mapping audit beside the upload file.

## Interpreting the displayed images

These bundled images illustrate correctly sized result panels; the app generates
or loads the corresponding outputs only for the active compatible dataset.

OMIX full-field tissue gate:

![OMIX tissue gate](www/omix016317_full_field_tissue_gate.png)

OMIX data-driven domains:

![OMIX domains](www/omix016317_tissue_only_domains.png)

Brain01 supplied-transform registration diagnostic:

![Brain01 registration](www/brain01_live_measurement_overlay.png)

Images and tables use bounded responsive containers so that a large source image
does not push controls below the viewport. Technical inventories are collapsed
or scrollable rather than overlaid on the image.

## Scientific boundaries

- Pixels and subregion tiles are not biological replicates.
- Data-driven domains and composition niches are not anatomical ground truth.
- LC-MS/MS evidence strength depends on whether the data are specimen matched,
  study matched, or an external reference.
- Accurate mass and pathway membership do not confirm molecular identity.
- Cell-type, cellular-niche and expression ligand--receptor analyses require
  suitable spatial resolution and matched biological inputs; function
  availability alone is not validation on MSI data.
