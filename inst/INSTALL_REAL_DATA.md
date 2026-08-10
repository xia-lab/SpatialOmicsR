# Installing and supplying the independent real datasets

Install with `R CMD INSTALL SpatialOmicsMSI_0.3.1.tar.gz`, then run
`SpatialOmicsMSI::run_spatial_app()`. Raw datasets are excluded from the tarball.

Expected server layout:

```
data_raw/mouse_brain_he_msi/metaspace_brain01/
  Brain01_Bregma-1-46_centroid.imzML
  Brain01_Bregma-1-46_centroid.ibd
  optical_brightfield.jpg
  optical_transform_api.json
  attribution_license.json
data_raw/full_tissue_mouse_brain/OMIX016317/
  OMIX016317-02.imzML
  OMIX016317-01.ibd
data_raw/mouse_brain_lcms/msv000090179/
  pos_mouse_female_brain_12w_1.mzML
  neg_mouse_female_brain_12w_1.mzML
  metadata_brain_pos.txt
  metadata_brain_neg.txt
```

- Dataset 1: METASPACE `2016-09-22_11h16m17s`, CC BY 4.0.
- Dataset 2: NGDC OMIX `OMIX016317`, sample `mbrain1_neg100`.
- Dataset 3: MassIVE `MSV000090179`, CC0, 12-week replicate-1 positive/negative mzML.

Download from each official accession, retain original filenames, verify provider
checksums when supplied, and record local MD5/SHA256. These are three independent
mouse-brain datasets for complementary technical validation, not a matched cohort
or same-sample cross-platform experiment.
