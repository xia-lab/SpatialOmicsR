#!/usr/bin/env python3
"""Export a public METASPACE kidney dataset as a pixel-by-feature CSV.

This uses 10% FDR HMDB annotations and the first isotope image for each
annotation. Intensities are restored to the METASPACE image intensity scale;
TIC normalization remains an explicit downstream Shiny processing choice.
"""

import argparse
import csv
import gzip
from pathlib import Path

import numpy as np
from metaspace import SMInstance


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("output", type=Path)
    parser.add_argument("--features", type=Path)
    parser.add_argument("--dataset-id", default="2023-12-19_13h02m39s")
    args = parser.parse_args()

    dataset = SMInstance().dataset(id=args.dataset_id)
    mask = dataset.diagnostic("IMZML_METADATA")["images"][0]["image"].astype(bool)
    metadata = dataset.annotations(
        fdr=0.1,
        return_vals=("sumFormula", "adduct", "neutralLoss", "chemMod"),
    )
    ion_images = dataset.all_annotation_images(
        fdr=0.1,
        only_first_isotope=True,
        scale_intensity=True,
    )

    rows, columns = np.nonzero(mask)
    feature_records = []
    feature_values = []
    seen = {}
    for annotation, images in zip(metadata, ion_images):
        if not images or images[0] is None or images[0].shape != mask.shape:
            continue
        mz = float(images.peak(0))
        mz_key = f"{mz:.6f}"
        seen[mz_key] = seen.get(mz_key, 0) + 1
        suffix = "" if seen[mz_key] == 1 else f"__{seen[mz_key]}"
        column_name = f"mz_{mz_key}{suffix}"
        values = np.asarray(images[0], dtype=np.float64)[mask]
        values[~np.isfinite(values)] = 0.0
        feature_values.append(values)
        feature_records.append(
            [column_name, mz, annotation[0], annotation[1], annotation[2], annotation[3]]
        )

    if not feature_values:
        raise RuntimeError("METASPACE returned no usable annotation images")
    matrix = np.column_stack(feature_values)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    opener = gzip.open if args.output.suffix == ".gz" else open
    with opener(args.output, "wt", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["pixel_id", "x", "y"] + [x[0] for x in feature_records])
        for index, (row, column, values) in enumerate(zip(rows, columns, matrix), start=1):
            # IMZML metadata reports x starting at 1 and y starting at 2.
            writer.writerow([index, column + 1, row + 2] + values.tolist())

    feature_path = args.features or args.output.with_name("kasarla_kidney_feature_metadata.csv")
    with open(feature_path, "w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["column_name", "mz", "sum_formula", "adduct", "neutral_loss", "chem_mod"])
        writer.writerows(feature_records)
    print(f"Wrote {matrix.shape[0]} pixels x {matrix.shape[1]} features to {args.output}")


if __name__ == "__main__":
    main()
