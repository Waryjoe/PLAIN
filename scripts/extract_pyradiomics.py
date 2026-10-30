#!/usr/bin/env python3
import argparse
import csv
import json
from pathlib import Path

import numpy as np
import pandas as pd
import SimpleITK as sitk
from radiomics import featureextractor


def primitive(value):
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, (tuple, list)):
        return json.dumps([primitive(v) for v in value], ensure_ascii=False)
    if hasattr(value, "tolist"):
        return value.tolist()
    return value


def image_info(path):
    img = sitk.ReadImage(str(path))
    arr = sitk.GetArrayFromImage(img)
    return {
        "path": str(path),
        "size_xyz": list(img.GetSize()),
        "spacing_xyz": list(img.GetSpacing()),
        "origin_xyz": list(img.GetOrigin()),
        "direction": list(img.GetDirection()),
        "pixel_type": img.GetPixelIDTypeAsString(),
        "min": float(np.nanmin(arr)),
        "max": float(np.nanmax(arr)),
        "nonzero_voxels": int(np.count_nonzero(arr)),
        "unique_values_first_20": [primitive(v) for v in np.unique(arr)[:20]],
    }


def resample_mask_to_image(mask, reference):
    return sitk.Resample(
        mask,
        reference,
        sitk.Transform(),
        sitk.sitkNearestNeighbor,
        0,
        sitk.sitkUInt8,
    )


def mask_voxel_summary(mask):
    arr = sitk.GetArrayFromImage(mask)
    voxels = int(np.count_nonzero(arr == 1))
    summary = {"label_1_voxels": voxels}
    nz = np.argwhere(arr == 1)
    if nz.size:
        summary["bbox_zyx_min"] = [int(v) for v in nz.min(axis=0)]
        summary["bbox_zyx_max"] = [int(v) for v in nz.max(axis=0)]
    return summary


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--case-id", required=True)
    parser.add_argument("--image", required=True)
    parser.add_argument("--stomach-mask", required=True)
    parser.add_argument("--tumor-mask", required=True)
    parser.add_argument("--params", required=True)
    parser.add_argument("--outdir", required=True)
    args = parser.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    aligned_dir = outdir / "aligned_masks"
    aligned_dir.mkdir(parents=True, exist_ok=True)

    image_path = Path(args.image)
    image = sitk.ReadImage(str(image_path))
    masks = {
        "stomach": Path(args.stomach_mask),
        "tumor": Path(args.tumor_mask),
    }

    extractor = featureextractor.RadiomicsFeatureExtractor(str(args.params))
    rows = []
    qc = {
        "case_id": args.case_id,
        "image": image_info(image_path),
        "params": str(Path(args.params)),
        "rois": {},
        "pyradiomics_version": getattr(__import__("radiomics"), "__version__", "unknown"),
    }

    for roi, mask_path in masks.items():
        source_mask = sitk.ReadImage(str(mask_path))
        aligned_mask = resample_mask_to_image(source_mask, image)
        aligned_path = aligned_dir / f"{args.case_id}_{roi}_mask_aligned_to_image.nii.gz"
        sitk.WriteImage(aligned_mask, str(aligned_path))

        qc["rois"][roi] = {
            "source_mask": image_info(mask_path),
            "source_label_1": mask_voxel_summary(source_mask),
            "aligned_mask": image_info(aligned_path),
            "aligned_label_1": mask_voxel_summary(aligned_mask),
        }

        result = extractor.execute(str(image_path), str(aligned_path), label=1)
        row = {
            "case_id": args.case_id,
            "roi": roi,
            "image_path": str(image_path),
            "source_mask_path": str(mask_path),
            "aligned_mask_path": str(aligned_path),
            "source_label_1_voxels": qc["rois"][roi]["source_label_1"]["label_1_voxels"],
            "aligned_label_1_voxels": qc["rois"][roi]["aligned_label_1"]["label_1_voxels"],
        }
        for key, value in result.items():
            row[key] = primitive(value)
        rows.append(row)

    wide = pd.DataFrame(rows)
    wide_path = outdir / f"{args.case_id}_radiomics_features_wide.csv"
    wide.to_csv(wide_path, index=False, quoting=csv.QUOTE_MINIMAL)

    feature_prefixes = (
        "original",
        "log-sigma",
        "wavelet",
        "square",
        "squareroot",
        "logarithm",
        "exponential",
        "gradient",
        "lbp",
    )
    feature_columns = [c for c in wide.columns if c.startswith(feature_prefixes)]
    features_only = wide[
        [
            "case_id",
            "roi",
            "image_path",
            "source_mask_path",
            "aligned_mask_path",
            "source_label_1_voxels",
            "aligned_label_1_voxels",
        ]
        + feature_columns
    ]
    features_only_path = outdir / f"{args.case_id}_radiomics_features_only_wide.csv"
    features_only.to_csv(features_only_path, index=False, quoting=csv.QUOTE_MINIMAL)

    long_rows = []
    metadata_columns = {
        "case_id",
        "roi",
        "image_path",
        "source_mask_path",
        "aligned_mask_path",
        "source_label_1_voxels",
        "aligned_label_1_voxels",
    }
    for row in rows:
        for key, value in row.items():
            if key in metadata_columns:
                continue
            long_rows.append(
                {
                    "case_id": row["case_id"],
                    "roi": row["roi"],
                    "feature_name": key,
                    "value": value,
                }
            )

    long = pd.DataFrame(long_rows)
    long_path = outdir / f"{args.case_id}_radiomics_features_long.csv"
    long.to_csv(long_path, index=False, quoting=csv.QUOTE_MINIMAL)

    features_only_long = long[long["feature_name"].isin(feature_columns)]
    features_only_long_path = outdir / f"{args.case_id}_radiomics_features_only_long.csv"
    features_only_long.to_csv(features_only_long_path, index=False, quoting=csv.QUOTE_MINIMAL)

    qc_path = outdir / f"{args.case_id}_radiomics_qc.json"
    qc["outputs"] = {
        "wide_csv": str(wide_path),
        "features_only_wide_csv": str(features_only_path),
        "long_csv": str(long_path),
        "features_only_long_csv": str(features_only_long_path),
        "aligned_mask_dir": str(aligned_dir),
    }
    qc_path.write_text(json.dumps(qc, ensure_ascii=False, indent=2), encoding="utf-8")

    print(json.dumps({
        "case_id": args.case_id,
        "n_roi_rows": len(wide),
        "n_feature_columns": len(feature_columns),
        "wide_csv": str(wide_path),
        "features_only_wide_csv": str(features_only_path),
        "long_csv": str(long_path),
        "features_only_long_csv": str(features_only_long_path),
        "qc_json": str(qc_path),
    }, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
