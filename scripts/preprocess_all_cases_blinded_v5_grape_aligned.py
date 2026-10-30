#!/usr/bin/env python3
"""Outcome-blinded GRAPE-aligned CT preprocessing for the CT2 cohort.

This is a radiomics sensitivity pipeline inspired by Hu et al., Nature Medicine
2025: resample to 0.717 x 0.717 x 5.0 mm and apply a training-set global
intensity standardisation. It is not claimed to reproduce GRAPE exactly: the
paper's original mean/SD and implementation are not released. Fixed source
masks are never edited; NN-resampled masks here are disposable working copies.
"""
from plain_paths import plain_path

import argparse
import csv
import json
import math
import os
import re
import shutil
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime
from pathlib import Path

import numpy as np
import openpyxl
import SimpleITK as sitk


PROJECT = Path(plain_path("PROJECT", ""))
SOURCE_MAPPING = PROJECT / "汇总/影像组学_训练验证集_当前提取版本_待一致化/影像组学训练验证集_病例映射.csv"
PLAIN_DICOM_ROOT = PROJECT / "CT2/平扫期"
DEFAULT_OUTPUT = PROJECT / "outputs/CT2_ALL_CASES_blinded_v5_grape_aligned"
TARGET_SPACING = (0.717, 0.717, 5.0)
MARGIN_XY, MARGIN_Z = 25, 3
EXCEPTION_CASE = "PRIVATE_CASE_001"


def same_geometry(first, second):
    return (
        first.GetSize() == second.GetSize()
        and np.allclose(first.GetSpacing(), second.GetSpacing())
        and np.allclose(first.GetOrigin(), second.GetOrigin())
        and np.allclose(first.GetDirection(), second.GetDirection())
    )


def label_voxels(mask):
    return int(np.count_nonzero(sitk.GetArrayViewFromImage(mask) == 1))


def physical_volume_ml(mask):
    return label_voxels(mask) * float(np.prod(mask.GetSpacing())) / 1000.0


def read_wide_paths(path, case_id):
    workbook = openpyxl.load_workbook(path, read_only=True, data_only=True)
    try:
        rows = {str(row[1]): row for row in workbook.active.iter_rows(min_row=2, max_row=3, values_only=True)}
    finally:
        workbook.close()
    if set(rows) != {"stomach", "tumor"} or any(str(row[0]) != case_id for row in rows.values()):
        raise ValueError(f"Unexpected wide metadata: {case_id}")
    return {
        "image": Path(str(rows["stomach"][2])),
        "stomach_source": Path(str(rows["stomach"][3])),
        "tumor_source": Path(str(rows["tumor"][3])),
        "stomach_aligned": Path(str(rows["stomach"][4])),
        "tumor_aligned": Path(str(rows["tumor"][4])),
    }


def read_dicom(folder):
    series = sitk.ImageSeriesReader.GetGDCMSeriesIDs(str(folder)) or []
    if len(series) != 1:
        raise ValueError(f"Expected one DICOM series, got {len(series)}")
    names = sitk.ImageSeriesReader.GetGDCMSeriesFileNames(str(folder), series[0])
    reader = sitk.ImageSeriesReader()
    reader.SetFileNames(names)
    return reader.Execute(), len(names)


def crop_region(stomach):
    z, y, x = np.where(sitk.GetArrayViewFromImage(stomach) == 1)
    if not len(x):
        raise ValueError("Empty fixed stomach mask")
    sx, sy, sz = stomach.GetSize()
    start = (max(0, int(x.min()) - MARGIN_XY), max(0, int(y.min()) - MARGIN_XY), max(0, int(z.min()) - MARGIN_Z))
    stop = (min(sx, int(x.max()) + 1 + MARGIN_XY), min(sy, int(y.max()) + 1 + MARGIN_XY), min(sz, int(z.max()) + 1 + MARGIN_Z))
    return start, tuple(stop[index] - start[index] for index in range(3))


def resample_to_spacing(image, interpolator, default_value):
    old_spacing = image.GetSpacing()
    old_size = image.GetSize()
    new_size = tuple(
        max(1, int(round(old_size[index] * old_spacing[index] / TARGET_SPACING[index])))
        for index in range(3)
    )
    return sitk.Resample(
        image,
        new_size,
        sitk.Transform(),
        interpolator,
        image.GetOrigin(),
        TARGET_SPACING,
        image.GetDirection(),
        default_value,
        sitk.sitkFloat32 if interpolator != sitk.sitkNearestNeighbor else image.GetPixelID(),
    )


def prepare_case(record, stage_root):
    case_id = record["case_id"]
    case_dir = stage_root / case_id
    targets = {"ct": case_dir / "CT_unnormalized.nii.gz", "stomach": case_dir / "stomach_mask.nii.gz", "tumor": case_dir / "tumor_mask.nii.gz"}
    if all(path.is_file() for path in targets.values()):
        return {"case_id": case_id, "status": "existing_valid", "split": record["split"]}

    paths = read_wide_paths(Path(record["wide_path"]), case_id)
    if case_id == EXCEPTION_CASE:
        image = sitk.ReadImage(str(paths["image"]))
        stomach = sitk.ReadImage(str(paths["stomach_aligned"]))
        tumor = sitk.ReadImage(str(paths["tumor_aligned"]))
        source = "exact_presegmentation_model_input_due_raw_dicom_coordinate_reset"
        slice_count = None
    else:
        stomach = sitk.ReadImage(str(paths["stomach_source"]))
        tumor = sitk.ReadImage(str(paths["tumor_source"]))
        if record["plain_dicom"]:
            image, slice_count = read_dicom(Path(record["plain_dicom"]))
            source = "raw_plain_dicom"
        else:
            image, slice_count = sitk.ReadImage(str(paths["image"])), None
            source = "raw_control_nifti"
        if not same_geometry(stomach, tumor):
            raise ValueError("Fixed stomach/tumor masks do not share geometry")
        if not same_geometry(image, stomach):
            image = sitk.Resample(image, stomach, sitk.Transform(), sitk.sitkBSpline, -1024.0, sitk.sitkFloat32)

    before_stomach_ml = physical_volume_ml(stomach)
    before_tumor_ml = physical_volume_ml(tumor)
    start, size = crop_region(stomach)
    image = sitk.RegionOfInterest(sitk.Cast(image, sitk.sitkFloat32), size=size, index=start)
    stomach = sitk.RegionOfInterest(stomach, size=size, index=start)
    tumor = sitk.RegionOfInterest(tumor, size=size, index=start)
    if label_voxels(stomach) == 0 or label_voxels(tumor) == 0:
        raise ValueError("ROI was truncated during crop")

    image = resample_to_spacing(image, sitk.sitkBSpline, -1024.0)
    stomach = resample_to_spacing(stomach, sitk.sitkNearestNeighbor, 0)
    tumor = resample_to_spacing(tumor, sitk.sitkNearestNeighbor, 0)
    if label_voxels(stomach) == 0 or label_voxels(tumor) == 0:
        raise ValueError("ROI became empty after working-grid resampling")
    case_dir.mkdir(parents=True, exist_ok=True)
    sitk.WriteImage(image, str(targets["ct"]), useCompression=True)
    sitk.WriteImage(stomach, str(targets["stomach"]), useCompression=True)
    sitk.WriteImage(tumor, str(targets["tumor"]), useCompression=True)

    values = np.asarray(sitk.GetArrayViewFromImage(image), dtype=np.float64)
    stomach_values = values[np.asarray(sitk.GetArrayViewFromImage(stomach)) == 1]
    metadata = {
        "case_id": case_id,
        "source_kind": source,
        "raw_dicom_slices": slice_count,
        "working_spacing_mm": TARGET_SPACING,
        "crop_source": "fixed stomach mask bounding box plus fixed margin",
        "mask_operation": "source masks unchanged; temporary working masks NN-resampled only",
        "stomach_source_ml": before_stomach_ml,
        "tumor_source_ml": before_tumor_ml,
        "stomach_working_ml": physical_volume_ml(stomach),
        "tumor_working_ml": physical_volume_ml(tumor),
    }
    (case_dir / "preprocessing_metadata.json").write_text(json.dumps(metadata, ensure_ascii=False, indent=2), encoding="utf-8")
    return {
        "case_id": case_id,
        "status": "completed",
        "split": record["split"],
        "stomach_n": int(stomach_values.size),
        "stomach_sum": float(stomach_values.sum()),
        "stomach_sumsq": float(np.square(stomach_values).sum()),
    }


def hard_link(source, destination):
    os.link(source, destination)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-root", default=str(DEFAULT_OUTPUT))
    parser.add_argument("--workers", type=int, default=3)
    parser.add_argument("--copy-script-to")
    args = parser.parse_args()
    output_root = Path(args.output_root)
    if output_root.exists():
        raise FileExistsError(f"Refusing to overwrite {output_root}")
    with SOURCE_MAPPING.open(encoding="utf-8-sig") as handle:
        source_rows = list(csv.DictReader(handle))
    dicom_dirs = {match.group(0): str(path) for path in PLAIN_DICOM_ROOT.iterdir() if path.is_dir() and (match := re.search(r"P\d+", path.name))}
    records = sorted(
        ({"case_id": row["case_id"], "wide_path": row["wide_path"], "split": row["split"], "plain_dicom": dicom_dirs.get(row["case_id"], "")} for row in source_rows),
        key=lambda row: row["case_id"],
    )
    if len(records) != 820 or len({row["case_id"] for row in records}) != 820:
        raise ValueError("Expected 820 unique cases")
    if {row["split"] for row in records} != {"训练集", "验证集"}:
        raise ValueError("Only the fixed train/validation split may be read; outcome is never read")

    stage_root = output_root / "stage_unnormalized"
    stage_root.mkdir(parents=True)
    batch_dir = output_root / f"batch_{datetime.now().strftime('%Y%m%d-%H%M%S')}"
    batch_dir.mkdir()
    results = []
    with ThreadPoolExecutor(max_workers=args.workers) as executor:
        futures = {executor.submit(prepare_case, record, stage_root): record["case_id"] for record in records}
        for index, future in enumerate(as_completed(futures), start=1):
            case_id = futures[future]
            try:
                result = future.result()
            except Exception as exc:
                result = {"case_id": case_id, "status": "failed", "error": str(exc)}
            results.append(result)
            print(f"[{index}/820] {case_id} {result['status']}", flush=True)
    failed = [row for row in results if row["status"] == "failed"]
    if failed:
        raise RuntimeError(f"Stage failed for {len(failed)} cases; v5 final data were not created")

    train_stats = [row for row in results if row["split"] == "训练集"]
    total_n = sum(row["stomach_n"] for row in train_stats)
    total_sum = sum(row["stomach_sum"] for row in train_stats)
    total_sumsq = sum(row["stomach_sumsq"] for row in train_stats)
    mean = total_sum / total_n
    variance = max(0.0, total_sumsq / total_n - mean * mean)
    standard_deviation = math.sqrt(variance)
    if not np.isfinite(mean) or not np.isfinite(standard_deviation) or standard_deviation <= 0:
        raise ValueError("Invalid training-only intensity statistics")

    final_rows = []
    for index, record in enumerate(records, start=1):
        case_id = record["case_id"]
        source_dir = stage_root / case_id
        target_dir = output_root / case_id
        target_dir.mkdir()
        image = sitk.ReadImage(str(source_dir / "CT_unnormalized.nii.gz"))
        normalized = sitk.ShiftScale(sitk.Cast(image, sitk.sitkFloat32), shift=-mean, scale=1.0 / standard_deviation)
        ct = target_dir / "CT.nii.gz"
        stomach = target_dir / "stomach_mask.nii.gz"
        tumor = target_dir / "tumor_mask.nii.gz"
        sitk.WriteImage(normalized, str(ct), useCompression=True)
        hard_link(source_dir / "stomach_mask.nii.gz", stomach)
        hard_link(source_dir / "tumor_mask.nii.gz", tumor)
        final_rows.append({"case_id": case_id, "ct": str(ct), "stomach_mask": str(stomach), "tumor_mask": str(tumor)})
        if index % 100 == 0:
            print(f"Normalized {index}/820", flush=True)
    with (output_root / "all_cases_manifest_blinded.csv").open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.DictWriter(handle, fieldnames=["case_id", "ct", "stomach_mask", "tumor_mask"])
        writer.writeheader()
        writer.writerows(final_rows)
    results.sort(key=lambda row: row["case_id"])
    fields = sorted({key for row in results for key in row})
    with (batch_dir / "stage_results.csv").open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(results)
    summary = {
        "protocol": "CT2-ALL-GRAPE-aligned-v5-sensitivity",
        "cases": 820,
        "outcome_not_read": True,
        "split_used_only_for_training_global_intensity_statistics": True,
        "target_spacing_mm": TARGET_SPACING,
        "training_stomach_roi_global_mean_HU": mean,
        "training_stomach_roi_global_sd_HU": standard_deviation,
        "source_masks_modified": False,
        "exception_case": EXCEPTION_CASE,
        "manifest": str(output_root / "all_cases_manifest_blinded.csv"),
    }
    (batch_dir / "summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")
    (output_root / "summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")
    current = Path(__file__).resolve()
    shutil.copy2(current, output_root / current.name)
    if args.copy_script_to:
        target = Path(args.copy_script_to)
        target.mkdir(parents=True, exist_ok=True)
        shutil.copy2(current, target / current.name)
    print(summary)


if __name__ == "__main__":
    main()
