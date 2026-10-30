#!/usr/bin/env python3
"""Batch-extract QC-validated CT2 v5 GRAPE-aligned radiomics features."""
from plain_paths import plain_path

import argparse
import csv
import hashlib
import json
import shutil
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd


PROJECT = Path(plain_path("PROJECT", ""))
MANIFEST = PROJECT / "outputs/CT2_ALL_CASES_blinded_v5_grape_aligned/all_cases_manifest_blinded.csv"
DEFAULT_EXCEL_ROOT = PROJECT / "CT2-EXCEL/GRAPE-v5"
DEFAULT_WORK_ROOT = PROJECT / "outputs/ct2_pyradiomics_v5_grape_aligned"
EXPECTED_FEATURES = 1595
METADATA = ["case_id", "roi", "image_path", "source_mask_path", "aligned_mask_path", "source_label_1_voxels", "aligned_label_1_voxels"]


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def validate(wide, long):
    if wide.shape[0] != 2 or set(wide["roi"]) != {"stomach", "tumor"}:
        raise ValueError("Expected exactly one stomach and one tumor row")
    if list(wide.columns[:len(METADATA)]) != METADATA:
        raise ValueError("Unexpected metadata schema")
    features = list(wide.columns[len(METADATA):])
    if len(features) != EXPECTED_FEATURES or any(name.startswith("diagnostics_") for name in features):
        raise ValueError(f"Expected {EXPECTED_FEATURES} pure feature columns")
    if not np.isfinite(wide[features].apply(pd.to_numeric, errors="coerce").to_numpy(float)).all():
        raise ValueError("Non-finite wide feature")
    if list(long.columns) != ["case_id", "roi", "feature_name", "value"] or long.shape[0] != EXPECTED_FEATURES * 2:
        raise ValueError("Unexpected long schema")
    if not np.isfinite(pd.to_numeric(long["value"], errors="coerce").to_numpy(float)).all():
        raise ValueError("Non-finite long feature")


def process(record, extractor, config, work_root, wide_dir, long_dir, logs_dir):
    case_id = record["case_id"]
    wide_xlsx = wide_dir / f"{case_id}_radiomics_features_GRAPE_v5_wide.xlsx"
    long_xlsx = long_dir / f"{case_id}_radiomics_features_GRAPE_v5_long.xlsx"
    if wide_xlsx.exists() or long_xlsx.exists():
        raise FileExistsError(f"Refusing to overwrite existing v5 output for {case_id}")
    workdir = work_root / case_id
    log_path = logs_dir / f"{case_id}.log"
    command = [
        sys.executable, str(extractor), "--case-id", case_id,
        "--image", record["ct"], "--stomach-mask", record["stomach_mask"], "--tumor-mask", record["tumor_mask"],
        "--params", str(config), "--outdir", str(workdir),
    ]
    with log_path.open("w", encoding="utf-8") as log:
        subprocess.run(command, check=True, stdout=log, stderr=subprocess.STDOUT)
    wide = pd.read_csv(workdir / f"{case_id}_radiomics_features_only_wide.csv")
    long = pd.read_csv(workdir / f"{case_id}_radiomics_features_only_long.csv")
    validate(wide, long)
    wide.to_excel(wide_xlsx, index=False, engine="openpyxl")
    long.to_excel(long_xlsx, index=False, engine="openpyxl")
    validate(pd.read_excel(wide_xlsx, engine="openpyxl"), pd.read_excel(long_xlsx, engine="openpyxl"))
    return {"case_id": case_id, "status": "completed", "wide_xlsx": str(wide_xlsx), "long_xlsx": str(long_xlsx), "features_per_roi": EXPECTED_FEATURES}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--excel-root", default=str(DEFAULT_EXCEL_ROOT))
    parser.add_argument("--work-root", default=str(DEFAULT_WORK_ROOT))
    parser.add_argument("--copy-script-to")
    args = parser.parse_args()
    project_root = Path(__file__).resolve().parent.parent
    config = project_root / "configs/ct2_pyradiomics_v5_grape_aligned_1595.yaml"
    extractor = project_root / "scripts/extract_pyradiomics.py"
    if not config.is_file() or not extractor.is_file(): raise FileNotFoundError("Missing v5 config or core extractor")
    with MANIFEST.open(encoding="utf-8-sig") as handle: records = list(csv.DictReader(handle))
    if len(records) != 820 or set(records[0]) != {"case_id", "ct", "stomach_mask", "tumor_mask"}: raise ValueError("Expected 820-case blind v5 manifest")
    excel_root, work_root = Path(args.excel_root), Path(args.work_root)
    wide_dir, long_dir, logs_dir = excel_root / "wide", excel_root / "long", work_root / "logs"
    for directory in (wide_dir, long_dir, logs_dir): directory.mkdir(parents=True, exist_ok=True)
    batch_dir = work_root / f"batch_{datetime.now().strftime('%Y%m%d-%H%M%S')}"; batch_dir.mkdir(parents=True)
    results = []
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        futures = {pool.submit(process, record, extractor, config, work_root, wide_dir, long_dir, logs_dir): record["case_id"] for record in records}
        for number, future in enumerate(as_completed(futures), start=1):
            case_id = futures[future]
            try: result = future.result()
            except Exception as error: result = {"case_id": case_id, "status": "failed", "error": str(error)}
            results.append(result); print(f"[{number}/820] {case_id} {result['status']}", flush=True)
    results.sort(key=lambda row: row["case_id"])
    fields = sorted({key for row in results for key in row})
    with (batch_dir / "extraction_results.csv").open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields); writer.writeheader(); writer.writerows(results)
    failed = [row for row in results if row["status"] == "failed"]
    summary = {"protocol": "CT2-PyRadiomics-v5-GRAPE-aligned", "cases": 820, "completed": 820-len(failed), "failed": len(failed), "features_per_roi": EXPECTED_FEATURES, "config_sha256": sha256(config), "extractor_sha256": sha256(extractor), "input_manifest": str(MANIFEST), "outcome_not_read": True}
    (batch_dir / "summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")
    shutil.copy2(config, batch_dir / config.name); shutil.copy2(Path(__file__), batch_dir / Path(__file__).name)
    if args.copy_script_to:
        destination = Path(args.copy_script_to); destination.mkdir(parents=True, exist_ok=True)
        shutil.copy2(Path(__file__), destination / Path(__file__).name); shutil.copy2(config, destination / config.name)
    print(summary)
    if failed: raise RuntimeError(f"Extraction failed for {len(failed)} cases")


if __name__ == "__main__": main()
