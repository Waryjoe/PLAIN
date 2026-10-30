#!/usr/bin/env python3
"""Create fixed train/validation matrices from blinded CT2 v5 radiomics output.

The v5 feature extraction was performed without the anastomotic-leakage
outcome.  This script performs the first outcome-aware operation: it joins the
already-fixed clinical split to per-patient wide files, prefixes features by
ROI, and validates every workbook before writing the analysis matrices.
"""

from __future__ import annotations
from plain_paths import plain_path

import csv
import json
import math
import shutil
from collections import Counter
from pathlib import Path

import openpyxl
from openpyxl import Workbook
from openpyxl.styles import Alignment, Font, PatternFill
from openpyxl.utils import get_column_letter


PROJECT = Path(plain_path("PROJECT", ""))
SOURCE_MAPPING = PROJECT / "汇总/影像组学_训练验证集_统一输入v2_fullroi/影像组学训练验证集_病例映射.csv"
WIDE_DIR = PROJECT / "CT2-EXCEL/GRAPE-v5/wide"
LONG_DIR = PROJECT / "CT2-EXCEL/GRAPE-v5/long"
BATCH_SUMMARY = PROJECT / "outputs/ct2_pyradiomics_v5_grape_aligned/batch_20260808-160958/summary.json"
DEFAULT_OUTPUT = PROJECT / "汇总/影像组学_训练验证集_v5_GRAPE对齐"

ROI_FEATURE_COUNT = 1595
METADATA_COLUMNS = [
    "case_id",
    "roi",
    "image_path",
    "source_mask_path",
    "aligned_mask_path",
    "source_label_1_voxels",
    "aligned_label_1_voxels",
]


def style_header(ws):
    fill = PatternFill("solid", fgColor="1F4E78")
    font = Font(color="FFFFFF", bold=True)
    for cell in ws[1]:
        cell.fill = fill
        cell.font = font
        cell.alignment = Alignment(horizontal="center", vertical="center")


def read_mapping():
    rows = []
    with SOURCE_MAPPING.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        expected = ["hospital_id", "case_id", "outcome", "split", "wide_path", "long_path"]
        if reader.fieldnames != expected:
            raise ValueError(f"Unexpected mapping schema: {reader.fieldnames}")
        for row in reader:
            rows.append({
                "hospital_id": str(row["hospital_id"]).zfill(8),
                "case_id": str(row["case_id"]),
                "outcome": int(row["outcome"]),
                "split": str(row["split"]),
            })
    if len(rows) != 820 or len({row["case_id"] for row in rows}) != 820:
        raise ValueError("Fixed mapping must contain 820 unique cases")
    observed = Counter((row["split"], row["outcome"]) for row in rows)
    expected = Counter({("训练集", 0): 530, ("训练集", 1): 44, ("验证集", 0): 227, ("验证集", 1): 19})
    if observed != expected:
        raise ValueError(f"Unexpected fixed split: {observed}")
    return rows


def source_paths(case_id):
    wide_matches = sorted(WIDE_DIR.glob(f"{case_id}_radiomics_features_GRAPE_v5_wide.xlsx"))
    long_matches = sorted(LONG_DIR.glob(f"{case_id}_radiomics_features_GRAPE_v5_long.xlsx"))
    if len(wide_matches) != 1 or len(long_matches) != 1:
        raise FileNotFoundError(
            f"Expected one v5 wide/long workbook for {case_id}; "
            f"wide={wide_matches}, long={long_matches}"
        )
    return wide_matches[0], long_matches[0]


def read_patient_features(case_id, expected_features=None):
    wide_path, _ = source_paths(case_id)
    wb = openpyxl.load_workbook(wide_path, read_only=True, data_only=True)
    try:
        if len(wb.sheetnames) != 1:
            raise ValueError(f"Unexpected sheets in {wide_path}: {wb.sheetnames}")
        ws = wb[wb.sheetnames[0]]
        header = list(next(ws.values))
        if header[:len(METADATA_COLUMNS)] != METADATA_COLUMNS:
            raise ValueError(f"Metadata schema mismatch: {wide_path}")
        features = header[len(METADATA_COLUMNS):]
        if len(features) != ROI_FEATURE_COUNT or any(not isinstance(x, str) for x in features):
            raise ValueError(f"Expected {ROI_FEATURE_COUNT} feature names in {wide_path}")
        if expected_features is not None and features != expected_features:
            raise ValueError(f"Feature header/order drift: {wide_path}")
        roi_values = {}
        for row in ws.iter_rows(min_row=2, values_only=True):
            if row[0] != case_id:
                raise ValueError(f"case_id mismatch in {wide_path}: {row[0]}")
            roi = row[1]
            values = list(row[len(METADATA_COLUMNS):])
            if roi not in {"stomach", "tumor"} or roi in roi_values:
                raise ValueError(f"Unexpected ROI layout in {wide_path}")
            if len(values) != ROI_FEATURE_COUNT or any(
                value is None or not isinstance(value, (int, float)) or not math.isfinite(value)
                for value in values
            ):
                raise ValueError(f"Invalid feature values in {wide_path}, ROI={roi}")
            roi_values[roi] = values
        if set(roi_values) != {"stomach", "tumor"}:
            raise ValueError(f"Both stomach and tumor rows required: {wide_path}")
        return features, roi_values["stomach"], roi_values["tumor"]
    finally:
        wb.close()


def write_matrix(path, rows):
    if path.exists():
        raise FileExistsError(f"Refusing to overwrite: {path}")
    features, _, _ = read_patient_features(rows[0]["case_id"])
    wb = Workbook(write_only=True)
    ws = wb.create_sheet("影像组学特征")
    header = [
        "住院号", "case_id", "吻合口", "数据集",
        *[f"stomach__{feature}" for feature in features],
        *[f"tumor__{feature}" for feature in features],
    ]
    ws.append(header)
    mapping_ws = wb.create_sheet("病例映射与输入")
    mapping_ws.append(["住院号", "case_id", "吻合口", "数据集", "wide文件", "long文件"])
    for row in rows:
        case_features, stomach, tumor = read_patient_features(row["case_id"], features)
        if case_features != features:
            raise ValueError("Feature schema changed during writing")
        wide_path, long_path = source_paths(row["case_id"])
        ws.append([row["hospital_id"], row["case_id"], row["outcome"], row["split"], *stomach, *tumor])
        mapping_ws.append([row["hospital_id"], row["case_id"], row["outcome"], row["split"], str(wide_path), str(long_path)])
    wb.save(path)
    return header


def write_qc(path, train, validation, train_header, validation_header, batch_summary):
    if path.exists():
        raise FileExistsError(f"Refusing to overwrite: {path}")
    if train_header != validation_header or len(train_header) != 4 + 2 * ROI_FEATURE_COUNT:
        raise ValueError("Training/validation feature schemas do not agree")
    if {row["case_id"] for row in train} & {row["case_id"] for row in validation}:
        raise ValueError("Fixed split contains overlapping case IDs")
    wb = Workbook()
    ws = wb.active
    ws.title = "拆分与质控"
    rows = [
        ("影像组学协议", batch_summary["protocol"]),
        ("总患者数", 820),
        ("训练集", "574例：530 non-leak，44 leak"),
        ("验证集", "246例：227 non-leak，19 leak"),
        ("训练验证患者交叉", 0),
        ("每ROI特征数", ROI_FEATURE_COUNT),
        ("每患者特征数", 2 * ROI_FEATURE_COUNT),
        ("ROI命名", "stomach__ 和 tumor__ 前缀；同名基础特征不会覆盖"),
        ("结局使用边界", "预处理及提取阶段未读取结局；本表仅在固定拆分后关联结局。"),
        ("提取完成", f"{batch_summary['completed']}/{batch_summary['cases']}，失败 {batch_summary['failed']}"),
        ("PyRadiomics配置SHA256", batch_summary["config_sha256"]),
        ("提取器SHA256", batch_summary["extractor_sha256"]),
        ("下一步", "仅在训练集进行方差/相关性筛选及10折LASSO；验证集全程封存。"),
        ("限制", "v5技术负对照AUC低于预设阈值，但不等同于已完全排除残余批次效应；建模前后均需继续审计。"),
    ]
    ws.append(["项目", "结果"])
    for row in rows:
        ws.append(row)
    style_header(ws)
    ws.column_dimensions["A"].width = 28
    ws.column_dimensions["B"].width = 120
    mapping_ws = wb.create_sheet("病例映射")
    mapping_ws.append(["住院号", "case_id", "吻合口", "数据集", "wide文件", "long文件"])
    for row in train + validation:
        wide_path, long_path = source_paths(row["case_id"])
        mapping_ws.append([row["hospital_id"], row["case_id"], row["outcome"], row["split"], str(wide_path), str(long_path)])
    style_header(mapping_ws)
    mapping_ws.freeze_panes = "A2"
    for column in range(1, 7):
        mapping_ws.column_dimensions[get_column_letter(column)].width = 25 if column < 5 else 105
    wb.save(path)


def readback(path, expected_rows, expected_columns):
    wb = openpyxl.load_workbook(path, read_only=True, data_only=True)
    try:
        ws = wb["影像组学特征"]
        header = list(next(ws.values))
        rows = sum(1 for _ in ws.iter_rows(min_row=2, values_only=True))
        if rows != expected_rows or len(header) != expected_columns:
            raise ValueError(f"Read-back mismatch for {path}: rows={rows}, columns={len(header)}")
        if not header[4].startswith("stomach__"):
            raise ValueError(f"Missing stomach prefix in {path}")
        if not header[-1].startswith("tumor__"):
            raise ValueError(f"Missing tumor prefix in {path}")
    finally:
        wb.close()


def main():
    output_dir = DEFAULT_OUTPUT
    output_dir.mkdir(parents=True, exist_ok=True)
    if not BATCH_SUMMARY.is_file():
        raise FileNotFoundError(BATCH_SUMMARY)
    batch_summary = json.loads(BATCH_SUMMARY.read_text(encoding="utf-8"))
    if (batch_summary.get("cases"), batch_summary.get("completed"), batch_summary.get("failed"), batch_summary.get("features_per_roi")) != (820, 820, 0, ROI_FEATURE_COUNT):
        raise ValueError(f"Invalid batch summary: {batch_summary}")
    rows = read_mapping()
    train = [row for row in rows if row["split"] == "训练集"]
    validation = [row for row in rows if row["split"] == "验证集"]
    train_path = output_dir / "影像组学v5_训练集_70%_574例_3190特征.xlsx"
    validation_path = output_dir / "影像组学v5_验证集_30%_246例_3190特征.xlsx"
    qc_path = output_dir / "影像组学v5_训练验证集_拆分与质控.xlsx"
    train_header = write_matrix(train_path, train)
    validation_header = write_matrix(validation_path, validation)
    write_qc(qc_path, train, validation, train_header, validation_header, batch_summary)
    readback(train_path, 574, 3194)
    readback(validation_path, 246, 3194)
    mapping_path = output_dir / "影像组学v5_训练验证集_病例映射.csv"
    with mapping_path.open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.writer(handle)
        writer.writerow(["hospital_id", "case_id", "outcome", "split", "wide_path", "long_path"])
        for row in rows:
            wide_path, long_path = source_paths(row["case_id"])
            writer.writerow([row["hospital_id"], row["case_id"], row["outcome"], row["split"], wide_path, long_path])
    shutil.copy2(Path(__file__).resolve(), output_dir / Path(__file__).name)
    archive = PROJECT / "代码与复现说明/scripts"
    archive.mkdir(parents=True, exist_ok=True)
    shutil.copy2(Path(__file__).resolve(), archive / Path(__file__).name)
    print(json.dumps({
        "training": str(train_path), "validation": str(validation_path, ), "qc": str(qc_path),
        "mapping": str(mapping_path), "train_rows": len(train), "validation_rows": len(validation),
        "columns": len(train_header),
    }, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
