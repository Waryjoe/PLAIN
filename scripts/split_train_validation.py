#!/usr/bin/env python3
from plain_paths import plain_path
import argparse
import re
from pathlib import Path

import numpy as np
import openpyxl
from openpyxl import Workbook
from openpyxl.styles import Alignment, Font, PatternFill
from openpyxl.utils import get_column_letter


DEFAULT_INPUT = Path(
    plain_path("PROJECT", "汇总/所有数据_匹配CT影像组学患者_排除M1_最终820例.xlsx")
)
DEFAULT_RADIOMICS_ROOT = Path(
    plain_path("PROJECT", "CT2-EXCEL")
)
DEFAULT_OUTPUT_DIR = Path(
    plain_path("PROJECT", "汇总/训练验证集")
)
DEFAULT_SEED = 20260729


def normalize_hospital_id(value):
    text = str(value or "").strip()
    if text.endswith(".0"):
        text = text[:-2]
    digits = "".join(re.findall(r"\d+", text))
    return digits.zfill(8) if digits else None


def style_header(ws, row=1):
    fill = PatternFill("solid", fgColor="1F4E78")
    font = Font(color="FFFFFF", bold=True)
    for cell in ws[row]:
        cell.fill = fill
        cell.font = font
        cell.alignment = Alignment(horizontal="center", vertical="center")


def create_split_workbook(
    output_path,
    split_name,
    split_rows,
    headers,
    header_index,
    case_id_by_hospital_id,
    radiomics_root,
    source_path,
    seed,
    train_fraction,
    source_counts,
):
    workbook = Workbook()
    data_ws = workbook.active
    data_ws.title = "临床数据"
    data_ws.append(headers)
    for row in split_rows:
        data_ws.append(row)
    style_header(data_ws)
    data_ws.freeze_panes = "A2"
    data_ws.auto_filter.ref = (
        f"A1:{get_column_letter(len(headers))}{len(split_rows) + 1}"
    )
    data_ws.column_dimensions["A"].width = 13
    data_ws.column_dimensions["B"].width = 4
    for column in range(3, len(headers) + 1):
        data_ws.column_dimensions[get_column_letter(column)].width = 15
    for cell in data_ws["A"]:
        cell.number_format = "@"

    mapping_ws = workbook.create_sheet("病例映射")
    mapping_ws.append(
        [
            "住院号",
            "影像组学case_id",
            "吻合口",
            "数据集",
            "wide文件",
            "long文件",
        ]
    )
    id_index = header_index["住院号"]
    outcome_index = header_index["吻合口"]
    for row in split_rows:
        hospital_id = normalize_hospital_id(row[id_index])
        case_id = case_id_by_hospital_id[hospital_id]
        outcome = int(float(row[outcome_index]))
        if outcome == 0:
            wide_path = (
                radiomics_root
                / "非漏"
                / "wide"
                / f"{case_id}_radiomics_features_expanded_wide.xlsx"
            )
            long_path = (
                radiomics_root
                / "非漏"
                / "long"
                / f"{case_id}_radiomics_features_expanded_long.xlsx"
            )
        else:
            wide_path = (
                radiomics_root
                / "wide"
                / f"{case_id}_radiomics_features_expanded_wide.xlsx"
            )
            long_path = (
                radiomics_root
                / "long"
                / f"{case_id}_radiomics_features_expanded_long.xlsx"
            )
        if not wide_path.is_file() or not long_path.is_file():
            raise FileNotFoundError(
                f"Missing radiomics files for {hospital_id} / {case_id}"
            )
        mapping_ws.append(
            [
                hospital_id,
                case_id,
                outcome,
                split_name,
                str(wide_path),
                str(long_path),
            ]
        )
    style_header(mapping_ws)
    mapping_ws.freeze_panes = "A2"
    mapping_ws.column_dimensions["A"].width = 13
    mapping_ws.column_dimensions["B"].width = 26
    mapping_ws.column_dimensions["D"].width = 14
    mapping_ws.column_dimensions["E"].width = 100
    mapping_ws.column_dimensions["F"].width = 100
    for cell in mapping_ws["A"]:
        cell.number_format = "@"

    outcomes = [int(float(row[outcome_index])) for row in split_rows]
    info_ws = workbook.create_sheet("分组说明")
    info_rows = [
        ("源文件", str(source_path)),
        ("数据集", split_name),
        ("患者数", len(split_rows)),
        ("非漏组", outcomes.count(0)),
        ("漏组", outcomes.count(1)),
        ("划分方法", "患者级、按吻合口漏结局分层随机划分"),
        ("训练集比例", train_fraction),
        ("验证集比例", 1 - train_fraction),
        ("随机种子", seed),
        ("排序规则", "各结局层内先按标准化住院号排序，再随机抽样"),
        ("源队列总数", source_counts["total"]),
        ("源队列非漏组", source_counts["nonleak"]),
        ("源队列漏组", source_counts["leak"]),
        ("验证类型", "内部独立留出验证；不是外部验证"),
        ("M1患者", "PRIVATE_CASE_002已在划分前排除"),
        ("影像组学协议", "CT2-PyRadiomics-v1；每个ROI 1595个特征"),
    ]
    info_ws.append(["项目", "内容"])
    for item, value in info_rows:
        info_ws.append([item, value])
    style_header(info_ws)
    info_ws.column_dimensions["A"].width = 28
    info_ws.column_dimensions["B"].width = 100

    output_path.parent.mkdir(parents=True, exist_ok=True)
    if output_path.exists():
        raise FileExistsError(f"Refusing to overwrite existing output: {output_path}")
    workbook.save(output_path)


def main():
    parser = argparse.ArgumentParser(
        description="Create a reproducible stratified 70/30 train-validation split"
    )
    parser.add_argument("--input", default=str(DEFAULT_INPUT))
    parser.add_argument("--radiomics-root", default=str(DEFAULT_RADIOMICS_ROOT))
    parser.add_argument("--output-dir", default=str(DEFAULT_OUTPUT_DIR))
    parser.add_argument("--seed", type=int, default=DEFAULT_SEED)
    parser.add_argument("--train-fraction", type=float, default=0.7)
    args = parser.parse_args()

    input_path = Path(args.input).expanduser().resolve()
    radiomics_root = Path(args.radiomics_root).expanduser().resolve()
    output_dir = Path(args.output_dir).expanduser().resolve()
    if not 0 < args.train_fraction < 1:
        raise ValueError("train-fraction must be between 0 and 1")

    source_wb = openpyxl.load_workbook(
        input_path, read_only=True, data_only=False
    )
    source_ws = source_wb["匹配后数据"]
    source_rows = source_ws.iter_rows(
        min_row=1,
        max_row=source_ws.max_row,
        min_col=1,
        max_col=137,
        values_only=True,
    )
    headers = list(next(source_rows))
    header_index = {
        value: index for index, value in enumerate(headers) if value is not None
    }
    id_index = header_index["住院号"]
    outcome_index = header_index["吻合口"]
    rows = [list(row) for row in source_rows]

    leak_mapping_ws = source_wb["漏组精确匹配63例"]
    leak_case_id_by_hospital_id = {
        normalize_hospital_id(row[0]): row[1]
        for row in leak_mapping_ws.iter_rows(min_row=2, values_only=True)
    }
    source_wb.close()

    ids = []
    outcomes = []
    row_by_id = {}
    for row in rows:
        hospital_id = normalize_hospital_id(row[id_index])
        if hospital_id is None or hospital_id in row_by_id:
            raise ValueError(f"Invalid or duplicate hospital ID: {hospital_id}")
        row[id_index] = hospital_id
        outcome = int(float(row[outcome_index]))
        if outcome not in (0, 1):
            raise ValueError(f"Invalid outcome for {hospital_id}: {outcome}")
        ids.append(hospital_id)
        outcomes.append(outcome)
        row_by_id[hospital_id] = row

    source_counts = {
        "total": len(ids),
        "nonleak": outcomes.count(0),
        "leak": outcomes.count(1),
    }
    if source_counts != {"total": 820, "nonleak": 757, "leak": 63}:
        raise ValueError(f"Unexpected source cohort: {source_counts}")

    rng = np.random.default_rng(args.seed)
    train_ids = set()
    validation_ids = set()
    for outcome in (0, 1):
        stratum_ids = sorted(
            hospital_id
            for hospital_id, value in zip(ids, outcomes)
            if value == outcome
        )
        shuffled_ids = list(rng.permutation(stratum_ids))
        train_count = int(round(len(stratum_ids) * args.train_fraction))
        train_ids.update(shuffled_ids[:train_count])
        validation_ids.update(shuffled_ids[train_count:])

    if train_ids & validation_ids or train_ids | validation_ids != set(ids):
        raise ValueError("Train-validation split is not disjoint and exhaustive")

    train_rows = [row_by_id[hospital_id] for hospital_id in ids if hospital_id in train_ids]
    validation_rows = [
        row_by_id[hospital_id]
        for hospital_id in ids
        if hospital_id in validation_ids
    ]
    train_outcomes = [int(float(row[outcome_index])) for row in train_rows]
    validation_outcomes = [
        int(float(row[outcome_index])) for row in validation_rows
    ]
    if (
        len(train_rows),
        train_outcomes.count(0),
        train_outcomes.count(1),
        len(validation_rows),
        validation_outcomes.count(0),
        validation_outcomes.count(1),
    ) != (574, 530, 44, 246, 227, 19):
        raise ValueError("Unexpected split counts")

    case_id_by_hospital_id = {
        hospital_id: f"gc_ZJCH_{hospital_id}" for hospital_id in ids
    }
    case_id_by_hospital_id.update(leak_case_id_by_hospital_id)
    if len(case_id_by_hospital_id) != 820:
        raise ValueError("Radiomics case-ID mapping is incomplete")

    train_path = output_dir / "训练集_70%_574例.xlsx"
    validation_path = output_dir / "验证集_30%_246例.xlsx"
    create_split_workbook(
        train_path,
        "训练集",
        train_rows,
        headers,
        header_index,
        case_id_by_hospital_id,
        radiomics_root,
        input_path,
        args.seed,
        args.train_fraction,
        source_counts,
    )
    create_split_workbook(
        validation_path,
        "验证集",
        validation_rows,
        headers,
        header_index,
        case_id_by_hospital_id,
        radiomics_root,
        input_path,
        args.seed,
        args.train_fraction,
        source_counts,
    )

    check_sets = {}
    for label, path, expected in (
        ("train", train_path, (574, 530, 44)),
        ("validation", validation_path, (246, 227, 19)),
    ):
        check_wb = openpyxl.load_workbook(path, read_only=True, data_only=True)
        check_ws = check_wb["临床数据"]
        check_rows = check_ws.iter_rows(
            min_row=1, max_col=137, values_only=True
        )
        check_headers = list(next(check_rows))
        check_index = {
            value: index
            for index, value in enumerate(check_headers)
            if value is not None
        }
        check_ids = set()
        check_outcomes = []
        for row in check_rows:
            check_ids.add(normalize_hospital_id(row[check_index["住院号"]]))
            check_outcomes.append(int(float(row[check_index["吻合口"]])))
        observed = (
            len(check_ids),
            check_outcomes.count(0),
            check_outcomes.count(1),
        )
        if observed != expected:
            raise ValueError(
                f"{label} read-back counts changed: {observed}; expected {expected}"
            )
        if check_wb["病例映射"].max_row != expected[0] + 1:
            raise ValueError(f"{label} radiomics mapping row count changed")
        check_sets[label] = check_ids
        check_wb.close()

    if (
        check_sets["train"] & check_sets["validation"]
        or check_sets["train"] | check_sets["validation"] != set(ids)
    ):
        raise ValueError("Read-back split integrity check failed")

    print(
        {
            "input": str(input_path),
            "seed": args.seed,
            "method": "patient-level outcome-stratified random split",
            "train": {
                "path": str(train_path),
                "total": 574,
                "nonleak": 530,
                "leak": 44,
            },
            "validation": {
                "path": str(validation_path),
                "total": 246,
                "nonleak": 227,
                "leak": 19,
            },
            "disjoint": True,
            "exhaustive": True,
            "radiomics_files_verified": True,
            "readback_valid": True,
        }
    )


if __name__ == "__main__":
    main()
