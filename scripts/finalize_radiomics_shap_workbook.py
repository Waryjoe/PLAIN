"""Rebuild the numerical workbook with standard OOXML and verify round trips."""

import json
import math
import sys
from pathlib import Path

import openpyxl
from lxml import etree
from openpyxl.styles import Alignment, Font, PatternFill
from openpyxl.utils import get_column_letter


def main(folder):
    folder = Path(folder)
    path = folder / "Radiomics_SHAP_SourceData_and_QC.xlsx"
    source = openpyxl.load_workbook(path, read_only=True, data_only=True)
    tables = {}
    for sheet in source:
        # Some openxlsx versions leave A1 dimensions and unused drawing links.
        sheet.reset_dimensions()
        tables[sheet.title] = list(sheet.values)
    source.close()
    workbook = openpyxl.Workbook()
    workbook.remove(workbook.active)
    for name, rows in tables.items():
        sheet = workbook.create_sheet(name)
        for row in rows:
            sheet.append(row)
        sheet.freeze_panes = "A2"
        sheet.auto_filter.ref = sheet.dimensions
        for cell in sheet[1]:
            cell.font = Font(name="Arial", bold=True, color="FFFFFF", size=10)
            cell.fill = PatternFill("solid", fgColor="24647A")
            cell.alignment = Alignment(vertical="center", wrap_text=True)
        sheet.row_dimensions[1].height = 30
        for column in range(1, sheet.max_column + 1):
            header = str(sheet.cell(1, column).value)
            width = 24 if len(header) < 27 else 32
            if header in {"Feature", "Model_column", "Feature_original"}:
                width = 76
            if header in {"Short_name", "Check", "Path"}:
                width = 52
            if header == "Detail":
                width = 112
            sheet.column_dimensions[get_column_letter(column)].width = width
        for row in sheet.iter_rows(min_row=2):
            for cell in row:
                cell.font = Font(name="Arial", size=10)
                cell.alignment = Alignment(vertical="top", wrap_text=name in {"Readme", "Input_provenance"})
                if isinstance(cell.value, float):
                    cell.number_format = "0.000000"
        if name == "Readme":
            for i in range(2, sheet.max_row + 1):
                sheet.row_dimensions[i].height = 45
    temp = path.with_name(path.stem + ".checked.xlsx")
    workbook.save(temp)
    check = openpyxl.load_workbook(temp, read_only=False, data_only=True)
    sizes = {}
    for name, expected in tables.items():
        actual = list(check[name].values)
        assert len(actual) == len(expected), name
        sizes[name] = [len(actual) - 1, len(actual[0])]
        for arow, erow in zip(actual, expected):
            for a, e in zip(arow, erow):
                if isinstance(e, (int, float)):
                    assert math.isclose(a, e, rel_tol=1e-13, abs_tol=1e-13), (name, a, e)
                else:
                    assert a == e, (name, a, e)
    assert sizes["Training_SHAP"][0] == 574
    assert sizes["Internal_SHAP"][0] == 246
    assert sizes["Plot_source_long"][0] == 5740
    assert all(row[2] == "PASS" for row in list(check["Verification"].values)[1:])
    errors = {}
    for sheet in ["Training_SHAP", "Internal_SHAP"]:
        errors[sheet] = max(abs(sum(row[4:]) - row[3]) for row in list(check[sheet].values)[1:])
        assert errors[sheet] < 1e-8
    check.close()
    temp.replace(path)
    svg_checks = {}
    for svg in folder.glob("*.svg"):
        doc = etree.parse(str(svg))
        counts = {name: len(doc.xpath(f"//*[local-name()='{name}']")) for name in ["circle", "text", "image"]}
        assert counts["circle"] == 5740 and counts["text"] > 10 and counts["image"] == 0
        svg_checks[svg.name] = counts
    result = {"status": "PASS", "sheet_dimensions": sizes, "SHAP_additivity_after_Excel_roundtrip": errors, "SVG_vector_checks": svg_checks}
    (folder / "Reproducibility/Delivery_QC.json").write_text(json.dumps(result, indent=2), encoding="utf-8")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main(sys.argv[1])
