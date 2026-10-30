#!/usr/bin/env python3
from plain_paths import plain_path
import argparse
import os
import subprocess
from pathlib import Path

import openpyxl
import pandas as pd
from docx import Document
from docx.enum.section import WD_ORIENT
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.enum.table import WD_CELL_VERTICAL_ALIGNMENT
from docx.oxml import OxmlElement
from docx.oxml.ns import qn
from docx.shared import Inches, Pt
from openpyxl.styles import Alignment, Border, Font, Side
from openpyxl.utils import get_column_letter


DEFAULT_INPUT = Path(
    plain_path("PROJECT", "汇总/所有数据_匹配CT影像组学患者_排除M1_最终820例.xlsx")
)
DEFAULT_OUTPUT_DIR = Path(
    plain_path("PROJECT", "汇总/三线表")
)
DEFAULT_R_SCRIPT = Path(__file__).with_name("generate_table1.R")

CONTINUOUS = [
    "年龄",
    "BMI",
    "肿瘤最长直径",
    "肿瘤次长径",
    "肿瘤短径",
    "CEA",
    "CA199",
    "CA125",
    "AFP",
    "CA242",
    "CA724",
    "红细胞比容Hct",
    "血红蛋白Hb",
    "红细胞RBC",
    "淋巴细胞Lymph",
    "血小板计数PLT",
    "白细胞WBC",
    "嗜碱性粒细胞百分比",
    "嗜酸性粒细胞百分比",
    "中性粒细胞百分比",
    "淋巴细胞百分比",
    "中性粒细胞Neut",
    "钾K",
    "钠Na",
    "钙Ca",
    "氯Cl",
    "总蛋白TP",
    "白蛋白ALB",
    "球蛋白GLO",
    "白蛋白球蛋白比值",
    "C反应蛋白CRP",
    "前白蛋白PA",
    "总胆红素TBIL",
    "ALT",
    "直接胆红素DBIL",
    "间接胆红素IBIL",
    "AST",
    "碱性磷酸酶ALP",
    "乳酸脱氢酶LDH",
    "甘油三酯TG",
    "总胆固醇TC",
    "高密度脂蛋白",
    "低密度脂蛋白",
    "凝血酶原时间PT",
    "纤维蛋白原Fbg",
    "凝血酶时间TT",
]

CATEGORICAL_SOURCE = [
    "性别",
    "家族史",
    "肿瘤部位2",
    "分化程度2",
    "T1",
    "N1",
    "TG2",
]


def set_cell_text(cell, text, bold=False, size=8.5):
    cell.text = ""
    paragraph = cell.paragraphs[0]
    paragraph.alignment = (
        WD_ALIGN_PARAGRAPH.LEFT
        if cell._tc.getparent().index(cell._tc) == 0
        else WD_ALIGN_PARAGRAPH.CENTER
    )
    run = paragraph.add_run(str(text))
    run.bold = bold
    run.font.name = "Times New Roman"
    run.font.size = Pt(size)
    run._element.rPr.rFonts.set(qn("w:eastAsia"), "宋体")
    cell.vertical_alignment = WD_CELL_VERTICAL_ALIGNMENT.CENTER


def set_table_borders(table):
    tbl_pr = table._tbl.tblPr
    borders = tbl_pr.first_child_found_in("w:tblBorders")
    if borders is None:
        borders = OxmlElement("w:tblBorders")
        tbl_pr.append(borders)
    for edge in ("top", "left", "bottom", "right", "insideH", "insideV"):
        element = borders.find(qn(f"w:{edge}"))
        if element is None:
            element = OxmlElement(f"w:{edge}")
            borders.append(element)
        if edge in ("top", "bottom"):
            element.set(qn("w:val"), "single")
            element.set(qn("w:sz"), "8")
            element.set(qn("w:color"), "000000")
        else:
            element.set(qn("w:val"), "nil")


def set_cell_bottom_border(cell, size=6):
    tc_pr = cell._tc.get_or_add_tcPr()
    borders = tc_pr.first_child_found_in("w:tcBorders")
    if borders is None:
        borders = OxmlElement("w:tcBorders")
        tc_pr.append(borders)
    bottom = OxmlElement("w:bottom")
    bottom.set(qn("w:val"), "single")
    bottom.set(qn("w:sz"), str(size))
    bottom.set(qn("w:color"), "000000")
    borders.append(bottom)


def repeat_table_header(row):
    tr_pr = row._tr.get_or_add_trPr()
    tbl_header = OxmlElement("w:tblHeader")
    tbl_header.set(qn("w:val"), "true")
    tr_pr.append(tbl_header)


def build_docx(table_rows, output_path, total_n, no_n, leak_n):
    document = Document()
    section = document.sections[0]
    section.orientation = WD_ORIENT.LANDSCAPE
    section.page_width, section.page_height = section.page_height, section.page_width
    section.top_margin = Inches(0.45)
    section.bottom_margin = Inches(0.45)
    section.left_margin = Inches(0.45)
    section.right_margin = Inches(0.45)

    title = document.add_paragraph()
    title.alignment = WD_ALIGN_PARAGRAPH.CENTER
    run = title.add_run("表1  按吻合口漏分组的基线临床特征")
    run.bold = True
    run.font.name = "Times New Roman"
    run.font.size = Pt(11)
    run._element.rPr.rFonts.set(qn("w:eastAsia"), "宋体")

    columns = [
        "变量",
        f"总体 (N={total_n})",
        f"非漏组 (N={no_n})",
        f"漏组 (N={leak_n})",
        "χ²/Z",
        "P",
    ]
    table = document.add_table(rows=1, cols=len(columns))
    table.autofit = False
    widths = [Inches(3.25), Inches(1.55), Inches(1.55), Inches(1.55), Inches(0.75), Inches(0.7)]
    for index, (cell, label) in enumerate(zip(table.rows[0].cells, columns)):
        cell.width = widths[index]
        set_cell_text(cell, label, bold=True, size=8.5)
        set_cell_bottom_border(cell)
    repeat_table_header(table.rows[0])
    set_table_borders(table)

    for record in table_rows.to_dict("records"):
        cells = table.add_row().cells
        is_variable = record["Row_type"] == "variable"
        values = [
            record["Variable"],
            record["Total"],
            record["No_leak"],
            record["Leak"],
            record["Statistic"],
            record["P"],
        ]
        for index, (cell, value) in enumerate(zip(cells, values)):
            cell.width = widths[index]
            set_cell_text(cell, value, bold=(is_variable and index == 0), size=8.2)
        p_numeric = record.get("P_numeric")
        if is_variable and pd.notna(p_numeric) and float(p_numeric) < 0.05:
            for run in cells[5].paragraphs[0].runs:
                run.bold = True

    notes = [
        "注：连续变量以中位数（Q1，Q3）表示，采用双侧Mann–Whitney U检验，统计量列报告Z值；",
        "分类变量以n（%）表示，采用Pearson χ²检验；当期望频数<5时采用Fisher精确检验或Fisher–Freeman–Halton检验。",
        "百分比按该变量非缺失病例计算；P值为未校正的探索性组间比较。M1远处转移患者PRIVATE_CASE_002已排除。",
        "肿瘤部位编码：1=胃上部，2=胃中部，3=胃下部，4=全胃；分化程度编码：1=低危，2=中危，3=高危。",
        "TG2按原数据库0/1/2编码展示；非漏组全部为0，存在数据来源或标签代理泄漏风险，定义核实前不得用于建模。",
    ]
    for note in notes:
        paragraph = document.add_paragraph()
        paragraph.paragraph_format.space_after = Pt(0)
        run = paragraph.add_run(note)
        run.font.name = "Times New Roman"
        run.font.size = Pt(8)
        run._element.rPr.rFonts.set(qn("w:eastAsia"), "宋体")

    document.save(output_path)


def build_xlsx(table_rows, output_path, total_n, no_n, leak_n, metadata):
    workbook = openpyxl.Workbook()
    ws = workbook.active
    ws.title = "Table 1"
    ws.merge_cells("A1:F1")
    ws["A1"] = "表1  按吻合口漏分组的基线临床特征"
    ws["A1"].font = Font(name="宋体", size=12, bold=True)
    ws["A1"].alignment = Alignment(horizontal="center")

    columns = [
        "变量",
        f"总体 (N={total_n})",
        f"非漏组 (N={no_n})",
        f"漏组 (N={leak_n})",
        "χ²/Z",
        "P",
    ]
    for column, value in enumerate(columns, start=1):
        cell = ws.cell(2, column, value)
        cell.font = Font(name="宋体", size=10, bold=True)
        cell.alignment = Alignment(horizontal="center", vertical="center")

    thin = Side(style="thin", color="000000")
    medium = Side(style="medium", color="000000")
    for cell in ws[2]:
        cell.border = Border(top=medium, bottom=thin)

    start_row = 3
    for row_index, record in enumerate(table_rows.to_dict("records"), start=start_row):
        values = [
            record["Variable"],
            record["Total"],
            record["No_leak"],
            record["Leak"],
            record["Statistic"],
            record["P"],
        ]
        is_variable = record["Row_type"] == "variable"
        for column, value in enumerate(values, start=1):
            cell = ws.cell(row_index, column, value)
            cell.font = Font(
                name="宋体",
                size=9,
                bold=(is_variable and column == 1)
                or (
                    column == 6
                    and is_variable
                    and pd.notna(record["P_numeric"])
                    and float(record["P_numeric"]) < 0.05
                ),
            )
            cell.alignment = Alignment(
                horizontal="left" if column == 1 else "center",
                vertical="center",
            )

    last_data_row = start_row + len(table_rows) - 1
    for cell in ws[last_data_row]:
        cell.border = Border(bottom=medium)

    notes = [
        "注：连续变量以中位数（Q1，Q3）表示，采用双侧Mann–Whitney U检验，统计量列报告Z值；分类变量以n（%）表示，采用Pearson χ²检验；期望频数<5时采用Fisher精确检验或Fisher–Freeman–Halton检验。",
        "百分比按该变量非缺失病例计算；P值为未校正的探索性组间比较。M1远处转移患者PRIVATE_CASE_002已排除。",
        "肿瘤部位：1=胃上部，2=胃中部，3=胃下部，4=全胃；分化程度：1=低危，2=中危，3=高危。TG2按原数据库0/1/2编码展示；非漏组全部为0，定义核实前不得用于建模。",
    ]
    for offset, note in enumerate(notes, start=2):
        row = last_data_row + offset
        ws.merge_cells(start_row=row, start_column=1, end_row=row, end_column=6)
        ws.cell(row, 1, note)
        ws.cell(row, 1).font = Font(name="宋体", size=8)
        ws.cell(row, 1).alignment = Alignment(wrap_text=True, vertical="top")

    widths = [42, 20, 20, 20, 12, 11]
    for index, width in enumerate(widths, start=1):
        ws.column_dimensions[get_column_letter(index)].width = width
    ws.freeze_panes = "A3"
    ws.sheet_view.showGridLines = False

    detail = workbook.create_sheet("统计明细")
    detail.append(list(table_rows.columns))
    for record in table_rows.itertuples(index=False, name=None):
        detail.append(list(record))
    for cell in detail[1]:
        cell.font = Font(bold=True)
    detail.freeze_panes = "A2"

    info = workbook.create_sheet("统计说明")
    info.append(["项目", "内容"])
    for key, value in metadata:
        info.append([key, value])
    for cell in info[1]:
        cell.font = Font(bold=True)
    info.column_dimensions["A"].width = 30
    info.column_dimensions["B"].width = 95

    exclusion = workbook.create_sheet("排除记录")
    exclusion.append(["住院号", "影像组学case_id", "姓名", "M", "排除原因"])
    exclusion.append(
        ["PRIVATE_RECORD_001", "PRIVATE_CASE_002", "PRIVATE_NAME_001", "M1", "明确远处转移，排除"]
    )
    for cell in exclusion[1]:
        cell.font = Font(bold=True)
    exclusion.column_dimensions["A"].width = 14
    exclusion.column_dimensions["B"].width = 22
    exclusion.column_dimensions["E"].width = 36

    workbook.save(output_path)


def main():
    parser = argparse.ArgumentParser(description="Build publication-style Table 1")
    parser.add_argument("--input", default=str(DEFAULT_INPUT))
    parser.add_argument("--output-dir", default=str(DEFAULT_OUTPUT_DIR))
    parser.add_argument("--r-script", default=str(DEFAULT_R_SCRIPT))
    args = parser.parse_args()

    input_path = Path(args.input).expanduser().resolve()
    output_dir = Path(args.output_dir).expanduser().resolve()
    r_script = Path(args.r_script).expanduser().resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    workbook = openpyxl.load_workbook(input_path, read_only=True, data_only=True)
    ws = workbook["匹配后数据"]
    rows = ws.iter_rows(min_row=1, max_col=137, values_only=True)
    headers = list(next(rows))
    header_index = {value: index for index, value in enumerate(headers) if value}
    required = ["吻合口", "M", *CONTINUOUS, *CATEGORICAL_SOURCE]
    missing = sorted(set(required) - set(header_index))
    if missing:
        raise ValueError(f"Required columns are missing: {missing}")

    records = []
    for row in rows:
        record = {column: row[header_index[column]] for column in required}
        m_value = str(record["M"] or "").strip().upper()
        if m_value == "M1" or m_value == "1":
            raise ValueError("M1 patient remains in cleaned clinical workbook")
        records.append(record)
    workbook.close()

    frame = pd.DataFrame(records)
    frame["吻合口"] = pd.to_numeric(frame["吻合口"], errors="raise").astype(int)
    total_n = len(frame)
    no_n = int((frame["吻合口"] == 0).sum())
    leak_n = int((frame["吻合口"] == 1).sum())
    if (total_n, no_n, leak_n) != (820, 757, 63):
        raise ValueError(
            f"Unexpected cohort counts: total={total_n}, no_leak={no_n}, leak={leak_n}"
        )

    source_csv = output_dir / "Table1_source_data_deidentified_820.csv"
    statistics_csv = output_dir / "Table1_statistics_820.csv"
    xlsx_path = output_dir / "Table1_临床特征_排除M1_820例_三线表.xlsx"
    docx_path = output_dir / "Table1_临床特征_排除M1_820例_三线表.docx"
    frame.drop(columns=["M"]).to_csv(source_csv, index=False, encoding="utf-8")

    r_environment = os.environ.copy()
    r_environment.update(
        {
            "LANG": "zh_CN.UTF-8",
            "LC_ALL": "zh_CN.UTF-8",
            "LC_CTYPE": "zh_CN.UTF-8",
        }
    )
    subprocess.run(
        ["Rscript", str(r_script), str(source_csv), str(statistics_csv)],
        check=True,
        env=r_environment,
    )
    table_rows = pd.read_csv(statistics_csv, encoding="utf-8", keep_default_na=False)
    table_rows["P_numeric"] = pd.to_numeric(
        table_rows["P_numeric"], errors="coerce"
    )
    if table_rows.empty or table_rows["P_numeric"].notna().sum() < 50:
        raise ValueError("R Table 1 output is unexpectedly incomplete")

    metadata = [
        ("临床数据源", str(input_path)),
        ("分析单位", "患者"),
        ("总样本量", total_n),
        ("非漏组", no_n),
        ("漏组", leak_n),
        ("明确M1排除", "住院号PRIVATE_RECORD_001 / PRIVATE_CASE_002 / PRIVATE_NAME_001"),
        ("连续变量", "中位数（Q1，Q3）；双侧Mann–Whitney U检验"),
        ("分类变量", "n（%）；Pearson χ²或Fisher类精确检验"),
        ("缺失值", "不参与该变量检验，单独报告"),
        ("P值", "双侧、未作多重比较校正，仅用于描述性比较"),
        (
            "TG2风险",
            "非漏组757例全部为0，可能是数据来源或结局代理变量；定义核实前不得用于建模",
        ),
        ("统计软件", "R 4.6.0（base stats）"),
    ]
    build_xlsx(table_rows, xlsx_path, total_n, no_n, leak_n, metadata)
    build_docx(table_rows, docx_path, total_n, no_n, leak_n)

    check_wb = openpyxl.load_workbook(xlsx_path, read_only=True, data_only=True)
    if "Table 1" not in check_wb.sheetnames or check_wb["Table 1"].max_row < 75:
        raise ValueError("Excel Table 1 read-back validation failed")
    check_wb.close()
    check_doc = Document(docx_path)
    if len(check_doc.tables) != 1 or len(check_doc.tables[0].rows) < 75:
        raise ValueError("Word Table 1 read-back validation failed")

    print(
        {
            "input": str(input_path),
            "cohort": {"total": total_n, "no_leak": no_n, "leak": leak_n},
            "table_rows": len(table_rows),
            "xlsx": str(xlsx_path),
            "docx": str(docx_path),
            "statistics_csv": str(statistics_csv),
            "source_csv": str(source_csv),
            "readback_valid": True,
        }
    )


if __name__ == "__main__":
    main()
