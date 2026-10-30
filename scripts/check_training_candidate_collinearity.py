#!/usr/bin/env python3
"""Check correlation and multicollinearity among priority clinical candidates.

This training-only diagnostic follows the preceding univariable screen. It does
not fit an outcome model and never reads the validation set.
"""
from plain_paths import plain_path

import argparse
import hashlib
import json
import sys
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd
from openpyxl import load_workbook
from openpyxl.styles import Alignment, Font, PatternFill
from openpyxl.utils import get_column_letter


DEFAULT_TRAINING_INPUT = Path(
    plain_path("PROJECT", "汇总/训练验证集/训练集_70%_574例.xlsx")
)
DEFAULT_CANDIDATE_INPUT = Path(
    plain_path("PROJECT", "汇总/训练集单因素候选变量筛选_排除肿瘤部位/训练集_单因素候选变量筛选_排除肿瘤部位.xlsx")
)
DEFAULT_OUTPUT_DIR = Path(
    plain_path("PROJECT", "汇总/训练集优先候选变量_相关性与VIF_排除肿瘤部位")
)
EXCLUDED_VARIABLES = {"肿瘤部位2", "TG2"}
PRIORITY_TIER = "优先候选：P<0.05"


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def vif_interpretation(value):
    if not np.isfinite(value):
        return "无法估计"
    if value < 2:
        return "低共线性（VIF<2）"
    if value < 5:
        return "轻度至中度共线性（2≤VIF<5）"
    if value < 10:
        return "需关注（5≤VIF<10）"
    return "高共线性（VIF≥10）"


def correlation_interpretation(rho):
    value = abs(rho)
    if value < 0.5:
        return "低相关（|rho|<0.50）"
    if value < 0.7:
        return "中等相关（0.50≤|rho|<0.70）"
    if value < 0.8:
        return "较高相关（0.70≤|rho|<0.80）"
    return "高度相关（|rho|≥0.80）"


def spearman_coefficient(left, right):
    """Compute Spearman rho as Pearson correlation of average ranks."""
    left_rank = pd.Series(left).rank(method="average").to_numpy(dtype=float)
    right_rank = pd.Series(right).rank(method="average").to_numpy(dtype=float)
    if left_rank.std(ddof=0) == 0 or right_rank.std(ddof=0) == 0:
        return np.nan
    return float(np.corrcoef(left_rank, right_rank)[0, 1])


def chi_square_statistic(table):
    observed = table.to_numpy(dtype=float)
    total = observed.sum()
    if total <= 0:
        return np.nan
    expected = np.outer(observed.sum(axis=1), observed.sum(axis=0)) / total
    if (expected <= 0).any():
        return np.nan
    return float(np.sum((observed - expected) ** 2 / expected))


def choose_priority_candidates(path):
    table = pd.read_excel(path, sheet_name="候选变量清单", engine="openpyxl")
    required = {"变量", "变量类型", "候选层级"}
    missing = required - set(table.columns)
    if missing:
        raise ValueError(f"候选工作表缺少列: {sorted(missing)}")
    result = table.loc[
        table["候选层级"].eq(PRIORITY_TIER)
        & ~table["变量"].isin(EXCLUDED_VARIABLES),
        [column for column in table.columns if column in {"变量", "变量类型", "角色", "肿瘤相关", "P值", "BH_FDR", "候选层级"}],
    ].copy()
    if result.empty:
        raise ValueError("未找到P<0.05的优先候选变量")
    return result.sort_values(["P值", "变量"], na_position="last").reset_index(drop=True)


def clean_candidate_data(clinical, continuous, categorical):
    cleaned = clinical.copy()
    zero_counts = []
    for variable in continuous:
        values = pd.to_numeric(cleaned[variable], errors="coerce")
        zero_count = int(values.eq(0).sum())
        values.loc[values.eq(0)] = np.nan
        cleaned[variable] = values
        zero_counts.append(
            {
                "变量": variable,
                "变量类型": "连续",
                "精确0按缺失处理数": zero_count,
                "清洗后缺失数": int(values.isna().sum()),
                "清洗后缺失率": float(values.isna().mean()),
            }
        )
    for variable in categorical:
        values = cleaned[variable]
        zero_counts.append(
            {
                "变量": variable,
                "变量类型": "分类",
                "精确0按缺失处理数": 0,
                "清洗后缺失数": int(values.isna().sum()),
                "清洗后缺失率": float(values.isna().mean()),
            }
        )
    return cleaned, pd.DataFrame(zero_counts)


def pairwise_spearman(cleaned, continuous):
    rows = []
    for index, left in enumerate(continuous):
        for right in continuous[index + 1 :]:
            values = cleaned[[left, right]].dropna()
            if len(values) < 3:
                rho = np.nan
            else:
                rho = spearman_coefficient(values[left], values[right])
            rows.append(
                {
                    "变量1": left,
                    "变量2": right,
                    "配对完整案例数": int(len(values)),
                    "Spearman_rho": rho,
                    "相关性P值": np.nan,
                    "相关性判断": correlation_interpretation(rho) if np.isfinite(rho) else "无法估计",
                    "需要处理": bool(np.isfinite(rho) and abs(rho) >= 0.70),
                }
            )
    return pd.DataFrame(rows).sort_values("Spearman_rho", key=lambda s: s.abs(), ascending=False)


def categorical_association(cleaned, categorical):
    rows = []
    for index, left in enumerate(categorical):
        for right in categorical[index + 1 :]:
            values = cleaned[[left, right]].dropna()
            table = pd.crosstab(values[left], values[right])
            if table.shape[0] < 2 or table.shape[1] < 2:
                chi2 = cramers_v = np.nan
            else:
                chi2 = chi_square_statistic(table)
                denom = len(values) * min(table.shape[0] - 1, table.shape[1] - 1)
                cramers_v = float(np.sqrt(chi2 / denom)) if denom > 0 else np.nan
            rows.append(
                {
                    "变量1": left,
                    "变量2": right,
                    "完整案例数": int(len(values)),
                    "卡方统计量": chi2,
                    "Cramers_V": cramers_v,
                }
            )
    return pd.DataFrame(rows)


def make_vif_design(cleaned, continuous, categorical):
    pieces = []
    imputation_rows = []
    term_map = {}
    for variable in continuous:
        values = pd.to_numeric(cleaned[variable], errors="coerce")
        fill_value = float(values.median())
        pieces.append(pd.DataFrame({variable: values.fillna(fill_value)}))
        term_map[variable] = [variable]
        imputation_rows.append(
            {
                "变量": variable,
                "变量类型": "连续",
                "缺失填补方式": "训练集内中位数填补（仅用于VIF设计矩阵）",
                "填补值或参考水平": fill_value,
                "填补数": int(values.isna().sum()),
            }
        )
    for variable in categorical:
        raw = cleaned[variable].copy()
        if raw.dropna().empty:
            raise ValueError(f"{variable} 全部缺失，无法构建VIF")
        mode = raw.dropna().mode().iloc[0]
        filled = raw.fillna(mode).astype(str)
        levels = sorted(filled.unique(), key=lambda x: float(x))
        reference = levels[0]
        dummies = pd.get_dummies(
            pd.Categorical(filled, categories=levels),
            prefix=variable,
            prefix_sep="=",
            drop_first=True,
            dtype=float,
        )
        pieces.append(dummies)
        term_map[variable] = list(dummies.columns)
        imputation_rows.append(
            {
                "变量": variable,
                "变量类型": "分类",
                "缺失填补方式": "训练集内众数填补（仅用于VIF设计矩阵）",
                "填补值或参考水平": f"填补众数={mode}; VIF参考水平={reference}",
                "填补数": int(raw.isna().sum()),
            }
        )
    design = pd.concat(pieces, axis=1).astype(float)
    if design.isna().any().any() or not np.isfinite(design.to_numpy()).all():
        raise ValueError("VIF设计矩阵仍含缺失或非有限值")
    return design, term_map, pd.DataFrame(imputation_rows)


def design_vif(design):
    rows = []
    values = design.to_numpy(dtype=float)
    for index, column in enumerate(design.columns):
        value = calculate_vif(values, index)
        rows.append(
            {
                "设计矩阵列": column,
                "VIF": value,
                "VIF判断": vif_interpretation(value),
            }
        )
    return pd.DataFrame(rows).sort_values("VIF", ascending=False)


def determinant(matrix):
    sign, logdet = np.linalg.slogdet(matrix)
    if sign <= 0:
        return np.nan
    return float(np.exp(logdet))


def calculate_vif(values, index):
    """Calculate classical VIF from R-squared without statsmodels."""
    target = values[:, index]
    others = np.delete(values, index, axis=1)
    design = np.column_stack([np.ones(len(target)), others])
    coefficients, _, _, _ = np.linalg.lstsq(design, target, rcond=None)
    fitted = design @ coefficients
    total_sum_squares = float(np.sum((target - target.mean()) ** 2))
    residual_sum_squares = float(np.sum((target - fitted) ** 2))
    if total_sum_squares <= 0:
        return np.nan
    r_squared = 1 - residual_sum_squares / total_sum_squares
    return float(1 / max(1 - r_squared, np.finfo(float).eps))


def style_workbook(path):
    workbook = load_workbook(path)
    header_fill = PatternFill("solid", fgColor="1F4E78")
    for worksheet in workbook.worksheets:
        worksheet.freeze_panes = "A2"
        worksheet.auto_filter.ref = worksheet.dimensions
        for cell in worksheet[1]:
            cell.font = Font(bold=True, color="FFFFFF")
            cell.fill = header_fill
            cell.alignment = Alignment(horizontal="center", vertical="center", wrap_text=True)
        for column_cells in worksheet.columns:
            letter = get_column_letter(column_cells[0].column)
            length = max(len(str(cell.value or "")) for cell in column_cells)
            worksheet.column_dimensions[letter].width = min(max(length + 2, 12), 40)
    workbook.save(path)


def termwise_gvif(design, term_map):
    correlation = np.corrcoef(design.to_numpy(), rowvar=False)
    all_indices = list(range(design.shape[1]))
    det_full = determinant(correlation)
    rows = []
    for term, columns in term_map.items():
        indices = [design.columns.get_loc(column) for column in columns]
        others = [index for index in all_indices if index not in indices]
        det_term = determinant(correlation[np.ix_(indices, indices)])
        det_others = determinant(correlation[np.ix_(others, others)]) if others else 1.0
        if not np.isfinite(det_full) or det_full <= 0:
            gvif = np.nan
        else:
            gvif = float(det_term * det_others / det_full)
        df = len(indices)
        adjusted = float(gvif ** (1 / (2 * df))) if np.isfinite(gvif) and gvif >= 0 else np.nan
        rows.append(
            {
                "变量": term,
                "自由度_df": df,
                "设计矩阵列": "；".join(columns),
                "GVIF": gvif,
                "GVIF^(1/(2df))": adjusted,
                "可比VIF": gvif if df == 1 else np.nan,
                "判断": vif_interpretation(gvif) if df == 1 else (
                    "低共线性（调整GVIF<1.4）" if adjusted < 1.4 else
                    "需关注（调整GVIF≥1.4）" if adjusted < 2 else
                    "高共线性（调整GVIF≥2）"
                ),
            }
        )
    return pd.DataFrame(rows).sort_values("GVIF^(1/(2df))", ascending=False)


def write_report(path, candidates, correlations, gvif):
    high_pairs = correlations.loc[correlations["需要处理"], ["变量1", "变量2", "Spearman_rho"]]
    high_vif = gvif.loc[
        (gvif["自由度_df"].eq(1) & gvif["GVIF"].ge(5))
        | (gvif["自由度_df"].gt(1) & gvif["GVIF^(1/(2df))"].ge(1.4))
    ]
    candidate_names = "、".join(candidates["变量"].tolist())
    text = [
        "# 训练集优先候选变量：相关性与VIF",
        "",
        "## 范围",
        "",
        f"- 仅分析训练集内单因素P<0.05的{len(candidates)}个优先候选：{candidate_names}。",
        "- 肿瘤部位2和TG2未纳入；扩展候选（0.05≤P<0.10）未纳入。",
        "- 本次是变量去冗余诊断，不拟合结局模型，验证集未读取。",
        "",
        "## 结果",
        "",
    ]
    if high_pairs.empty:
        text.append("- 连续变量间未见|Spearman rho|≥0.70的较高相关配对。")
    else:
        text.append("- 需关注的连续变量相关配对（|rho|≥0.70）：")
        text.extend(
            f"  - {row['变量1']} vs {row['变量2']}: rho={row['Spearman_rho']:.3f}"
            for _, row in high_pairs.iterrows()
        )
    if high_vif.empty:
        text.append("- 未见需关注的VIF/GVIF信号；候选变量可暂保留进入下一轮受限多变量建模。")
    else:
        text.append("- 以下变量需要在下一步以临床优先级择一或合并处理：")
        text.extend(
            f"  - {row['变量']}: {row['判断']}"
            for _, row in high_vif.iterrows()
        )
    text.extend([
        "",
        "## 方法边界",
        "",
        "- 连续变量相关性使用Spearman秩相关，按每一对变量的完整案例计算；本诊断以|rho|阈值判定，不报告相关性P值。",
        "- VIF设计矩阵中的连续缺失值以训练集内中位数填补、分类缺失以训练集内众数填补；填补只用于诊断，不是最终建模缺失值策略。",
        "- 分类变量以哑变量进入设计矩阵；总体报告GVIF及其按自由度调整后的GVIF^(1/(2df))，因此不将T1或分化程度2的编码直接视为连续数值。",
        "- VIF/GVIF反映预测变量之间的冗余，不能证明变量与吻合口漏的独立关联。",
        "",
    ])
    path.write_text("\n".join(text), encoding="utf-8")


def main():
    parser = argparse.ArgumentParser(description="Training-only candidate correlation and VIF")
    parser.add_argument("--training-input", default=str(DEFAULT_TRAINING_INPUT))
    parser.add_argument("--candidate-input", default=str(DEFAULT_CANDIDATE_INPUT))
    parser.add_argument("--output-dir", default=str(DEFAULT_OUTPUT_DIR))
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()

    training_path = Path(args.training_input).expanduser().resolve()
    candidate_path = Path(args.candidate_input).expanduser().resolve()
    output_dir = Path(args.output_dir).expanduser().resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    workbook_path = output_dir / "训练集_优先候选变量_相关性与VIF_排除肿瘤部位.xlsx"
    report_path = output_dir / "训练集_优先候选变量_相关性与VIF_结果解读.md"
    config_path = output_dir / "训练集_优先候选变量_相关性与VIF_运行配置.json"
    manifest_path = output_dir / "SHA256SUMS.txt"
    if workbook_path.exists() and not args.overwrite:
        raise FileExistsError(f"{workbook_path} exists; pass --overwrite")

    candidates = choose_priority_candidates(candidate_path)
    continuous = candidates.loc[candidates["变量类型"].eq("连续"), "变量"].tolist()
    categorical = candidates.loc[candidates["变量类型"].eq("分类"), "变量"].tolist()
    clinical = pd.read_excel(training_path, sheet_name="临床数据", engine="openpyxl")
    missing_columns = set(candidates["变量"]) - set(clinical.columns)
    if missing_columns:
        raise ValueError(f"训练集缺少候选变量: {sorted(missing_columns)}")
    cleaned, missing_qc = clean_candidate_data(clinical, continuous, categorical)
    correlations = pairwise_spearman(cleaned, continuous)
    cat_associations = categorical_association(cleaned, categorical)
    design, term_map, imputation = make_vif_design(cleaned, continuous, categorical)
    column_vif = design_vif(design)
    gvif = termwise_gvif(design, term_map)
    notes = pd.DataFrame(
        [
            ("分析目的", "仅对训练集单因素P<0.05的优先候选变量进行去冗余诊断"),
            ("候选变量", "、".join(candidates["变量"].tolist())),
            ("未纳入", "扩展候选；肿瘤部位2；TG2"),
            ("分析集", f"训练集{len(clinical)}例；验证集未读取"),
            ("连续变量相关", "Spearman秩相关；每一对变量按配对完整案例计算。诊断以|rho|阈值判定，不报告相关性P值"),
            ("分类变量处理", "以哑变量进入VIF；分类变量本身报告GVIF及调整GVIF"),
            ("缺失处理", "仅VIF设计矩阵：连续变量训练集内中位数、分类变量训练集内众数填补；不作为最终模型策略"),
            ("VIF阈值", "连续变量：VIF<2低，2-<5轻中度，5-<10需关注，≥10高；分类变量看调整GVIF"),
            ("解释边界", "共线性诊断不代表与结局的独立关联，也不进行变量最终定稿"),
        ],
        columns=["项目", "说明"],
    )
    with pd.ExcelWriter(workbook_path, engine="openpyxl") as writer:
        notes.to_excel(writer, sheet_name="分析说明", index=False)
        candidates.to_excel(writer, sheet_name="纳入的优先候选", index=False)
        correlations.to_excel(writer, sheet_name="连续变量_Spearman相关", index=False)
        cat_associations.to_excel(writer, sheet_name="分类变量_CramersV", index=False)
        gvif.to_excel(writer, sheet_name="变量层面_GVIF", index=False)
        column_vif.to_excel(writer, sheet_name="哑变量层面_VIF", index=False)
        imputation.to_excel(writer, sheet_name="VIF缺失处理", index=False)
        missing_qc.to_excel(writer, sheet_name="候选变量缺失QC", index=False)
        design.to_excel(writer, sheet_name="VIF设计矩阵", index=False)
    style_workbook(workbook_path)
    write_report(report_path, candidates, correlations, gvif)

    config = {
        "created_at": datetime.now().isoformat(timespec="seconds"),
        "analysis_version": "training-priority-candidate-collinearity-v1",
        "training_input": str(training_path),
        "training_input_sha256": sha256(training_path),
        "candidate_input": str(candidate_path),
        "candidate_input_sha256": sha256(candidate_path),
        "analysis_set": "training set only; validation not read",
        "candidate_rule": "only univariable raw P<0.05 priority candidates",
        "excluded_variables": sorted(EXCLUDED_VARIABLES),
        "continuous_zero_rule": "all exact zeros treated as missing",
        "vif_missing_strategy": "median for continuous and mode for categorical, fitted within training data, diagnostic use only",
        "correlation_method": "pairwise complete-case Spearman for continuous variables",
        "categorical_vif_strategy": "reference-cell dummy coding; report termwise GVIF and GVIF^(1/(2df))",
        "output": str(workbook_path),
        "output_sha256": sha256(workbook_path),
        "software": {
            "python": sys.version.split()[0],
            "numpy": np.__version__,
            "pandas": pd.__version__,
            "dependencies": "numpy, pandas, openpyxl (no scipy/statsmodels required)",
        },
    }
    config_path.write_text(json.dumps(config, ensure_ascii=False, indent=2), encoding="utf-8")
    manifest_path.write_text(
        "\n".join(
            f"{sha256(path)}  {path.name}"
            for path in [workbook_path, report_path, config_path]
        ) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(config, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
