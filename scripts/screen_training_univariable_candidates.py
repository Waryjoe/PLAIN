#!/usr/bin/env python3
"""Screen clinical candidate variables using training-set univariable analysis.

This is a candidate-screening step, not a final multivariable prediction model.
"""
from plain_paths import plain_path

import argparse
import hashlib
import json
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm

import screen_training_features as common


DEFAULT_INPUT = Path(
    plain_path("PROJECT", "汇总/训练验证集/训练集_70%_574例.xlsx")
)
DEFAULT_OUTPUT_DIR = Path(
    plain_path("PROJECT", "汇总/训练集单因素候选变量筛选_排除肿瘤部位")
)
EXCLUDED_VARIABLES = {"肿瘤部位2", "TG2"}

ALL_CONTINUOUS = (
    common.CONTINUOUS_PREOPERATIVE
    + common.CONTINUOUS_POSTOPERATIVE_OR_UNCLEAR
)
ALL_CATEGORICAL = (
    common.CATEGORICAL_PREOPERATIVE
    + common.CATEGORICAL_POSTOPERATIVE_OR_UNCLEAR
)
TUMOR_RELATED = {
    "分化程度2",
    "肿瘤最长直径",
    "肿瘤次长径",
    "肿瘤短径",
    "T1",
    "N1",
    "CEA",
    "CA199",
    "CA125",
    "AFP",
    "CA242",
    "CA724",
}
CATEGORY_LABELS = {
    "性别2": {"1": "男", "2": "女"},
    "家族史2": {"1": "有", "2": "无"},
    "分化程度2": {"1": "低危", "2": "中危", "3": "高危"},
    "T1": {"1": "T1", "2": "T2", "3": "T3", "4": "T4"},
    "N1": {"0": "N0", "1": "N1", "2": "N2", "3": "N3"},
}


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def code_text(value):
    value = str(value)
    return value[:-2] if value.endswith(".0") else value


def code_label(variable, value):
    code = code_text(value)
    label = CATEGORY_LABELS.get(variable, {}).get(code)
    return f"{code} ({label})" if label else code


def candidate_tier(p_value):
    if not np.isfinite(p_value):
        return "无法评估"
    if p_value < 0.05:
        return "优先候选：P<0.05"
    if p_value < 0.10:
        return "扩展候选：0.05≤P<0.10"
    return "暂不进入候选池：P≥0.10"


def variable_role(variable):
    if variable in common.CONTINUOUS_PREOPERATIVE:
        return "术前候选"
    if variable in common.CATEGORICAL_PREOPERATIVE:
        return "术前候选"
    return "术后或时间点不明确：仅关联探索"


def run_continuous_univariable(cleaned):
    outcome = pd.to_numeric(cleaned["吻合口"], errors="raise").astype(int)
    rows = []
    for variable in ALL_CONTINUOUS:
        if variable in EXCLUDED_VARIABLES:
            continue
        values = pd.to_numeric(cleaned[variable], errors="coerce")
        keep = values.notna() & outcome.notna()
        x = values.loc[keep].astype(float)
        y = outcome.loc[keep].astype(int)
        row = {
            "变量": variable,
            "变量类型": "连续",
            "角色": variable_role(variable),
            "肿瘤相关": variable in TUMOR_RELATED,
            "完整案例数": int(keep.sum()),
            "事件数": int(y.sum()),
            "缺失数": int((~keep).sum()),
            "非漏组中位数(Q1,Q3)": common.fmt_median_iqr(
                values.loc[outcome.eq(0)]
            ),
            "漏组中位数(Q1,Q3)": common.fmt_median_iqr(
                values.loc[outcome.eq(1)]
            ),
        }
        if len(x) < 20 or y.nunique() < 2 or x.std(ddof=0) == 0:
            row.update(
                {
                    "每1SD_OR": np.nan,
                    "OR_95CI下限": np.nan,
                    "OR_95CI上限": np.nan,
                    "P值": np.nan,
                }
            )
        else:
            z = (x - x.mean()) / x.std(ddof=0)
            fit = sm.GLM(
                y.to_numpy(),
                sm.add_constant(z.to_numpy()),
                family=sm.families.Binomial(),
            ).fit()
            beta = float(fit.params[1])
            se = float(fit.bse[1])
            row.update(
                {
                    "每1SD_OR": float(np.exp(beta)),
                    "OR_95CI下限": float(np.exp(beta - 1.96 * se)),
                    "OR_95CI上限": float(np.exp(beta + 1.96 * se)),
                    "P值": float(fit.pvalues[1]),
                }
            )
        rows.append(row)
    return pd.DataFrame(rows)


def categorical_level_or(cleaned, variables):
    outcome = pd.to_numeric(cleaned["吻合口"], errors="raise").astype(int)
    rows = []
    for variable in variables:
        values = cleaned[variable]
        keep = values.notna() & outcome.notna()
        observed = values.loc[keep].astype(str)
        levels = sorted(observed.unique(), key=lambda value: float(value))
        if len(levels) < 2:
            continue
        reference = levels[0]
        dummies = pd.get_dummies(
            pd.Categorical(observed, categories=levels),
            prefix=variable,
            prefix_sep="=",
            drop_first=True,
            dtype=float,
        )
        design = sm.add_constant(dummies, has_constant="add")
        try:
            fit = sm.GLM(
                outcome.loc[keep].to_numpy(),
                design.to_numpy(dtype=float),
                family=sm.families.Binomial(),
            ).fit()
        except Exception as error:
            rows.append(
                {
                    "变量": variable,
                    "比较水平": "模型未收敛",
                    "参考水平": code_label(variable, reference),
                    "错误": str(error),
                }
            )
            continue
        for index, term in enumerate(design.columns[1:], start=1):
            beta = float(fit.params[index])
            se = float(fit.bse[index])
            rows.append(
                {
                    "变量": variable,
                    "比较水平": code_label(variable, term.split("=", 1)[1]),
                    "参考水平": code_label(variable, reference),
                    "完整案例数": int(keep.sum()),
                    "事件数": int(outcome.loc[keep].sum()),
                    "OR": float(np.exp(beta)),
                    "OR_95CI下限": float(np.exp(beta - 1.96 * se)),
                    "OR_95CI上限": float(np.exp(beta + 1.96 * se)),
                    "P值": float(fit.pvalues[index]),
                }
            )
    return pd.DataFrame(rows)


def add_fdr_and_candidates(continuous, categorical):
    tests = pd.concat(
        [
            continuous[["变量", "变量类型", "角色", "肿瘤相关", "P值"]],
            categorical[["变量", "变量类型", "角色", "肿瘤相关", "P值"]],
        ],
        ignore_index=True,
    )
    tests["BH_FDR"] = common.bh_adjust(tests["P值"])
    tests["候选层级"] = tests["P值"].map(candidate_tier)
    tests["进入下一轮候选"] = tests["P值"].lt(0.10)
    tests["肿瘤相关处理建议"] = np.where(
        tests["肿瘤相关"],
        np.where(
            tests["角色"].eq("术前候选"),
            "若确属术前可得信息，可进入后续预测候选池",
            "若为术后病理，仅用于关联探索，不用于术前预测",
        ),
        "不适用",
    )
    return tests.sort_values(["P值", "变量"], na_position="last")


def build_notes(counts, candidates):
    priority = candidates["候选层级"].eq("优先候选：P<0.05").sum()
    extended = candidates["候选层级"].eq("扩展候选：0.05≤P<0.10").sum()
    tumor = candidates.loc[
        candidates["肿瘤相关"] & candidates["进入下一轮候选"], "变量"
    ].tolist()
    return pd.DataFrame(
        [
            ("分析目的", "训练集内单因素筛选广义临床候选变量；未建立最终预测模型"),
            ("分析集", f"仅训练集：{counts['total']}例；非漏{counts['negative']}例；漏{counts['positive']}例"),
            ("排除变量", "肿瘤部位2按研究目的排除；TG2定义不清且高度提示结局，排除"),
            ("连续变量0值", "所有连续测量值的精确0按用户既定规则处理为缺失"),
            ("分类变量0值", "保留有临床含义的编码0，例如N1=0（N0）与吻合口=0"),
            ("连续变量检验", "单因素Logistic回归；OR按每增加1个标准差报告"),
            ("分类变量检验", "分类变量采用整体卡方/Fisher/Monte-Carlo检验；另报告水平OR"),
            ("候选规则", "P<0.05为优先候选；0.05≤P<0.10为扩展候选；BH-FDR仅作敏感性信息"),
            ("优先候选数量", int(priority)),
            ("扩展候选数量", int(extended)),
            ("肿瘤相关候选", "、".join(tumor) if tumor else "无"),
            ("解释边界", "单因素关联不等于独立危险因素或最终预测因子；后续需做共线性和受限的多变量建模"),
            ("验证集", "未读取，也未参与筛选"),
        ],
        columns=["项目", "说明"],
    )


def write_report(path, notes, candidates):
    priority = candidates.loc[
        candidates["候选层级"].eq("优先候选：P<0.05")
    ]
    extended = candidates.loc[
        candidates["候选层级"].eq("扩展候选：0.05≤P<0.10")
    ]

    def line_list(frame):
        return "、".join(frame["变量"].tolist()) if len(frame) else "无"

    text = [
        "# 训练集单因素候选变量筛选：排除肿瘤部位",
        "",
        "## 筛选结果",
        "",
        f"- 优先候选（P<0.05）：{line_list(priority)}。",
        f"- 扩展候选（0.05≤P<0.10）：{line_list(extended)}。",
        "- 肿瘤相关变量需要按可获得时间点决定用途：术后病理变量仅用于关联探索，不能用于术前预测。",
        "",
        "## 使用边界",
        "",
        "- 本结果仅用于构建下一轮候选池，不把单因素结果称为独立危险因素。",
        "- 验证集未被读取；后续变量去冗余、惩罚回归与模型锁定均必须只在训练集内完成。",
        "- BH-FDR列用于评估多重检验后的稳健性，不替代本次宽松候选筛选规则。",
        "",
    ]
    path.write_text("\n".join(text), encoding="utf-8")


def main():
    parser = argparse.ArgumentParser(
        description="Training-only univariable clinical candidate screening"
    )
    parser.add_argument("--input", default=str(DEFAULT_INPUT))
    parser.add_argument("--output-dir", default=str(DEFAULT_OUTPUT_DIR))
    parser.add_argument("--seed", type=int, default=20260729)
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()

    input_path = Path(args.input).expanduser().resolve()
    output_dir = Path(args.output_dir).expanduser().resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    workbook_path = output_dir / "训练集_单因素候选变量筛选_排除肿瘤部位.xlsx"
    report_path = output_dir / "训练集_单因素候选变量筛选_结果解读.md"
    config_path = output_dir / "训练集_单因素候选变量筛选_运行配置.json"
    manifest_path = output_dir / "SHA256SUMS.txt"
    if workbook_path.exists() and not args.overwrite:
        raise FileExistsError(f"{workbook_path} exists; pass --overwrite")

    clinical = pd.read_excel(input_path, sheet_name="临床数据", engine="openpyxl")
    outcome = pd.to_numeric(clinical["吻合口"], errors="raise").astype(int)
    if set(outcome.unique()) != {0, 1}:
        raise ValueError("吻合口 must contain both 0 and 1")

    common.ZERO_AS_MISSING = set(ALL_CONTINUOUS)
    common.CATEGORICAL_PREOPERATIVE = [
        variable
        for variable in common.CATEGORICAL_PREOPERATIVE
        if variable not in EXCLUDED_VARIABLES
    ]
    common.CATEGORICAL_POSTOPERATIVE_OR_UNCLEAR = [
        variable
        for variable in common.CATEGORICAL_POSTOPERATIVE_OR_UNCLEAR
        if variable not in EXCLUDED_VARIABLES
    ]
    cleaned, qc = common.clean_clinical(clinical)
    continuous = run_continuous_univariable(cleaned)
    categorical_raw, categorical_levels = common.categorical_univariable(
        cleaned, args.seed
    )
    categorical = categorical_raw.rename(columns={"角色": "角色"}).copy()
    categorical["变量类型"] = "分类"
    categorical["肿瘤相关"] = categorical["变量"].isin(TUMOR_RELATED)
    categorical_levels_or = categorical_level_or(
        cleaned, categorical["变量"].tolist()
    )
    candidates = add_fdr_and_candidates(continuous, categorical)
    notes = build_notes(
        {
            "total": int(len(outcome)),
            "negative": int(outcome.eq(0).sum()),
            "positive": int(outcome.eq(1).sum()),
        },
        candidates,
    )

    with pd.ExcelWriter(workbook_path, engine="openpyxl") as writer:
        notes.to_excel(writer, sheet_name="分析说明", index=False)
        candidates.to_excel(writer, sheet_name="候选变量清单", index=False)
        candidates.loc[candidates["肿瘤相关"]].to_excel(
            writer, sheet_name="肿瘤相关变量", index=False
        )
        continuous.sort_values("P值", na_position="last").to_excel(
            writer, sheet_name="单因素_连续变量", index=False
        )
        categorical.sort_values("P值", na_position="last").to_excel(
            writer, sheet_name="单因素_分类变量整体", index=False
        )
        categorical_levels.to_excel(
            writer, sheet_name="分类变量构成", index=False
        )
        categorical_levels_or.to_excel(
            writer, sheet_name="分类变量水平OR", index=False
        )
        qc.to_excel(writer, sheet_name="缺失值与0值QC", index=False)
    common.style_workbook(workbook_path)
    write_report(report_path, notes, candidates)

    config = {
        "created_at": datetime.now().isoformat(timespec="seconds"),
        "analysis_version": "training-univariable-candidate-v1",
        "input": str(input_path),
        "input_sha256": sha256(input_path),
        "output": str(workbook_path),
        "output_sha256": sha256(workbook_path),
        "report": str(report_path),
        "seed": args.seed,
        "analysis_set": "training set only; validation not read",
        "excluded_variables": sorted(EXCLUDED_VARIABLES),
        "continuous_zero_rule": "all exact zeros treated as missing",
        "categorical_zero_rule": "retain coded zeros with clinical meaning",
        "candidate_rule": "priority raw P<0.05; extended raw P<0.10; FDR sensitivity only",
        "software": {
            "python": __import__("sys").version.split()[0],
            "numpy": np.__version__,
            "pandas": pd.__version__,
            "statsmodels": sm.__version__,
        },
    }
    config_path.write_text(json.dumps(config, ensure_ascii=False, indent=2), encoding="utf-8")
    manifest_path.write_text(
        "\n".join(
            f"{sha256(path)}  {path.name}"
            for path in [workbook_path, report_path, config_path]
        )
        + "\n",
        encoding="utf-8",
    )
    print(json.dumps(config, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
