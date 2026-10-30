#!/usr/bin/env python3
from plain_paths import plain_path
import argparse
import hashlib
import json
import math
from collections import Counter, defaultdict
from datetime import datetime
from pathlib import Path

import joblib
import numpy as np
import openpyxl
import pandas as pd
import SimpleITK as sitk
from openpyxl.styles import Alignment, Font, PatternFill
from openpyxl.utils import get_column_letter
from scipy.stats import (
    MonteCarloMethod,
    chi2_contingency,
    fisher_exact,
    mannwhitneyu,
)
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import make_scorer, average_precision_score
from sklearn.model_selection import GridSearchCV, RepeatedStratifiedKFold
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
import statsmodels.api as sm


DEFAULT_INPUT = Path(
    plain_path("PROJECT", "汇总/训练验证集/训练集_70%_574例.xlsx")
)
DEFAULT_OUTPUT_DIR = Path(
    plain_path("PROJECT", "汇总/训练集变量筛选")
)
DEFAULT_SEED = 20260729

CONTINUOUS_PREOPERATIVE = [
    "年龄",
    "BMI",
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

CONTINUOUS_POSTOPERATIVE_OR_UNCLEAR = [
    "肿瘤最长直径",
    "肿瘤次长径",
    "肿瘤短径",
]

CATEGORICAL_PREOPERATIVE = [
    "性别2",
    "家族史2",
    "肿瘤部位2",
]

CATEGORICAL_POSTOPERATIVE_OR_UNCLEAR = [
    "分化程度2",
    "T1",
    "N1",
]

# Exact zero cannot represent a valid measurement for these variables. It is
# treated as a missing-value placeholder. Exact zeros are retained for tumor
# markers and basophil/eosinophil percentages because zero may be reportable.
ZERO_AS_MISSING = {
    "BMI",
    "肿瘤最长直径",
    "肿瘤次长径",
    "肿瘤短径",
    "红细胞比容Hct",
    "血红蛋白Hb",
    "红细胞RBC",
    "淋巴细胞Lymph",
    "血小板计数PLT",
    "白细胞WBC",
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
    "前白蛋白PA",
    "ALT",
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
}

NEGATIVE_AS_MISSING = {"C反应蛋白CRP"}


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def bh_adjust(p_values):
    values = np.asarray(p_values, dtype=float)
    adjusted = np.full(values.shape, np.nan)
    valid = np.isfinite(values)
    if not valid.any():
        return adjusted
    p = values[valid]
    order = np.argsort(p)
    ranked = p[order]
    q = ranked * len(ranked) / np.arange(1, len(ranked) + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    q = np.clip(q, 0, 1)
    restored = np.empty_like(q)
    restored[order] = q
    adjusted[valid] = restored
    return adjusted


def fmt_median_iqr(values):
    values = pd.to_numeric(values, errors="coerce").dropna().to_numpy()
    if len(values) == 0:
        return ""
    q1, median, q3 = np.quantile(values, [0.25, 0.5, 0.75])
    return f"{median:.3g} ({q1:.3g}, {q3:.3g})"


def clean_clinical(clinical):
    cleaned = clinical.copy()
    continuous = (
        CONTINUOUS_PREOPERATIVE + CONTINUOUS_POSTOPERATIVE_OR_UNCLEAR
    )
    qc_rows = []
    y = pd.to_numeric(cleaned["吻合口"], errors="raise").astype(int)
    for variable in continuous:
        raw = pd.to_numeric(cleaned[variable], errors="coerce")
        value = raw.copy()
        zero_replaced = pd.Series(False, index=value.index)
        negative_replaced = pd.Series(False, index=value.index)
        if variable in ZERO_AS_MISSING:
            zero_replaced = value.eq(0)
            value.loc[zero_replaced] = np.nan
        if variable in NEGATIVE_AS_MISSING:
            negative_replaced = value.lt(0)
            value.loc[negative_replaced] = np.nan
        cleaned[variable] = value
        qc_rows.append(
            {
                "变量": variable,
                "角色": (
                    "术前候选"
                    if variable in CONTINUOUS_PREOPERATIVE
                    else "术后或时间点不明确，仅探索"
                ),
                "原始缺失数": int(raw.isna().sum()),
                "原始精确0数": int(raw.eq(0).sum()),
                "按规则转缺失数": int(
                    (zero_replaced | negative_replaced).sum()
                ),
                "清洗后缺失数": int(value.isna().sum()),
                "清洗后缺失率": float(value.isna().mean()),
                "非漏组清洗后缺失率": float(value[y.eq(0)].isna().mean()),
                "漏组清洗后缺失率": float(value[y.eq(1)].isna().mean()),
                "清洗后最小值": (
                    float(value.min()) if value.notna().any() else np.nan
                ),
                "清洗后最大值": (
                    float(value.max()) if value.notna().any() else np.nan
                ),
            }
        )
    return cleaned, pd.DataFrame(qc_rows)


def univariable_logistic_or(values, y):
    values = pd.to_numeric(values, errors="coerce")
    keep = values.notna() & y.notna()
    x = values[keep].astype(float)
    outcome = y[keep].astype(int)
    if len(x) < 20 or outcome.nunique() < 2 or x.std(ddof=0) == 0:
        return np.nan, np.nan, np.nan
    z = (x - x.mean()) / x.std(ddof=0)
    design = sm.add_constant(z.to_numpy())
    try:
        fit = sm.GLM(
            outcome.to_numpy(),
            design,
            family=sm.families.Binomial(),
        ).fit()
        beta = fit.params[1]
        se = fit.bse[1]
        return (
            float(np.exp(beta)),
            float(np.exp(beta - 1.96 * se)),
            float(np.exp(beta + 1.96 * se)),
        )
    except Exception:
        return np.nan, np.nan, np.nan


def continuous_univariable(cleaned):
    y = pd.to_numeric(cleaned["吻合口"], errors="raise").astype(int)
    rows = []
    for variable in (
        CONTINUOUS_PREOPERATIVE + CONTINUOUS_POSTOPERATIVE_OR_UNCLEAR
    ):
        values = pd.to_numeric(cleaned[variable], errors="coerce")
        negative = values[y.eq(0)].dropna()
        positive = values[y.eq(1)].dropna()
        if len(negative) and len(positive):
            test = mannwhitneyu(
                positive,
                negative,
                alternative="two-sided",
                method="asymptotic",
            )
            u_statistic = float(test.statistic)
            p_value = float(test.pvalue)
            rank_biserial = (
                2 * u_statistic / (len(positive) * len(negative)) - 1
            )
        else:
            u_statistic = p_value = rank_biserial = np.nan
        odds_ratio, ci_low, ci_high = univariable_logistic_or(values, y)
        rows.append(
            {
                "变量": variable,
                "角色": (
                    "术前候选"
                    if variable in CONTINUOUS_PREOPERATIVE
                    else "术后或时间点不明确，仅探索"
                ),
                "总例数": int(values.notna().sum()),
                "非漏组例数": int(negative.size),
                "漏组例数": int(positive.size),
                "非漏组_中位数(Q1,Q3)": fmt_median_iqr(negative),
                "漏组_中位数(Q1,Q3)": fmt_median_iqr(positive),
                "Mann-Whitney_U": u_statistic,
                "秩二列相关": rank_biserial,
                "每1SD_OR": odds_ratio,
                "OR_95CI下限": ci_low,
                "OR_95CI上限": ci_high,
                "P值": p_value,
            }
        )
    result = pd.DataFrame(rows)
    result["BH_FDR"] = bh_adjust(result["P值"])
    result["原始P<0.05"] = result["P值"].lt(0.05)
    result["FDR<0.05"] = result["BH_FDR"].lt(0.05)
    return result.sort_values(["BH_FDR", "P值", "变量"])


def categorical_univariable(cleaned, seed=DEFAULT_SEED):
    y = pd.to_numeric(cleaned["吻合口"], errors="raise").astype(int)
    global_rows = []
    level_rows = []
    variables = (
        CATEGORICAL_PREOPERATIVE + CATEGORICAL_POSTOPERATIVE_OR_UNCLEAR
    )
    for variable_index, variable in enumerate(variables):
        values = cleaned[variable].copy()
        values = values.where(values.notna(), "缺失").astype(str)
        table = pd.crosstab(y, values)
        table = table.reindex(index=[0, 1], fill_value=0)
        expected = None
        if table.shape[1] >= 2:
            chi2, chi_p, _, expected = chi2_contingency(table)
            if table.shape == (2, 2) and (expected < 5).any():
                statistic, p_value = fisher_exact(table.to_numpy())
                test_name = "Fisher精确检验"
            elif (expected < 5).any():
                monte_carlo = MonteCarloMethod(
                    n_resamples=99999,
                    rng=np.random.default_rng(seed + variable_index),
                )
                statistic, p_value, _, _ = chi2_contingency(
                    table,
                    correction=False,
                    method=monte_carlo,
                )
                test_name = "Pearson卡方Monte-Carlo检验(99999次)"
            else:
                statistic, p_value = chi2, chi_p
                test_name = "Pearson卡方检验"
            total = table.to_numpy().sum()
            min_dim = min(table.shape) - 1
            cramers_v = (
                math.sqrt(float(chi2) / (total * min_dim))
                if total > 0 and min_dim > 0
                else np.nan
            )
        else:
            statistic = p_value = cramers_v = np.nan
            test_name = "不可估计"
        global_rows.append(
            {
                "变量": variable,
                "角色": (
                    "术前候选"
                    if variable in CATEGORICAL_PREOPERATIVE
                    else "术后或时间点不明确，仅探索"
                ),
                "检验": test_name,
                "统计量": statistic,
                "Cramers_V": cramers_v,
                "P值": p_value,
                "水平数": int(table.shape[1]),
                "任一期望频数<5": (
                    bool((expected < 5).any()) if expected is not None else True
                ),
            }
        )
        for level in table.columns:
            n0 = int(table.loc[0, level])
            n1 = int(table.loc[1, level])
            level_rows.append(
                {
                    "变量": variable,
                    "水平": level,
                    "非漏组_n": n0,
                    "非漏组_%": 100 * n0 / int(y.eq(0).sum()),
                    "漏组_n": n1,
                    "漏组_%": 100 * n1 / int(y.eq(1).sum()),
                }
            )
    global_result = pd.DataFrame(global_rows)
    global_result["BH_FDR"] = bh_adjust(global_result["P值"])
    global_result["原始P<0.05"] = global_result["P值"].lt(0.05)
    global_result["FDR<0.05"] = global_result["BH_FDR"].lt(0.05)
    return (
        global_result.sort_values(["BH_FDR", "P值", "变量"]),
        pd.DataFrame(level_rows),
    )


def encode_clinical(cleaned, preoperative_only):
    continuous = list(CONTINUOUS_PREOPERATIVE)
    categorical = list(CATEGORICAL_PREOPERATIVE)
    if not preoperative_only:
        continuous += CONTINUOUS_POSTOPERATIVE_OR_UNCLEAR
        categorical += CATEGORICAL_POSTOPERATIVE_OR_UNCLEAR

    pieces = []
    base_by_column = {}
    for variable in continuous:
        column = pd.to_numeric(cleaned[variable], errors="coerce").astype(float)
        pieces.append(column.rename(variable))
        base_by_column[variable] = variable
    for variable in categorical:
        values = cleaned[variable].where(
            cleaned[variable].notna(), "缺失"
        ).astype(str)
        dummies = pd.get_dummies(
            values,
            prefix=variable,
            prefix_sep="=",
            drop_first=True,
            dtype=float,
        )
        pieces.append(dummies)
        for column in dummies.columns:
            base_by_column[column] = variable
    matrix = pd.concat(pieces, axis=1)
    return matrix, base_by_column


def stratified_subsample(y, fraction, rng):
    indices = []
    for label in sorted(np.unique(y)):
        label_indices = np.flatnonzero(y == label)
        size = max(2, int(round(len(label_indices) * fraction)))
        indices.extend(rng.choice(label_indices, size=size, replace=False))
    indices = np.asarray(indices)
    rng.shuffle(indices)
    return indices


def tune_lasso(X, y, seed, n_jobs):
    pipeline = Pipeline(
        [
            ("imputer", SimpleImputer(strategy="median")),
            ("scaler", StandardScaler()),
            (
                "model",
                LogisticRegression(
                    penalty="l1",
                    solver="liblinear",
                    class_weight="balanced",
                    max_iter=5000,
                    tol=1e-5,
                    random_state=seed,
                ),
            ),
        ]
    )
    cv = RepeatedStratifiedKFold(
        n_splits=5,
        n_repeats=3,
        random_state=seed,
    )
    search = GridSearchCV(
        pipeline,
        param_grid={
            "model__C": np.logspace(-3, 1, 17),
        },
        scoring=make_scorer(average_precision_score, response_method="predict_proba"),
        cv=cv,
        n_jobs=n_jobs,
        refit=True,
        error_score="raise",
    )
    search.fit(X, y)
    return search


def stability_selection(
    X,
    y,
    base_by_column,
    best_params,
    seed,
    repetitions,
    n_jobs,
):
    rng = np.random.default_rng(seed)
    transformed_counter = Counter()
    base_counter = Counter()
    coefficient_sums = defaultdict(float)
    model = Pipeline(
        [
            ("imputer", SimpleImputer(strategy="median")),
            ("scaler", StandardScaler()),
            (
                "model",
                LogisticRegression(
                    penalty="l1",
                    solver="liblinear",
                    class_weight="balanced",
                    max_iter=5000,
                    tol=1e-5,
                    random_state=seed,
                    C=float(best_params["model__C"]),
                ),
            ),
        ]
    )
    for repetition in range(repetitions):
        indices = stratified_subsample(y, 0.8, rng)
        fitted = model.fit(X.iloc[indices], y[indices])
        coefficients = fitted.named_steps["model"].coef_[0]
        selected = np.flatnonzero(np.abs(coefficients) > 1e-10)
        selected_bases = set()
        for index in selected:
            column = X.columns[index]
            transformed_counter[column] += 1
            coefficient_sums[column] += float(coefficients[index])
            selected_bases.add(base_by_column[column])
        for base in selected_bases:
            base_counter[base] += 1

    transformed_rows = []
    for column in X.columns:
        count = transformed_counter[column]
        transformed_rows.append(
            {
                "派生特征": column,
                "原始变量": base_by_column[column],
                "入选次数": count,
                "稳定选择频率": count / repetitions,
                "入选时平均标准化系数": (
                    coefficient_sums[column] / count if count else 0.0
                ),
            }
        )
    base_rows = [
        {
            "原始变量": base,
            "入选次数": base_counter[base],
            "稳定选择频率": base_counter[base] / repetitions,
            "达到0.60阈值": base_counter[base] / repetitions >= 0.60,
        }
        for base in sorted(set(base_by_column.values()))
    ]
    return (
        pd.DataFrame(transformed_rows).sort_values(
            ["稳定选择频率", "派生特征"], ascending=[False, True]
        ),
        pd.DataFrame(base_rows).sort_values(
            ["稳定选择频率", "原始变量"], ascending=[False, True]
        ),
    )


def load_radiomics(mapping, cache_path):
    if cache_path.is_file():
        cache = np.load(cache_path, allow_pickle=True)
        return (
            cache["X"].astype(float),
            cache["case_ids"].astype(str),
            cache["feature_names"].astype(str),
            pd.DataFrame(json.loads(str(cache["qc_json"].item()))),
        )

    matrices = []
    case_ids = []
    feature_names = None
    qc_rows = []
    metadata = {
        "case_id",
        "roi",
        "image_path",
        "source_mask_path",
        "aligned_mask_path",
        "source_label_1_voxels",
        "aligned_label_1_voxels",
    }
    for row_number, row in mapping.reset_index(drop=True).iterrows():
        path = Path(row["wide文件"])
        if not path.is_file():
            raise FileNotFoundError(path)
        wide = pd.read_excel(path, engine="openpyxl")
        if set(wide["roi"]) != {"stomach", "tumor"}:
            raise ValueError(f"Unexpected ROI rows in {path}")
        pure_features = [column for column in wide.columns if column not in metadata]
        wide = wide.set_index("roi").loc[["stomach", "tumor"]]
        names = np.asarray(
            [
                f"{roi}__{feature}"
                for roi in ["stomach", "tumor"]
                for feature in pure_features
            ],
            dtype=object,
        )
        values = np.concatenate(
            [
                pd.to_numeric(wide.loc[roi, pure_features], errors="coerce")
                .to_numpy(dtype=float)
                for roi in ["stomach", "tumor"]
            ]
        )
        if feature_names is None:
            feature_names = names
        elif not np.array_equal(feature_names, names):
            raise ValueError(f"Radiomics feature drift in {path}")
        matrices.append(values)
        case_ids.append(str(row["影像组学case_id"]))
        image_path = Path(str(wide.iloc[0]["image_path"]))
        image = sitk.ReadImage(str(image_path))
        size = image.GetSize()
        spacing = image.GetSpacing()
        qc_rows.append(
            {
                "case_id": str(row["影像组学case_id"]),
                "wide文件": str(path),
                "image_path": str(image_path),
                "image_pixel_type": image.GetPixelIDTypeAsString(),
                "size_x": int(size[0]),
                "size_y": int(size[1]),
                "size_z": int(size[2]),
                "spacing_x": float(spacing[0]),
                "spacing_y": float(spacing[1]),
                "spacing_z": float(spacing[2]),
                "ROI数": int(wide.shape[0]),
                "每ROI特征数": int(len(pure_features)),
                "总特征数": int(len(values)),
                "非有限值数": int((~np.isfinite(values)).sum()),
            }
        )
        if (row_number + 1) % 50 == 0:
            print(
                f"Loaded radiomics {row_number + 1}/{len(mapping)}",
                flush=True,
            )
    X = np.vstack(matrices)
    qc = pd.DataFrame(qc_rows)
    np.savez_compressed(
        cache_path,
        X=X,
        case_ids=np.asarray(case_ids, dtype=object),
        feature_names=feature_names,
        qc_json=json.dumps(qc.to_dict(orient="records"), ensure_ascii=False),
    )
    return X, np.asarray(case_ids), feature_names.astype(str), qc


def imaging_batch_summary(radiomics_case_qc, mapping):
    qc = radiomics_case_qc.copy()
    qc["吻合口"] = pd.to_numeric(mapping["吻合口"], errors="raise").astype(int)
    rows = []
    for label, name in [(0, "非漏组"), (1, "漏组")]:
        subset = qc.loc[qc["吻合口"].eq(label)]
        rows.append(
            {
                "分组": name,
                "例数": int(len(subset)),
                "像素类型及例数": "; ".join(
                    f"{key}: {value}"
                    for key, value in subset["image_pixel_type"]
                    .value_counts()
                    .items()
                ),
                "size_x中位数": float(subset["size_x"].median()),
                "size_y中位数": float(subset["size_y"].median()),
                "size_z中位数": float(subset["size_z"].median()),
                "spacing_x中位数": float(subset["spacing_x"].median()),
                "spacing_y中位数": float(subset["spacing_y"].median()),
                "spacing_z中位数": float(subset["spacing_z"].median()),
                "不同spacing组合数": int(
                    subset[
                        ["spacing_x", "spacing_y", "spacing_z"]
                    ].drop_duplicates().shape[0]
                ),
            }
        )
    return pd.DataFrame(rows)


def radiomics_univariable(X, y, feature_names):
    rows = []
    negative = X[y == 0]
    positive = X[y == 1]
    for index, feature in enumerate(feature_names):
        x0 = negative[:, index]
        x1 = positive[:, index]
        finite0 = x0[np.isfinite(x0)]
        finite1 = x1[np.isfinite(x1)]
        if (
            len(finite0) == 0
            or len(finite1) == 0
            or np.nanstd(X[:, index]) == 0
        ):
            u_statistic = p_value = effect = auc_directional = np.nan
        else:
            test = mannwhitneyu(
                finite1,
                finite0,
                alternative="two-sided",
                method="asymptotic",
            )
            u_statistic = float(test.statistic)
            p_value = float(test.pvalue)
            auc_directional = u_statistic / (len(finite1) * len(finite0))
            effect = 2 * auc_directional - 1
        roi, feature_name = feature.split("__", 1)
        rows.append(
            {
                "完整特征名": feature,
                "ROI": roi,
                "PyRadiomics特征名": feature_name,
                "非漏组中位数": float(np.nanmedian(finite0)),
                "漏组中位数": float(np.nanmedian(finite1)),
                "Mann-Whitney_U": u_statistic,
                "方向性AUC": auc_directional,
                "秩二列相关": effect,
                "P值": p_value,
            }
        )
    result = pd.DataFrame(rows)
    result["BH_FDR"] = bh_adjust(result["P值"])
    result["原始P<0.05"] = result["P值"].lt(0.05)
    result["FDR<0.05"] = result["BH_FDR"].lt(0.05)
    return result.sort_values(
        ["BH_FDR", "P值", "完整特征名"]
    ).reset_index(drop=True)


def correlation_filter(X, univariable, feature_names, threshold):
    variances = np.nanvar(X, axis=0)
    valid_indices = np.flatnonzero(np.isfinite(variances) & (variances > 0))
    values = X[:, valid_indices]
    medians = np.nanmedian(values, axis=0)
    missing_rows, missing_cols = np.where(~np.isfinite(values))
    values = values.copy()
    values[missing_rows, missing_cols] = medians[missing_cols]
    ranks = pd.DataFrame(values).rank(method="average").to_numpy(dtype=float)
    rank_nonconstant = np.nanstd(ranks, axis=0) > 0
    valid_indices = valid_indices[rank_nonconstant]
    ranks = ranks[:, rank_nonconstant]
    order_table = univariable.set_index("完整特征名").loc[
        feature_names[valid_indices]
    ]
    order = valid_indices[
        np.argsort(
            -np.nan_to_num(
                np.abs(order_table["秩二列相关"].to_numpy()),
                nan=-1,
            )
        )
    ]
    correlation = np.corrcoef(ranks, rowvar=False)
    position = {original: local for local, original in enumerate(valid_indices)}

    kept = []
    dropped_rows = []
    for original_index in order:
        local_index = position[original_index]
        if not kept:
            kept.append(original_index)
            continue
        kept_local = [position[index] for index in kept]
        correlations = np.abs(correlation[local_index, kept_local])
        max_position = int(np.nanargmax(correlations))
        max_correlation = float(correlations[max_position])
        if max_correlation >= threshold:
            representative = kept[max_position]
            dropped_rows.append(
                {
                    "被剔除特征": feature_names[original_index],
                    "保留代表特征": feature_names[representative],
                    "Spearman绝对相关系数": max_correlation,
                }
            )
        else:
            kept.append(original_index)
    return np.asarray(kept), pd.DataFrame(dropped_rows), valid_indices


def radiomics_stability(
    X,
    y,
    feature_names,
    seed,
    repetitions,
    n_jobs,
):
    frame = pd.DataFrame(X, columns=feature_names)
    base_by_column = {column: column for column in frame.columns}
    search = tune_lasso(frame, y, seed, n_jobs)
    transformed, _ = stability_selection(
        frame,
        y,
        base_by_column,
        search.best_params_,
        seed + 101,
        repetitions,
        n_jobs,
    )
    transformed = transformed.rename(columns={"派生特征": "完整特征名"})
    transformed["达到0.60阈值"] = transformed["稳定选择频率"].ge(0.60)
    return search, transformed


def write_dataframe(writer, dataframe, sheet_name):
    dataframe.to_excel(writer, sheet_name=sheet_name, index=False)


def style_workbook(path):
    workbook = openpyxl.load_workbook(path)
    header_fill = PatternFill("solid", fgColor="1F4E78")
    header_font = Font(color="FFFFFF", bold=True)
    for worksheet in workbook.worksheets:
        worksheet.freeze_panes = "A2"
        worksheet.auto_filter.ref = worksheet.dimensions
        for cell in worksheet[1]:
            cell.fill = header_fill
            cell.font = header_font
            cell.alignment = Alignment(
                horizontal="center", vertical="center", wrap_text=True
            )
        for column_cells in worksheet.columns:
            letter = get_column_letter(column_cells[0].column)
            max_length = max(
                len(str(cell.value)) if cell.value is not None else 0
                for cell in list(column_cells)[:300]
            )
            worksheet.column_dimensions[letter].width = min(
                max(max_length + 2, 10), 42
            )
    workbook.save(path)


def analysis_notes(args, clinical_counts, summary):
    rows = [
        (
            "当前结论",
            "方法学阻断：病例来源与吻合口结局完全重合。当前结果只能用于"
            "数据审计和探索，不能据此锁定最终影像组学特征或报告模型性能",
        ),
        ("分析对象", "仅训练集；验证集未读取、未参与任何筛选或调参"),
        (
            "样本与事件",
            f"{clinical_counts['total']}例；非漏{clinical_counts['negative']}例；"
            f"吻合口漏{clinical_counts['positive']}例",
        ),
        (
            "临床差异检验",
            "连续变量Mann-Whitney U；分类变量Pearson卡方或Fisher；"
            "同时报告BH-FDR和效应量",
        ),
        (
            "临床缺失处理",
            "生理上不可能的精确0按缺失占位符处理；模型插补在每次训练内完成",
        ),
        (
            "临床稳定筛选",
            f"LASSO Logistic；5折交叉验证重复3次调参；"
            f"{args.stability_repetitions}次分层80%子抽样；稳定阈值0.60",
        ),
        (
            "影像组学初筛",
            f"常量过滤；Spearman |rho|>={args.correlation_threshold}冗余过滤；"
            f"按训练集单因素效应排序后最多保留{args.radiomics_top_k}项进入稳定筛选",
        ),
        (
            "影像组学稳定筛选",
            f"LASSO Logistic；5折交叉验证重复3次调参；"
            f"{args.stability_repetitions}次分层80%子抽样；稳定阈值0.60",
        ),
        ("类别不平衡", "所有惩罚Logistic使用class_weight=balanced"),
        ("固定随机种子", str(args.seed)),
        (
            "TG2",
            "定义不清且与结局高度重合，作为疑似标签泄漏变量排除",
        ),
        (
            "时间点风险",
            "分化程度2、T1、N1和病理肿瘤径仅作探索；若目标是术前预测，"
            "除非证明术前可获得，否则不得进入术前模型",
        ),
        (
            "ICC限制",
            "当前没有重复勾画/重复扫描数据，无法完成ICC稳定性过滤；"
            "入模前仍需补做分割重现性评估",
        ),
        (
            "重要解释",
            "P值表描述组间差异；稳定筛选频率描述候选变量可重复性。"
            "两者均不是独立验证性能，也不证明因果关系",
        ),
        (
            "批次处理限制",
            "当数据来源与结局100%共线时，ComBat或普通协变量校正无法区分"
            "真实结局效应与来源效应；需要同一来源中同时包含漏与非漏病例",
        ),
    ]
    for key, value in summary.items():
        rows.append((key, value))
    return pd.DataFrame(rows, columns=["项目", "说明"])


def leakage_and_source_audit(clinical, mapping):
    audit_rows = []
    case_ids = mapping["影像组学case_id"].astype(str)
    sources = np.select(
        [
            case_ids.str.startswith("P"),
            case_ids.str.startswith("gc_ZJCH_"),
        ],
        ["P系列原漏组影像来源", "gc_ZJCH系列原非漏组影像来源"],
        default="其他来源",
    )
    outcome = pd.to_numeric(mapping["吻合口"], errors="raise").astype(int)
    source_table = pd.crosstab(
        pd.Series(sources, name="数据来源"),
        outcome.rename("吻合口"),
    ).reindex(columns=[0, 1], fill_value=0)
    source_accuracy = (
        float(source_table.max(axis=1).sum() / source_table.to_numpy().sum())
        if source_table.to_numpy().sum()
        else np.nan
    )
    for source, row in source_table.iterrows():
        audit_rows.append(
            {
                "审计项目": "病例来源与结局",
                "分层": source,
                "非漏组_n": int(row.get(0, 0)),
                "漏组_n": int(row.get(1, 0)),
                "说明": "仅凭来源预测结局的训练集准确率"
                f"={source_accuracy:.3f}",
            }
        )

    baseline_dates = pd.to_datetime(
        clinical["基线-首诊时间"], errors="coerce"
    )
    clinical_outcome = pd.to_numeric(
        clinical["吻合口"], errors="raise"
    ).astype(int)
    for label, name in [(0, "非漏组"), (1, "漏组")]:
        values = baseline_dates[clinical_outcome.eq(label)].dropna()
        audit_rows.append(
            {
                "审计项目": "基线日期",
                "分层": name,
                "非漏组_n": (
                    int(values.size) if label == 0 else np.nan
                ),
                "漏组_n": int(values.size) if label == 1 else np.nan,
                "说明": (
                    f"可用{values.size}例；缺失"
                    f"{int(clinical_outcome.eq(label).sum() - values.size)}例；"
                    f"范围{values.min()}至{values.max()}；"
                    f"中位数{values.median()}"
                    if values.size
                    else "无可用日期"
                ),
            }
        )

    tg2_table = pd.crosstab(
        clinical_outcome.rename("吻合口"),
        clinical["TG2"].rename("TG2"),
    ).reindex(index=[0, 1], fill_value=0)
    for label, row in tg2_table.iterrows():
        audit_rows.append(
            {
                "审计项目": "TG2疑似标签泄漏",
                "分层": "非漏组" if label == 0 else "漏组",
                "非漏组_n": int(row.sum()) if label == 0 else np.nan,
                "漏组_n": int(row.sum()) if label == 1 else np.nan,
                "说明": "; ".join(
                    f"TG2={level}: {int(count)}"
                    for level, count in row.items()
                ),
            }
        )
    return pd.DataFrame(audit_rows), source_accuracy


def main():
    parser = argparse.ArgumentParser(
        description="Leakage-safe feature screening using training data only"
    )
    parser.add_argument("--input", default=str(DEFAULT_INPUT))
    parser.add_argument("--output-dir", default=str(DEFAULT_OUTPUT_DIR))
    parser.add_argument("--seed", type=int, default=DEFAULT_SEED)
    parser.add_argument("--correlation-threshold", type=float, default=0.90)
    parser.add_argument("--radiomics-top-k", type=int, default=200)
    parser.add_argument("--stability-repetitions", type=int, default=100)
    parser.add_argument("--n-jobs", type=int, default=-1)
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()

    input_path = Path(args.input).expanduser().resolve()
    output_dir = Path(args.output_dir).expanduser().resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    workbook_path = output_dir / "训练集变量筛选_临床与影像组学.xlsx"
    cache_path = output_dir / "训练集影像组学矩阵_cache.npz"
    model_path = (
        output_dir / "训练集变量筛选_探索性模型参数_禁止直接使用.joblib"
    )
    config_path = output_dir / "训练集变量筛选_运行配置.json"
    if workbook_path.exists() and not args.overwrite:
        raise FileExistsError(
            f"{workbook_path} exists; use --overwrite after confirmation"
        )

    clinical = pd.read_excel(
        input_path, sheet_name="临床数据", engine="openpyxl"
    )
    mapping = pd.read_excel(
        input_path, sheet_name="病例映射", engine="openpyxl"
    )
    y = pd.to_numeric(clinical["吻合口"], errors="raise").astype(int).to_numpy()
    if len(clinical) != len(mapping):
        raise ValueError("Clinical and mapping row counts differ")
    if not np.array_equal(
        pd.to_numeric(mapping["吻合口"], errors="raise").astype(int).to_numpy(),
        y,
    ):
        raise ValueError("Outcome order differs between clinical and mapping sheets")
    if set(np.unique(y)) != {0, 1}:
        raise ValueError("Outcome must contain 0 and 1")
    leakage_audit, source_accuracy = leakage_and_source_audit(
        clinical, mapping
    )

    cleaned, clinical_qc = clean_clinical(clinical)
    continuous_results = continuous_univariable(cleaned)
    categorical_results, categorical_levels = categorical_univariable(
        cleaned, args.seed
    )
    clinical_joint_p = pd.concat(
        [
            continuous_results[["变量", "P值"]].assign(类型="连续"),
            categorical_results[["变量", "P值"]].assign(类型="分类"),
        ],
        ignore_index=True,
    )
    clinical_joint_p["临床全变量_BH_FDR"] = bh_adjust(
        clinical_joint_p["P值"]
    )
    continuous_fdr = clinical_joint_p.loc[
        clinical_joint_p["类型"].eq("连续"),
        ["变量", "临床全变量_BH_FDR"],
    ]
    categorical_fdr = clinical_joint_p.loc[
        clinical_joint_p["类型"].eq("分类"),
        ["变量", "临床全变量_BH_FDR"],
    ]
    continuous_results = continuous_results.drop(
        columns=["BH_FDR", "FDR<0.05"]
    ).merge(continuous_fdr, on="变量", how="left")
    categorical_results = categorical_results.drop(
        columns=["BH_FDR", "FDR<0.05"]
    ).merge(categorical_fdr, on="变量", how="left")
    continuous_results = continuous_results.rename(
        columns={"临床全变量_BH_FDR": "BH_FDR"}
    )
    categorical_results = categorical_results.rename(
        columns={"临床全变量_BH_FDR": "BH_FDR"}
    )
    continuous_results["FDR<0.05"] = continuous_results["BH_FDR"].lt(0.05)
    categorical_results["FDR<0.05"] = categorical_results["BH_FDR"].lt(0.05)
    continuous_results = continuous_results.sort_values(
        ["BH_FDR", "P值", "变量"]
    )
    categorical_results = categorical_results.sort_values(
        ["BH_FDR", "P值", "变量"]
    )

    clinical_models = {}
    clinical_stability_sheets = {}
    for label, preoperative_only in [
        ("术前模型", True),
        ("探索性全临床模型", False),
    ]:
        matrix, base_by_column = encode_clinical(
            cleaned, preoperative_only=preoperative_only
        )
        search = tune_lasso(matrix, y, args.seed, args.n_jobs)
        transformed, base = stability_selection(
            matrix,
            y,
            base_by_column,
            search.best_params_,
            args.seed + (11 if preoperative_only else 29),
            args.stability_repetitions,
            args.n_jobs,
        )
        transformed.insert(0, "模型", label)
        base.insert(0, "模型", label)
        clinical_stability_sheets[label] = (transformed, base)
        clinical_models[label] = search.best_estimator_
        clinical_models[f"{label}_best_params"] = search.best_params_

    X_radiomics, case_ids, feature_names, radiomics_case_qc = load_radiomics(
        mapping, cache_path
    )
    imaging_batch_qc = imaging_batch_summary(
        radiomics_case_qc, mapping.reset_index(drop=True)
    )
    expected_case_ids = mapping["影像组学case_id"].astype(str).to_numpy()
    if not np.array_equal(case_ids, expected_case_ids):
        raise ValueError("Radiomics cache case order mismatch")
    if X_radiomics.shape != (len(clinical), 3190):
        raise ValueError(
            f"Expected radiomics shape {(len(clinical), 3190)}, "
            f"observed {X_radiomics.shape}"
        )
    radiomics_univ = radiomics_univariable(X_radiomics, y, feature_names)
    kept_indices, correlated_drops, nonconstant_indices = correlation_filter(
        X_radiomics,
        radiomics_univ,
        feature_names,
        args.correlation_threshold,
    )
    univ_rank = radiomics_univ.set_index("完整特征名")
    kept_names = feature_names[kept_indices]
    kept_order = sorted(
        kept_names,
        key=lambda name: (
            univ_rank.loc[name, "BH_FDR"],
            univ_rank.loc[name, "P值"],
            -abs(univ_rank.loc[name, "秩二列相关"]),
            name,
        ),
    )
    selected_for_lasso = np.asarray(
        kept_order[: min(args.radiomics_top_k, len(kept_order))]
    )
    name_to_index = {name: index for index, name in enumerate(feature_names)}
    selected_indices = np.asarray(
        [name_to_index[name] for name in selected_for_lasso]
    )
    radiomics_search, radiomics_stable = radiomics_stability(
        X_radiomics[:, selected_indices],
        y,
        selected_for_lasso,
        args.seed,
        args.stability_repetitions,
        args.n_jobs,
    )
    radiomics_stable = radiomics_stable.merge(
        radiomics_univ[
            [
                "完整特征名",
                "ROI",
                "PyRadiomics特征名",
                "秩二列相关",
                "P值",
                "BH_FDR",
            ]
        ],
        on="完整特征名",
        how="left",
    )
    radiomics_stable = radiomics_stable.sort_values(
        ["稳定选择频率", "BH_FDR", "完整特征名"],
        ascending=[False, True, True],
    )
    radiomics_prefilter = radiomics_univ.copy()
    radiomics_prefilter["非常量"] = radiomics_prefilter["完整特征名"].isin(
        feature_names[nonconstant_indices]
    )
    radiomics_prefilter["相关过滤后保留"] = radiomics_prefilter[
        "完整特征名"
    ].isin(kept_names)
    radiomics_prefilter["进入LASSO候选集"] = radiomics_prefilter[
        "完整特征名"
    ].isin(selected_for_lasso)

    clinical_base = pd.concat(
        [value[1] for value in clinical_stability_sheets.values()],
        ignore_index=True,
    )
    clinical_derived = pd.concat(
        [value[0] for value in clinical_stability_sheets.values()],
        ignore_index=True,
    )
    final_candidates = pd.concat(
        [
            clinical_base.loc[
                clinical_base["达到0.60阈值"],
                ["模型", "原始变量", "稳定选择频率"],
            ].rename(columns={"原始变量": "候选变量"}),
            radiomics_stable.loc[
                radiomics_stable["达到0.60阈值"],
                ["完整特征名", "稳定选择频率"],
            ]
            .rename(columns={"完整特征名": "候选变量"})
            .assign(模型="影像组学模型"),
        ],
        ignore_index=True,
    )
    final_candidates["说明"] = (
        "仅为当前混杂训练集的探索性候选；来源与结局完全重合，禁止直接进入"
        "最终模型或用于论文结论"
    )
    final_candidates["可直接进入最终模型"] = False

    summary = {
        "临床连续变量原始P<0.05数": str(
            int(continuous_results["原始P<0.05"].sum())
        ),
        "临床连续变量FDR<0.05数": str(
            int(continuous_results["FDR<0.05"].sum())
        ),
        "临床分类变量原始P<0.05数": str(
            int(categorical_results["原始P<0.05"].sum())
        ),
        "临床分类变量FDR<0.05数": str(
            int(categorical_results["FDR<0.05"].sum())
        ),
        "影像组学原始P<0.05数": str(
            int(radiomics_univ["原始P<0.05"].sum())
        ),
        "影像组学FDR<0.05数": str(
            int(radiomics_univ["FDR<0.05"].sum())
        ),
        "影像组学相关过滤后特征数": str(len(kept_indices)),
        "影像组学稳定频率>=0.60数": str(
            int(radiomics_stable["达到0.60阈值"].sum())
        ),
        "仅凭病例来源预测结局准确率": f"{source_accuracy:.3f}",
        "方法学状态": (
            "BLOCKED_来源与结局完全混杂"
            if source_accuracy == 1.0
            else "需进一步审计"
        ),
        "非漏组影像像素类型": str(
            imaging_batch_qc.loc[
                imaging_batch_qc["分组"].eq("非漏组"), "像素类型及例数"
            ].iloc[0]
        ),
        "漏组影像像素类型": str(
            imaging_batch_qc.loc[
                imaging_batch_qc["分组"].eq("漏组"), "像素类型及例数"
            ].iloc[0]
        ),
    }
    counts = {
        "total": int(len(y)),
        "negative": int((y == 0).sum()),
        "positive": int((y == 1).sum()),
    }
    notes = analysis_notes(args, counts, summary)

    with pd.ExcelWriter(workbook_path, engine="openpyxl") as writer:
        write_dataframe(writer, notes, "分析说明")
        write_dataframe(writer, leakage_audit, "关键泄漏与来源审计")
        write_dataframe(writer, clinical_qc, "临床_QC")
        write_dataframe(writer, continuous_results, "临床_连续变量差异")
        write_dataframe(writer, categorical_results, "临床_分类变量全局")
        write_dataframe(writer, categorical_levels, "临床_分类变量水平")
        write_dataframe(writer, clinical_base, "临床_稳定筛选_变量")
        write_dataframe(writer, clinical_derived, "临床_稳定筛选_派生项")
        write_dataframe(writer, radiomics_case_qc, "影像组学_病例QC")
        write_dataframe(writer, imaging_batch_qc, "影像组学_批次QC汇总")
        write_dataframe(writer, radiomics_prefilter, "影像组学_差异与初筛")
        write_dataframe(writer, correlated_drops, "影像组学_相关剔除")
        write_dataframe(writer, radiomics_stable, "影像组学_稳定筛选")
        write_dataframe(writer, final_candidates, "探索候选_禁止直接入模")
    style_workbook(workbook_path)

    joblib.dump(
        {
            "clinical_models": clinical_models,
            "radiomics_best_estimator": radiomics_search.best_estimator_,
            "radiomics_feature_names": selected_for_lasso,
            "radiomics_stability": radiomics_stable,
            "seed": args.seed,
        },
        model_path,
    )
    config = {
        "created_at": datetime.now().isoformat(timespec="seconds"),
        "input": str(input_path),
        "input_sha256": sha256(input_path),
        "output_workbook": str(workbook_path),
        "output_workbook_sha256": sha256(workbook_path),
        "cache": str(cache_path),
        "model": str(model_path),
        "seed": args.seed,
        "correlation_threshold": args.correlation_threshold,
        "radiomics_top_k": args.radiomics_top_k,
        "stability_repetitions": args.stability_repetitions,
        "class_counts": counts,
        "summary": summary,
        "source_only_outcome_accuracy": source_accuracy,
        "clinical_best_params": {
            key: value
            for key, value in clinical_models.items()
            if key.endswith("_best_params")
        },
        "radiomics_best_params": radiomics_search.best_params_,
        "software": {
            "python": __import__("sys").version.split()[0],
            "numpy": np.__version__,
            "pandas": pd.__version__,
            "scipy": __import__("scipy").__version__,
            "scikit_learn": __import__("sklearn").__version__,
            "statsmodels": sm.__version__,
            "SimpleITK": sitk.Version_VersionString(),
        },
    }
    config_path.write_text(
        json.dumps(config, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )
    print(json.dumps(config, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
