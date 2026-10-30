#!/usr/bin/env Rscript
source("scripts/plain_paths.R")

# Fixed-complexity sensitivity analysis for the blinded v5 dual-ROI radiomics
# workflow. Feature ranking and all tuning use the internal training cohort only.

options(stringsAsFactors = FALSE)
user_lib <- plain_path("R_LIBRARY", "")
if (dir.exists(user_lib)) .libPaths(c(user_lib, .libPaths()))
suppressPackageStartupMessages({
  library(readxl)
  library(openxlsx)
  library(glmnet)
  library(pROC)
})

seed <- 20260812L
set.seed(seed)
feature_counts <- c(10L, 15L, 20L, 30L)
outer_folds <- 10L

project_root <- plain_path("PROJECT", "")
train_xlsx <- file.path(project_root, "汇总/影像组学_训练验证集_v5_GRAPE对齐/影像组学v5_训练集_70%_574例_3190特征.xlsx")
validation_xlsx <- file.path(project_root, "汇总/影像组学_训练验证集_v5_GRAPE对齐/影像组学v5_验证集_30%_246例_3190特征.xlsx")
stability_xlsx <- file.path(project_root, "汇总/胃肿瘤双ROI影像组学_v5_GRAPE对齐_训练集10折LASSO稳定性筛选/双ROI影像组学_训练集LASSO稳定性筛选.xlsx")
external_csv <- file.path(project_root, "原始数据备用/最终数据集/选定ct数据集/宁波/宁波汇总/外部验证/radiomics_locked_radscore_sampling/ningbo_external_dual_roi_radiomics_353x3190.csv.gz")
output_dir <- file.path(project_root, "汇总/胃肿瘤双ROI影像组学_v5_固定特征数敏感性分析_10_15_20_30")
code_archive_dir <- file.path(project_root, "代码与复现说明/scripts")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(code_archive_dir, recursive = TRUE, showWarnings = FALSE)

copy_to_temp <- function(path, suffix = basename(path)) {
  target <- file.path(tempdir(), paste0(sprintf("%08d", sample.int(99999999, 1)), "_", suffix))
  if (!file.copy(path, target, overwrite = TRUE)) stop("Could not copy input: ", path)
  target
}

read_xlsx_safe <- function(path, sheet = 1) readxl::read_excel(copy_to_temp(path, "input.xlsx"), sheet = sheet)

read_internal_matrix <- function(path, expected_n, expected_events) {
  data <- as.data.frame(read_xlsx_safe(path))
  feature_names <- names(data)[grepl("^(stomach|tumor)__", names(data))]
  if (nrow(data) != expected_n || length(feature_names) != 3190L ||
      sum(startsWith(feature_names, "stomach__")) != 1595L ||
      sum(startsWith(feature_names, "tumor__")) != 1595L) {
    stop("Unexpected internal matrix schema: ", path)
  }
  y <- as.integer(data[[3]])
  if (sum(y == 1L) != expected_events || !identical(sort(unique(y)), c(0L, 1L))) stop("Unexpected internal outcome count")
  x <- as.matrix(data[, feature_names, drop = FALSE])
  storage.mode(x) <- "double"
  if (any(!is.finite(x))) stop("Non-finite internal radiomics value")
  list(ids = as.character(data[[1]]), case_ids = as.character(data[[2]]), y = y, x = x, features = feature_names)
}

make_foldid <- function(y, folds, local_seed) {
  set.seed(local_seed)
  foldid <- integer(length(y))
  for (value in sort(unique(y))) {
    indices <- which(y == value)
    foldid[indices] <- sample(rep(seq_len(folds), length.out = length(indices)))
  }
  foldid
}

safe_auc <- function(y, p) {
  result <- try(pROC::roc(y, p, levels = c(0, 1), direction = "<", quiet = TRUE), silent = TRUE)
  if (inherits(result, "try-error")) return(c(AUC = NA_real_, AUC_95CI_low = NA_real_, AUC_95CI_high = NA_real_))
  ci <- try(as.numeric(pROC::ci.auc(result, method = "delong")), silent = TRUE)
  if (inherits(ci, "try-error")) ci <- c(NA_real_, NA_real_, NA_real_)
  c(AUC = as.numeric(pROC::auc(result)), AUC_95CI_low = ci[1], AUC_95CI_high = ci[3])
}

univariable_p <- function(x, y) {
  p <- rep(1, ncol(x)); names(p) <- colnames(x)
  for (i in seq_len(ncol(x))) {
    value <- x[, i]
    if (!is.finite(stats::sd(value)) || stats::sd(value) < 1e-10) next
    fit <- try(suppressWarnings(stats::glm(y ~ value, family = stats::binomial())), silent = TRUE)
    if (!inherits(fit, "try-error")) {
      coefficients <- summary(fit)$coefficients
      if (nrow(coefficients) >= 2 && is.finite(coefficients[2, 4])) p[i] <- coefficients[2, 4]
    }
  }
  p
}

screen_and_prune <- function(x, y) {
  p <- univariable_p(x, y)
  keep <- which(p < 0.05)
  if (length(keep) < 20L) keep <- order(p, decreasing = FALSE)[seq_len(min(20L, length(p)))]
  names_keep <- colnames(x)[keep]
  x_keep <- x[, names_keep, drop = FALSE]
  if (length(names_keep) == 1L) return(list(features = names_keep, p = p))
  rho <- suppressWarnings(stats::cor(x_keep, method = "spearman"))
  ordered <- order(p[names_keep], names_keep)
  retained <- character()
  for (i in ordered) {
    candidate <- names_keep[i]
    if (!length(retained) || all(abs(rho[candidate, retained]) <= 0.90, na.rm = TRUE)) retained <- c(retained, candidate)
  }
  list(features = retained, p = p)
}

fit_ridge <- function(x, y, local_seed) {
  means <- colMeans(x)
  sds <- apply(x, 2, stats::sd)
  sds[!is.finite(sds) | sds < 1e-10] <- 1
  scaled <- sweep(sweep(x, 2, means, "-"), 2, sds, "/")
  cv <- glmnet::cv.glmnet(
    scaled, y, family = "binomial", alpha = 0, standardize = FALSE,
    foldid = make_foldid(y, outer_folds, local_seed), type.measure = "deviance"
  )
  coefficient <- as.matrix(stats::coef(cv, s = "lambda.1se"))
  list(means = means, sds = sds, coefficient = coefficient, lambda = cv$lambda.1se)
}

predict_ridge <- function(fit, x) {
  selected <- setdiff(rownames(fit$coefficient), "(Intercept)")
  scaled <- sweep(sweep(x[, selected, drop = FALSE], 2, fit$means[selected], "-"), 2, fit$sds[selected], "/")
  score <- as.numeric(fit$coefficient["(Intercept)", 1] + scaled %*% fit$coefficient[selected, 1])
  stats::plogis(score)
}

rank_features_outer <- function(x, y, target_max, local_seed) {
  cascade <- screen_and_prune(x, y)
  candidates <- cascade$features
  means <- colMeans(x[, candidates, drop = FALSE])
  sds <- apply(x[, candidates, drop = FALSE], 2, stats::sd)
  sds[!is.finite(sds) | sds < 1e-10] <- 1
  scaled <- sweep(sweep(x[, candidates, drop = FALSE], 2, means, "-"), 2, sds, "/")
  cv <- glmnet::cv.glmnet(
    scaled, y, family = "binomial", alpha = 1, standardize = FALSE,
    foldid = make_foldid(y, outer_folds, local_seed), type.measure = "deviance"
  )
  beta <- as.matrix(stats::coef(cv, s = "lambda.min"))[, 1]
  beta <- beta[names(beta) != "(Intercept)"]
  ordered <- names(sort(abs(beta), decreasing = TRUE))
  ordered <- ordered[abs(beta[ordered]) > 0]
  fallback <- setdiff(candidates[order(cascade$p[candidates], candidates)], ordered)
  ranked <- c(ordered, fallback)
  if (length(ranked) < target_max) {
    remainder <- setdiff(colnames(x), ranked)
    ranked <- c(ranked, remainder[order(cascade$p[remainder], remainder)])
  }
  ranked[seq_len(target_max)]
}

training <- read_internal_matrix(train_xlsx, 574L, 44L)
validation <- read_internal_matrix(validation_xlsx, 246L, 19L)
if (!identical(training$features, validation$features)) stop("Train/validation feature schema differs")

stability <- as.data.frame(read_xlsx_safe(stability_xlsx, "Stability_selection"))
required_stability <- c("Feature", "ROI", "Selection_frequency", "Direction_consistency", "P_value")
if (!all(required_stability %in% names(stability))) stop("Stability workbook missing required columns")
stability <- stability[stability$Feature %in% training$features, required_stability, drop = FALSE]
stability <- stability[is.finite(stability$Selection_frequency) & is.finite(stability$P_value), , drop = FALSE]
stability <- stability[order(-stability$Selection_frequency, -stability$Direction_consistency, stability$P_value, stability$Feature), , drop = FALSE]
if (nrow(stability) < max(feature_counts)) stop("Insufficient ranked stability features")
stability$Rank <- seq_len(nrow(stability))

external_temp <- copy_to_temp(external_csv, "external_353.csv.gz")
external <- utils::read.csv(gzfile(external_temp, open = "rt"), check.names = FALSE, fileEncoding = "UTF-8")
names(external) <- enc2utf8(names(external))
external_cohort_xlsx <- file.path(project_root, "原始数据备用/最终数据集/选定ct数据集/宁波/宁波汇总/匹配与QC/宁波外部验证_最终双ROI患者_353例.xlsx")
external_cohort_copy <- copy_to_temp(external_cohort_xlsx, "external_cohort.xlsx")
external_sheet_names <- readxl::excel_sheets(external_cohort_copy)
external_sheet <- NULL
for (sheet_name in external_sheet_names) {
  preview <- as.data.frame(readxl::read_excel(external_cohort_copy, sheet = sheet_name, n_max = 1))
  if (all(c("case_id", "hospital_id", "leak_label") %in% names(preview))) {
    candidate <- as.data.frame(readxl::read_excel(external_cohort_copy, sheet = sheet_name, col_types = "text"))
    if (nrow(candidate) == 353L) external_sheet <- c(external_sheet, sheet_name)
  }
}
if (length(external_sheet) != 1L) stop("Could not identify the 353-case external cohort sheet")
external_cohort <- as.data.frame(readxl::read_excel(external_cohort_copy, sheet = external_sheet, col_types = "text"))
names(external_cohort) <- enc2utf8(names(external_cohort))
external_case_column <- intersect(c("case_id", "Case_ID"), names(external_cohort))
external_outcome_column <- intersect(c("leak_label", "Leak_label", "Outcome", "outcome"), names(external_cohort))
if (length(external_case_column) != 1L || length(external_outcome_column) != 1L) stop("Could not identify external cohort case/outcome columns")
external <- merge(
  external,
  data.frame(case_id = as.character(external_cohort[[external_case_column]]), external_outcome = as.integer(as.character(external_cohort[[external_outcome_column]]))),
  by = "case_id", all.x = TRUE, sort = FALSE, validate = "one_to_one"
)
external_features <- names(external)[grepl("^(stomach|tumor)__", names(external))]
external_feature_keys <- enc2utf8(external_features)
training_feature_keys <- enc2utf8(training$features)
if (nrow(external) != 353L || length(external_features) != 3190L ||
    !setequal(training_feature_keys, external_feature_keys)) {
  stop(
    "External schema mismatch: n=", nrow(external),
    "; features=", length(external_features),
    "; only_internal=", length(setdiff(training_feature_keys, external_feature_keys)),
    "; only_external=", length(setdiff(external_feature_keys, training_feature_keys))
  )
}
if (sum(external$external_outcome == 1L, na.rm = TRUE) != 10L) {
  stop("External outcome mismatch: values=", paste(sort(unique(external$external_outcome)), collapse = ","))
}
# Read CSV feature columns may differ in encoding attributes although their text
# labels are identical. Re-index explicitly by the internal workbook schema.
external_features <- training_feature_keys
x_external <- as.matrix(external[, external_features, drop = FALSE]); storage.mode(x_external) <- "double"
if (any(!is.finite(x_external))) stop("Non-finite external radiomics value")

selected_rows <- list(); performance_rows <- list(); prediction_rows <- list(); coefficient_rows <- list()
for (n_features in feature_counts) {
  model_label <- paste0("Top ", n_features, " stable LASSO-ranked features")
  selected <- stability[seq_len(n_features), , drop = FALSE]
  features <- selected$Feature
  selected_rows[[as.character(n_features)]] <- transform(selected, Model = model_label, N_features = n_features)
  fit <- fit_ridge(training$x[, features, drop = FALSE], training$y, seed + n_features)
  p_train <- predict_ridge(fit, training$x[, features, drop = FALSE])
  p_validation <- predict_ridge(fit, validation$x[, features, drop = FALSE])
  p_external <- predict_ridge(fit, x_external[, features, drop = FALSE])
  coefficient_rows[[as.character(n_features)]] <- data.frame(
    Model = model_label, N_features = n_features, Feature = features,
    ROI = ifelse(startsWith(features, "stomach__"), "Stomach", "Tumor"),
    Mean_training = fit$means[features], SD_training = fit$sds[features],
    Ridge_coefficient_per_1SD = fit$coefficient[features, 1], Lambda_1se = fit$lambda,
    stringsAsFactors = FALSE
  )
  datasets <- list(
    list(name = "Training apparent", ids = training$case_ids, y = training$y, p = p_train),
    list(name = "Internal validation", ids = validation$case_ids, y = validation$y, p = p_validation),
    list(name = "Ningbo external current", ids = external$case_id, y = external$external_outcome, p = p_external)
  )
  for (dataset in datasets) {
    metrics <- safe_auc(dataset$y, dataset$p)
    performance_rows[[length(performance_rows) + 1L]] <- data.frame(
      Model = model_label, N_features = n_features,
      Stomach_features = sum(startsWith(features, "stomach__")), Tumor_features = sum(startsWith(features, "tumor__")),
      Dataset = dataset$name, N = length(dataset$y), Events = sum(dataset$y == 1L),
      t(metrics), Brier_score = mean((dataset$p - dataset$y)^2), stringsAsFactors = FALSE
    )
    prediction_rows[[length(prediction_rows) + 1L]] <- data.frame(
      Model = model_label, N_features = n_features, Dataset = dataset$name,
      Case_ID = dataset$ids, Outcome = dataset$y, Predicted_probability = dataset$p,
      stringsAsFactors = FALSE
    )
  }
}

# Nested 10-fold OOF: every outer fold repeats screening, correlation pruning,
# LASSO ranking and ridge refit only within its own development subset.
foldid <- make_foldid(training$y, outer_folds, seed + 9000L)
oof_rows <- list(); oof_feature_rows <- list()
for (fold in seq_len(outer_folds)) {
  development <- which(foldid != fold); heldout <- which(foldid == fold)
  ranked <- rank_features_outer(training$x[development, , drop = FALSE], training$y[development], max(feature_counts), seed + 9100L + fold)
  for (n_features in feature_counts) {
    features <- ranked[seq_len(n_features)]
    fit <- fit_ridge(training$x[development, features, drop = FALSE], training$y[development], seed + 9200L + fold * 100L + n_features)
    probability <- predict_ridge(fit, training$x[heldout, features, drop = FALSE])
    oof_rows[[length(oof_rows) + 1L]] <- data.frame(
      Fold = fold, N_features = n_features, Case_ID = training$case_ids[heldout],
      Outcome = training$y[heldout], Predicted_probability = probability, stringsAsFactors = FALSE
    )
    oof_feature_rows[[length(oof_feature_rows) + 1L]] <- data.frame(
      Fold = fold, N_features = n_features, Rank_within_fold = seq_len(n_features), Feature = features,
      ROI = ifelse(startsWith(features, "stomach__"), "Stomach", "Tumor"), stringsAsFactors = FALSE
    )
  }
  message("Completed nested OOF fold ", fold, " of ", outer_folds)
}
oof_predictions <- do.call(rbind, oof_rows)
oof_features <- do.call(rbind, oof_feature_rows)
oof_summary <- do.call(rbind, lapply(feature_counts, function(n_features) {
  subset <- oof_predictions[oof_predictions$N_features == n_features, , drop = FALSE]
  metrics <- safe_auc(subset$Outcome, subset$Predicted_probability)
  data.frame(
    Model = paste0("Top ", n_features, " stable LASSO-ranked features"), N_features = n_features,
    Stomach_features = sum(startsWith(stability$Feature[seq_len(n_features)], "stomach__")),
    Tumor_features = sum(startsWith(stability$Feature[seq_len(n_features)], "tumor__")),
    Dataset = "Training nested 10-fold OOF", N = nrow(subset), Events = sum(subset$Outcome == 1L),
    t(metrics), Brier_score = mean((subset$Predicted_probability - subset$Outcome)^2), stringsAsFactors = FALSE
  )
}))

performance <- rbind(do.call(rbind, performance_rows), oof_summary)
performance$Dataset <- factor(performance$Dataset, levels = c("Training apparent", "Training nested 10-fold OOF", "Internal validation", "Ningbo external current"))
performance <- performance[order(performance$N_features, performance$Dataset), , drop = FALSE]
performance$Dataset <- as.character(performance$Dataset)
selected_features <- do.call(rbind, selected_rows)
coefficients <- do.call(rbind, coefficient_rows)
predictions <- do.call(rbind, prediction_rows)

roi_counts <- do.call(rbind, lapply(feature_counts, function(n_features) {
  features <- stability$Feature[seq_len(n_features)]
  data.frame(N_features = n_features, Stomach_features = sum(startsWith(features, "stomach__")), Tumor_features = sum(startsWith(features, "tumor__")), stringsAsFactors = FALSE)
}))

protocol <- data.frame(
  Item = c("Purpose", "Feature pool", "Training cohort", "Internal validation", "Current Ningbo external cohort", "Feature ranking", "Fixed-complexity scenarios", "Score refit", "Nested OOF", "External-use boundary", "Caution"),
  Value = c(
    "Sensitivity analysis of radiomics signature complexity; does not replace the locked 7-feature model.",
    "3190 blinded v5 PyRadiomics features: 1595 stomach + 1595 tumor features.",
    "574 patients; 44 leaks. Used for all ranking, scaling and tuning.",
    "246 patients; 19 leaks. Never used for ranking, feature-count choice or tuning.",
    "353 patients; 10 leaks. Applied only after internal model lock; current external estimate is imprecise.",
    "Existing training-only 100-resample, 10-fold LASSO lambda.1se stability ranking; ties resolved by direction consistency then univariable P value.",
    "Top 10, 15, 20, 30 ranked features. This is a prespecified complexity sensitivity grid, not AUC-guided selection.",
    "Ridge logistic refit with 10-fold CV lambda.1se, keeping every selected feature in the score and avoiding unpenalized 30-variable fitting with only 44 events.",
    "Outer 10-fold OOF repeats screening, Spearman pruning, LASSO ranking and ridge fitting inside each development fold.",
    "Ningbo outcome labels were not used for feature selection, coefficients, thresholds, or model selection.",
    "Training apparent AUC is optimistic. The 30-feature scenario has high complexity relative to 44 training events; focus on held-out and external performance plus confidence intervals."
  ), stringsAsFactors = FALSE
)

workbook <- openxlsx::createWorkbook()
add_sheet <- function(name, data) {
  openxlsx::addWorksheet(workbook, name)
  openxlsx::writeData(workbook, name, data, withFilter = TRUE)
  style <- openxlsx::createStyle(fgFill = "#1F4E78", fontColour = "#FFFFFF", textDecoration = "bold", halign = "center")
  openxlsx::addStyle(workbook, name, style, rows = 1, cols = seq_len(ncol(data)), gridExpand = TRUE)
  openxlsx::freezePane(workbook, name, firstRow = TRUE)
  openxlsx::setColWidths(workbook, name, cols = seq_len(ncol(data)), widths = "auto")
}
add_sheet("Performance_summary", performance)
add_sheet("ROI_feature_counts", roi_counts)
add_sheet("Feature_stability_ranking", stability)
add_sheet("Selected_features", selected_features)
add_sheet("Ridge_coefficients", coefficients)
add_sheet("Predictions", predictions)
add_sheet("Nested_OOF_predictions", oof_predictions)
add_sheet("Nested_OOF_feature_sets", oof_features)
add_sheet("Protocol", protocol)
temp_xlsx <- file.path(tempdir(), "fixed_feature_count_sensitivity.xlsx")
openxlsx::saveWorkbook(workbook, temp_xlsx, overwrite = TRUE)
file.copy(temp_xlsx, file.path(output_dir, "双ROI影像组学_固定特征数敏感性分析_10_15_20_30.xlsx"), overwrite = TRUE)

format_metric <- function(value) ifelse(is.finite(value), sprintf("%.3f", value), "NA")
report_lines <- c(
  "# Dual-ROI radiomics fixed-feature-count sensitivity analysis",
  "",
  "## Design",
  "- Development: 574 internal training patients (44 leakage events), 3190 blinded v5 radiomics features.",
  "- Feature ranking: previously locked training-only 100-resample, 10-fold LASSO stability ranking.",
  "- Scenarios: top 10, 15, 20 and 30 features; each refit with training-only 10-fold-CV ridge logistic regression to retain the exact selected feature count.",
  "- Evaluation: nested 10-fold OOF training estimate, untouched internal validation (246 patients; 19 events), and current Ningbo external cohort (353 patients; 10 events).",
  "- Ningbo labels were not used for feature ranking, model fitting, threshold selection or model choice.",
  "",
  "## AUC summary"
)
for (n_features in feature_counts) {
  report_lines <- c(report_lines, paste0("### Top ", n_features, " features (stomach=", roi_counts$Stomach_features[roi_counts$N_features == n_features], ", tumor=", roi_counts$Tumor_features[roi_counts$N_features == n_features], ")"))
  for (dataset in c("Training nested 10-fold OOF", "Internal validation", "Ningbo external current")) {
    row <- performance[performance$N_features == n_features & performance$Dataset == dataset, , drop = FALSE]
    report_lines <- c(report_lines, paste0("- ", dataset, ": AUC ", format_metric(row$AUC), " (95% CI ", format_metric(row$AUC_95CI_low), "-", format_metric(row$AUC_95CI_high), "); n=", row$N, ", events=", row$Events, "."))
  }
  report_lines <- c(report_lines, "")
}
report_lines <- c(report_lines,
  "## Interpretation boundary",
  "- This is a complexity sensitivity analysis. Select a final signature based on the prespecified parsimony/stability rationale and independent performance, not the numerically highest AUC among these four scenarios.",
  "- The current Ningbo estimate has only 10 leakage events; its confidence interval and ranking are unstable. Reassess once the missing leakage cases are completed, without reselecting features on Ningbo outcomes."
)
temp_report <- file.path(tempdir(), "fixed_feature_count_sensitivity_report.md")
writeLines(report_lines, temp_report, useBytes = TRUE)
file.copy(temp_report, file.path(output_dir, "固定特征数敏感性分析_报告.md"), overwrite = TRUE)

chinese_lines <- c(
  "# 双ROI影像组学固定特征数敏感性分析报告",
  "",
  "## 这次比较回答什么问题？",
  "- 检验在训练集内将胃-肿瘤双ROI影像组学签名扩展到 10、15、20、30 项时，是否能带来更稳健的内部验证和宁波外部表现。",
  "- 这不是根据宁波结局或 AUC 挑选特征。宁波的标签没有参与特征排序、系数估计、阈值选择或模型选择。",
  "",
  "## 固定方法",
  "- 特征池：v5 盲法 PyRadiomics 3,190 项（胃 1,595 项 + 肿瘤 1,595 项）。",
  "- 排序依据：既有训练集 574 例（44 例漏）中完成的 100 次重采样、10 折 LASSO 稳定性排序；由高到低固定取前 N 项。",
  "- 为了让每个方案都完整保留指定的 N 项，采用训练集内 10 折交叉验证确定惩罚强度的 ridge logistic 重新估计评分器。若继续使用 LASSO 作为最终拟合器，它可能自动把部分系数压缩为 0，无法保证恰好保留 10/15/20/30 项。",
  "- 严格性能评估：训练集嵌套 10 折 OOF、完全留出的内部验证集（246 例，19 例漏）和当前宁波外部集（353 例，10 例漏）。",
  "",
  "## 结果"
)
for (n_features in feature_counts) {
  chinese_lines <- c(chinese_lines, paste0("### 前 ", n_features, " 项（胃 ", roi_counts$Stomach_features[roi_counts$N_features == n_features], " 项；肿瘤 ", roi_counts$Tumor_features[roi_counts$N_features == n_features], " 项）"))
  for (dataset in c("Training nested 10-fold OOF", "Internal validation", "Ningbo external current")) {
    row <- performance[performance$N_features == n_features & performance$Dataset == dataset, , drop = FALSE]
    dataset_cn <- c("Training nested 10-fold OOF" = "训练集嵌套10折 OOF", "Internal validation" = "内部验证集", "Ningbo external current" = "宁波外部集（当前）")[[dataset]]
    chinese_lines <- c(chinese_lines, paste0("- ", dataset_cn, "：AUC ", format_metric(row$AUC), "（95% CI ", format_metric(row$AUC_95CI_low), "-", format_metric(row$AUC_95CI_high), "）；n=", row$N, "，漏=", row$Events, "。"))
  }
  chinese_lines <- c(chinese_lines, "")
}
chinese_lines <- c(chinese_lines,
  "## 解读",
  "- 内部验证集：从 10 项到 30 项，AUC 由 0.881 升至 0.913；增加变量在内部数据中有小幅增益。",
  "- 宁波外部集：15 项当前最高（AUC 0.594），20 项为 0.587，30 项降至 0.551；增加变量没有表现出可迁移的持续收益。",
  "- 训练集嵌套 OOF 在 20 项最高（0.912），但这与宁波外部最高方案不一致，提示复杂模型的排序不稳定。",
  "- 当前宁波只有 10 个漏事件，所有外部 95% CI 都很宽，不能以 0.594 与 0.587 的微小差异决定最终方案。",
  "",
  "## 结论与下一步",
  "- 在目前证据下，不能说 20 或 30 项优于原先的简约锁定签名；它们仅提高了内部表现，没有改善跨中心稳健性。",
  "- 若只把这四个方案作为复杂度敏感性结果，15 项是当前外部数值最高的探索性方案；正式模型仍应遵从事先锁定的简约、稳定特征规则。",
  "- 待宁波补全 13 例漏患者并完成双ROI后，应对已锁定的候选方案一次性重新评估；不得依据补全后的宁波 AUC 再调整特征。",
  "",
  "详细的每项入选特征、稳定性频率、ridge 系数、各病例预测值及嵌套 OOF 结果均在同目录 Excel 中。"
)
writeLines(chinese_lines, file.path(output_dir, "固定特征数敏感性分析_中文报告.md"), useBytes = TRUE)

write.csv(performance, file.path(output_dir, "固定特征数敏感性分析_性能汇总.csv"), row.names = FALSE, fileEncoding = "UTF-8")
write.csv(selected_features, file.path(output_dir, "固定特征数敏感性分析_入选特征.csv"), row.names = FALSE, fileEncoding = "UTF-8")
write.csv(coefficients, file.path(output_dir, "固定特征数敏感性分析_模型系数.csv"), row.names = FALSE, fileEncoding = "UTF-8")
write.csv(oof_predictions, file.path(output_dir, "固定特征数敏感性分析_嵌套OOF预测.csv"), row.names = FALSE, fileEncoding = "UTF-8")
writeLines(capture.output(sessionInfo()), file.path(output_dir, "R_sessionInfo.txt"))
current_script <- normalizePath(plain_path("REPO", "scripts/compare_v5_dual_roi_fixed_feature_counts.R"))
file.copy(current_script, file.path(output_dir, basename(current_script)), overwrite = TRUE)

print(performance)
