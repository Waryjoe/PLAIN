#!/usr/bin/env Rscript
source("scripts/plain_paths.R")

# Benchmark the same 12 classifier families used for the nine-variable clinical
# model after fixed, training-only v5 dual-ROI LASSO stability feature ranking.
# Internal validation and Ningbo outcomes never take part in feature ranking,
# tuning, model selection, or fusion fitting.

options(stringsAsFactors = FALSE)
user_lib <- plain_path("R_LIBRARY", "")
if (dir.exists(user_lib)) .libPaths(c(user_lib, .libPaths()))
suppressPackageStartupMessages({
  library(readxl)
  library(openxlsx)
  library(caret)
  library(pROC)
  library(Boruta)
  library(lightgbm)
  library(xgboost)
})

seed <- 20260812L
set.seed(seed)
feature_counts <- c(7L, 10L, 15L, 20L, 30L)
model_order <- c("BPNN", "Boruta", "DT", "GBDT", "GNB", "KNN", "LightGBM", "Logistic", "PLS-DA", "RF", "SVM", "XGBoost")
project_root <- plain_path("PROJECT", "")
output_dir <- file.path(project_root, "汇总/胃肿瘤双ROI影像组学_v5_LASSO特征数_机器学习比较与临床融合")
code_archive_dir <- file.path(project_root, "代码与复现说明/scripts")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(code_archive_dir, recursive = TRUE, showWarnings = FALSE)

train_xlsx <- file.path(project_root, "汇总/影像组学_训练验证集_v5_GRAPE对齐/影像组学v5_训练集_70%_574例_3190特征.xlsx")
validation_xlsx <- file.path(project_root, "汇总/影像组学_训练验证集_v5_GRAPE对齐/影像组学v5_验证集_30%_246例_3190特征.xlsx")
stability_xlsx <- file.path(project_root, "汇总/胃肿瘤双ROI影像组学_v5_GRAPE对齐_训练集10折LASSO稳定性筛选/双ROI影像组学_训练集LASSO稳定性筛选.xlsx")
external_matrix_csv <- file.path(project_root, "原始数据备用/最终数据集/选定ct数据集/宁波/宁波汇总/外部验证/radiomics_locked_radscore_sampling/ningbo_external_dual_roi_radiomics_353x3190.csv.gz")
external_cohort_xlsx <- file.path(project_root, "原始数据备用/最终数据集/选定ct数据集/宁波/宁波汇总/匹配与QC/宁波外部验证_最终双ROI患者_353例.xlsx")
clinical_oof_csv <- file.path(project_root, "汇总/9个术前候选变量_训练集OOF模型比较。训练集/9个术前候选变量_训练集重复交叉验证OOF预测概率.csv")
clinical_validation_csv <- file.path(project_root, "汇总/9个术前候选变量_临床模型比较--内部验证集/9个术前候选变量_验证集预测概率.csv")

copy_to_temp <- function(path, suffix = basename(path)) {
  destination <- file.path(tempdir(), paste0(sprintf("%08d", sample.int(99999999L, 1L)), "_", suffix))
  if (!file.copy(path, destination, overwrite = TRUE)) stop("Could not copy input: ", path)
  destination
}

read_xlsx_safe <- function(path, sheet = 1) readxl::read_excel(copy_to_temp(path, "input.xlsx"), sheet = sheet)
clip_probability <- function(x) pmin(pmax(as.numeric(x), 1e-6), 1 - 1e-6)
logit_probability <- function(x) qlogis(clip_probability(x))

read_radiomics_matrix <- function(path, expected_n, expected_events) {
  data <- as.data.frame(read_xlsx_safe(path))
  feature_names <- names(data)[grepl("^(stomach|tumor)__", names(data))]
  if (nrow(data) != expected_n || length(feature_names) != 3190L ||
      sum(startsWith(feature_names, "stomach__")) != 1595L ||
      sum(startsWith(feature_names, "tumor__")) != 1595L) stop("Unexpected radiomics input schema")
  outcome <- as.integer(data[[3]])
  if (!identical(sort(unique(outcome)), c(0L, 1L)) || sum(outcome == 1L) != expected_events) stop("Unexpected outcome count")
  x <- as.matrix(data[, feature_names, drop = FALSE]); storage.mode(x) <- "double"
  if (any(!is.finite(x))) stop("Non-finite radiomics feature")
  list(hospital_id = as.character(data[[1]]), case_id = as.character(data[[2]]), outcome = outcome, x = x, features = feature_names)
}

make_foldid <- function(y, folds, local_seed) {
  set.seed(local_seed)
  output <- integer(length(y))
  for (value in sort(unique(y))) {
    indices <- sample(which(y == value))
    output[indices] <- rep(seq_len(folds), length.out = length(indices))
  }
  output
}

aggregate_repeated_oof <- function(row_index, outcome, probability, n_expected) {
  raw <- data.frame(RowIndex = row_index, Outcome = as.character(outcome), Probability_Leak = probability)
  averaged <- stats::aggregate(Probability_Leak ~ RowIndex + Outcome, data = raw, FUN = mean)
  averaged <- averaged[order(averaged$RowIndex), , drop = FALSE]
  observed <- sort(unique(as.integer(averaged$RowIndex)))
  missing_rows <- setdiff(seq_len(n_expected), observed)
  if (nrow(averaged) != n_expected || length(missing_rows) > 0L || anyDuplicated(observed)) {
    stop("OOF aggregation failed: expected ", n_expected, " unique rows; observed ", length(observed),
      "; missing rows: ", paste(utils::head(missing_rows, 10L), collapse = ","))
  }
  averaged$Outcome <- factor(averaged$Outcome, levels = c("Leak", "NoLeak"))
  averaged
}

youden_threshold <- function(y, probability) {
  candidates <- sort(unique(c(0, 1, probability)))
  scores <- vapply(candidates, function(threshold) {
    predicted <- ifelse(probability >= threshold, 1L, 0L)
    sensitivity <- sum(predicted[y == 1L] == 1L) / sum(y == 1L)
    specificity <- sum(predicted[y == 0L] == 0L) / sum(y == 0L)
    sensitivity + specificity - 1
  }, numeric(1))
  candidates[which.max(scores)]
}

auc_metrics <- function(y, probability, threshold = NULL) {
  if (is.null(threshold)) threshold <- youden_threshold(y, probability)
  roc <- pROC::roc(y, probability, levels = c(0, 1), direction = "<", quiet = TRUE)
  ci <- as.numeric(pROC::ci.auc(roc, method = "delong"))
  predicted <- as.integer(probability >= threshold)
  tp <- sum(predicted == 1L & y == 1L); fp <- sum(predicted == 1L & y == 0L)
  fn <- sum(predicted == 0L & y == 1L); tn <- sum(predicted == 0L & y == 0L)
  sensitivity <- tp / (tp + fn); specificity <- tn / (tn + fp)
  ppv <- ifelse(tp + fp == 0L, NA_real_, tp / (tp + fp))
  npv <- ifelse(tn + fn == 0L, NA_real_, tn / (tn + fn))
  f1 <- ifelse(is.na(ppv) || ppv + sensitivity == 0, NA_real_, 2 * ppv * sensitivity / (ppv + sensitivity))
  data.frame(
    AUC = as.numeric(pROC::auc(roc)), AUC_95CI_low = ci[1], AUC_95CI_high = ci[3],
    Threshold_from_training = threshold, Accuracy = mean(predicted == y), F1 = f1, PPV = ppv, NPV = npv,
    Sensitivity = sensitivity, Specificity = specificity, Brier_score = mean((probability - y)^2),
    TP = tp, FP = fp, FN = fn, TN = tn, check.names = FALSE
  )
}

prepare_numeric <- function(train_x, new_x) {
  means <- colMeans(train_x)
  sds <- apply(train_x, 2, stats::sd)
  sds[!is.finite(sds) | sds < 1e-10] <- 1
  list(train = sweep(sweep(as.matrix(train_x), 2, means, "-"), 2, sds, "/"),
       new = sweep(sweep(as.matrix(new_x), 2, means, "-"), 2, sds, "/"))
}

fold_auc <- function(y, probability) as.numeric(pROC::auc(pROC::roc(y, probability, levels = c("NoLeak", "Leak"), direction = "<", quiet = TRUE)))

run_native_cv <- function(x, y_factor, fold_index, tuning_grid, fit_predict, local_seed) {
  tuned <- vector("list", nrow(tuning_grid))
  for (grid_index in seq_len(nrow(tuning_grid))) {
    fold_scores <- numeric(length(fold_index)); all_probability <- numeric(); all_y <- factor(character(), levels = c("Leak", "NoLeak")); all_index <- integer()
    for (fold_position in seq_along(fold_index)) {
      development <- fold_index[[fold_position]]; heldout <- setdiff(seq_len(nrow(x)), development)
      prepared <- prepare_numeric(x[development, , drop = FALSE], x[heldout, , drop = FALSE])
      probability <- fit_predict(prepared$train, y_factor[development], prepared$new, tuning_grid[grid_index, , drop = FALSE], local_seed + grid_index * 1000L + fold_position)
      fold_scores[fold_position] <- fold_auc(y_factor[heldout], probability)
      all_probability <- c(all_probability, probability); all_y <- c(all_y, y_factor[heldout]); all_index <- c(all_index, heldout)
    }
    tuned[[grid_index]] <- list(ROC = mean(fold_scores), ROCSD = stats::sd(fold_scores), probability = all_probability, outcome = all_y, row_index = all_index)
  }
  values <- vapply(tuned, `[[`, numeric(1), "ROC"); best <- which.max(values)
  list(best = tuning_grid[best, , drop = FALSE], result = tuned[[best]])
}

xgb_predict <- function(training_x, training_y, new_x, parameters, local_seed) {
  set.seed(local_seed)
  model <- xgboost::xgb.train(
    params = list(objective = "binary:logistic", eval_metric = "auc", max_depth = parameters$max_depth,
                  eta = parameters$eta, gamma = parameters$gamma, subsample = parameters$subsample,
                  colsample_bytree = parameters$colsample_bytree, min_child_weight = parameters$min_child_weight,
                  seed = local_seed, nthread = 1L),
    data = xgboost::xgb.DMatrix(data = training_x, label = as.integer(training_y == "Leak")),
    nrounds = parameters$nrounds, verbose = 0
  )
  as.numeric(predict(model, new_x))
}

lgb_predict <- function(training_x, training_y, new_x, parameters, local_seed) {
  set.seed(local_seed)
  model <- lightgbm::lgb.train(
    params = list(objective = "binary", metric = "auc", learning_rate = parameters$learning_rate,
                  num_leaves = parameters$num_leaves, min_data_in_leaf = parameters$min_data_in_leaf,
                  feature_fraction = parameters$feature_fraction, bagging_fraction = parameters$bagging_fraction,
                  bagging_freq = 1L, verbosity = -1L, seed = local_seed, num_threads = 1L),
    data = lightgbm::lgb.Dataset(data = training_x, label = as.integer(training_y == "Leak")),
    nrounds = parameters$nrounds, verbose = -1
  )
  as.numeric(predict(model, new_x))
}

boruta_predict <- function(training_x, training_y, new_x, parameters, local_seed) {
  set.seed(local_seed)
  fit <- Boruta::Boruta(x = as.data.frame(training_x), y = training_y, maxRuns = parameters$maxRuns, doTrace = 0)
  selected <- Boruta::getSelectedAttributes(fit, withTentative = FALSE)
  if (!length(selected)) { fit <- Boruta::TentativeRoughFix(fit); selected <- Boruta::getSelectedAttributes(fit, withTentative = FALSE) }
  if (!length(selected)) selected <- colnames(training_x)
  logistic <- stats::glm(y ~ ., data = data.frame(y = as.integer(training_y == "Leak"), training_x[, selected, drop = FALSE], check.names = FALSE), family = stats::binomial())
  as.numeric(stats::predict(logistic, newdata = as.data.frame(new_x[, selected, drop = FALSE], check.names = FALSE), type = "response"))
}

fit_native_full <- function(model, x, y, parameter, local_seed) {
  prepared <- prepare_numeric(x, x)
  if (model == "LightGBM") {
    set.seed(local_seed)
    object <- lightgbm::lgb.train(params = list(objective = "binary", metric = "auc", learning_rate = parameter$learning_rate,
      num_leaves = parameter$num_leaves, min_data_in_leaf = parameter$min_data_in_leaf, feature_fraction = parameter$feature_fraction,
      bagging_fraction = parameter$bagging_fraction, bagging_freq = 1L, verbosity = -1L, seed = local_seed, num_threads = 1L),
      data = lightgbm::lgb.Dataset(data = prepared$train, label = as.integer(y == "Leak")), nrounds = parameter$nrounds, verbose = -1)
    return(list(predict = function(new_x) as.numeric(predict(object, prepare_numeric(x, new_x)$new))))
  }
  if (model == "XGBoost") {
    set.seed(local_seed)
    object <- xgboost::xgb.train(params = list(objective = "binary:logistic", eval_metric = "auc", max_depth = parameter$max_depth,
      eta = parameter$eta, gamma = parameter$gamma, subsample = parameter$subsample, colsample_bytree = parameter$colsample_bytree,
      min_child_weight = parameter$min_child_weight, seed = local_seed, nthread = 1L),
      data = xgboost::xgb.DMatrix(data = prepared$train, label = as.integer(y == "Leak")), nrounds = parameter$nrounds, verbose = 0)
    return(list(predict = function(new_x) as.numeric(predict(object, prepare_numeric(x, new_x)$new))))
  }
  set.seed(local_seed)
  object <- Boruta::Boruta(x = as.data.frame(prepared$train), y = y, maxRuns = parameter$maxRuns, doTrace = 0)
  selected <- Boruta::getSelectedAttributes(object, withTentative = FALSE)
  if (!length(selected)) { object <- Boruta::TentativeRoughFix(object); selected <- Boruta::getSelectedAttributes(object, withTentative = FALSE) }
  if (!length(selected)) selected <- colnames(prepared$train)
  logistic <- stats::glm(y ~ ., data = data.frame(y = as.integer(y == "Leak"), prepared$train[, selected, drop = FALSE], check.names = FALSE), family = stats::binomial())
  list(predict = function(new_x) as.numeric(stats::predict(logistic, newdata = as.data.frame(prepare_numeric(x, new_x)$new[, selected, drop = FALSE], check.names = FALSE), type = "response")))
}

train <- read_radiomics_matrix(train_xlsx, 574L, 44L)
validation <- read_radiomics_matrix(validation_xlsx, 246L, 19L)
if (!identical(train$features, validation$features)) stop("Internal feature schema mismatch")

stability <- as.data.frame(read_xlsx_safe(stability_xlsx, "Stability_selection"))
required_rank_columns <- c("Feature", "ROI", "Selection_frequency", "Direction_consistency", "P_value")
if (!all(required_rank_columns %in% names(stability))) stop("Stability table schema mismatch")
stability <- stability[stability$Feature %in% train$features, required_rank_columns, drop = FALSE]
stability <- stability[order(-stability$Selection_frequency, -stability$Direction_consistency, stability$P_value, stability$Feature), , drop = FALSE]
stability$Rank <- seq_len(nrow(stability))
if (nrow(stability) < max(feature_counts)) stop("Feature stability ranking too short")

external_matrix_temp <- copy_to_temp(external_matrix_csv, "external.csv.gz")
external <- utils::read.csv(gzfile(external_matrix_temp, "rt"), check.names = FALSE, fileEncoding = "UTF-8")
names(external) <- enc2utf8(names(external))
external_features <- names(external)[grepl("^(stomach|tumor)__", names(external))]
if (nrow(external) != 353L || length(external_features) != 3190L || !setequal(enc2utf8(train$features), enc2utf8(external_features))) stop("External radiomics schema mismatch")
external <- external[, c("case_id", enc2utf8(train$features)), drop = FALSE]
external_cohort_copy <- copy_to_temp(external_cohort_xlsx, "external_cohort.xlsx")
sheet_candidates <- readxl::excel_sheets(external_cohort_copy)
selected_sheet <- NULL
for (sheet in sheet_candidates) {
  preview <- as.data.frame(readxl::read_excel(external_cohort_copy, sheet = sheet, n_max = 1))
  if (all(c("case_id", "leak_label") %in% names(preview))) {
    body <- as.data.frame(readxl::read_excel(external_cohort_copy, sheet = sheet, col_types = "text"))
    if (nrow(body) == 353L) selected_sheet <- c(selected_sheet, sheet)
  }
}
if (length(selected_sheet) != 1L) stop("Could not identify external 353-case cohort sheet")
external_outcomes <- as.data.frame(readxl::read_excel(external_cohort_copy, sheet = selected_sheet, col_types = "text"))
external <- merge(external, data.frame(case_id = external_outcomes$case_id, outcome = as.integer(external_outcomes$leak_label)), by = "case_id", all.x = TRUE, sort = FALSE)
if (sum(external$outcome == 1L, na.rm = TRUE) != 10L || anyNA(external$outcome)) stop("External outcome matching failed")
external_x <- as.matrix(external[, enc2utf8(train$features), drop = FALSE]); storage.mode(external_x) <- "double"

# The same folds are reused across feature counts and algorithms, as in the
# clinical benchmark. They are tuning/OOF-development folds only.
y_train_factor <- factor(ifelse(train$outcome == 1L, "Leak", "NoLeak"), levels = c("Leak", "NoLeak"))
fold_index <- caret::createMultiFolds(y_train_factor, k = 5, times = 5)
control <- caret::trainControl(method = "repeatedcv", number = 5, repeats = 5, index = fold_index,
  classProbs = TRUE, summaryFunction = caret::twoClassSummary, savePredictions = "final", allowParallel = FALSE)
caret_specs <- list(
  list(name = "Logistic", method = "glm", tune = NULL, args = list()),
  list(name = "BPNN", method = "nnet", tune = expand.grid(size = c(1, 3), decay = c(0, 0.05)), args = list(trace = FALSE, MaxNWts = 10000, maxit = 500)),
  list(name = "DT", method = "rpart", tune = data.frame(cp = c(0.001, 0.01, 0.03, 0.05)), args = list()),
  list(name = "GBDT", method = "gbm", tune = expand.grid(n.trees = c(50, 100), interaction.depth = c(1, 2), shrinkage = 0.05, n.minobsinnode = 10), args = list(verbose = FALSE)),
  list(name = "GNB", method = "nb", tune = expand.grid(fL = c(0, 1), usekernel = c(TRUE, FALSE), adjust = 1), args = list()),
  list(name = "KNN", method = "knn", tune = data.frame(k = c(3, 5, 7, 9, 11)), args = list()),
  list(name = "PLS-DA", method = "pls", tune = data.frame(ncomp = 1:5), args = list()),
  list(name = "RF", method = "rf", tune = data.frame(mtry = c(2, 4, 6)), args = list(ntree = 500)),
  list(name = "SVM", method = "svmRadial", tune = expand.grid(sigma = c(0.03, 0.1), C = c(0.5, 2)), args = list())
)
native_specs <- list(
  list(name = "Boruta", method = "Boruta + Logistic", tune = data.frame(maxRuns = 50), predict = boruta_predict),
  list(name = "LightGBM", method = "lightgbm::lgb.train", tune = expand.grid(nrounds = c(50, 100), num_leaves = c(3, 7), learning_rate = 0.05, min_data_in_leaf = 10, feature_fraction = 0.8, bagging_fraction = 0.8), predict = lgb_predict),
  list(name = "XGBoost", method = "xgboost::xgb.train", tune = expand.grid(nrounds = c(50, 100), max_depth = c(1, 2), eta = 0.05, gamma = 0, colsample_bytree = 0.8, min_child_weight = 1, subsample = 0.8), predict = xgb_predict)
)

performance_rows <- list(); tuning_rows <- list(); oof_rows <- list(); internal_rows <- list(); external_rows <- list(); failure_rows <- list(); feature_rows <- list(); models <- list()
for (n_features in feature_counts) {
  features <- stability$Feature[seq_len(n_features)]
  scenario <- paste0("LASSO top ", n_features)
  feature_rows[[as.character(n_features)]] <- data.frame(Feature_count = n_features, Rank = seq_len(n_features), Feature = features,
    ROI = ifelse(startsWith(features, "stomach__"), "Stomach", "Tumor"), Selection_frequency = stability$Selection_frequency[seq_len(n_features)],
    Direction_consistency = stability$Direction_consistency[seq_len(n_features)], P_value = stability$P_value[seq_len(n_features)], stringsAsFactors = FALSE)
  x_train <- as.data.frame(train$x[, features, drop = FALSE]); x_validation <- as.data.frame(validation$x[, features, drop = FALSE]); x_external <- as.data.frame(external_x[, features, drop = FALSE])
  # Some R model interfaces alter '-' in PyRadiomics feature names. Use an
  # equivalent, unique internal schema consistently across train/validation/external.
  internal_feature_names <- make.names(features, unique = TRUE)
  colnames(x_train) <- internal_feature_names; colnames(x_validation) <- internal_feature_names; colnames(x_external) <- internal_feature_names
  for (spec_index in seq_along(caret_specs)) {
    spec <- caret_specs[[spec_index]]
    message("Running ", scenario, " / ", spec$name)
    arguments <- c(list(x = x_train, y = y_train_factor, method = spec$method, metric = "ROC", trControl = control, preProcess = c("center", "scale")),
                   if (is.null(spec$tune)) list(tuneLength = 1) else list(tuneGrid = spec$tune), spec$args)
    fit <- tryCatch(do.call(caret::train, arguments), error = function(e) e)
    if (inherits(fit, "error")) {
      failure_rows[[length(failure_rows) + 1L]] <- data.frame(Feature_count = n_features, Model = spec$name, Error = conditionMessage(fit)); next
    }
    # Compare tuning parameters using exact best-tune columns without relying on row order.
    matches <- rep(TRUE, nrow(fit$pred))
    for (column in names(fit$bestTune)) matches <- matches & fit$pred[[column]] == fit$bestTune[[column]][1]
    oof <- tryCatch(
      aggregate_repeated_oof(fit$pred$rowIndex[matches], fit$pred$obs[matches], fit$pred$Leak[matches], nrow(x_train)),
      error = function(e) e
    )
    if (inherits(oof, "error")) {
      failure_rows[[length(failure_rows) + 1L]] <- data.frame(Feature_count = n_features, Model = spec$name,
        Error = paste0("Incomplete or undefined repeated-CV OOF probabilities: ", conditionMessage(oof)))
      next
    }
    threshold <- youden_threshold(train$outcome, oof$Probability_Leak)
    p_validation <- as.numeric(predict(fit, newdata = x_validation, type = "prob")[, "Leak"])
    p_external <- as.numeric(predict(fit, newdata = x_external, type = "prob")[, "Leak"])
    performance_rows[[length(performance_rows) + 1L]] <- rbind(
      cbind(Feature_count = n_features, Stomach_features = sum(startsWith(features, "stomach__")), Tumor_features = sum(startsWith(features, "tumor__")), Model = spec$name, Dataset = "Training repeated 5x5 OOF", auc_metrics(train$outcome, oof$Probability_Leak, threshold)),
      cbind(Feature_count = n_features, Stomach_features = sum(startsWith(features, "stomach__")), Tumor_features = sum(startsWith(features, "tumor__")), Model = spec$name, Dataset = "Internal validation", auc_metrics(validation$outcome, p_validation, threshold)),
      cbind(Feature_count = n_features, Stomach_features = sum(startsWith(features, "stomach__")), Tumor_features = sum(startsWith(features, "tumor__")), Model = spec$name, Dataset = "Ningbo external current", auc_metrics(external$outcome, p_external, threshold))
    )
    tuning_rows[[length(tuning_rows) + 1L]] <- data.frame(Feature_count = n_features, Model = spec$name, Method = spec$method,
      CV_ROC = max(fit$results$ROC), CV_ROC_SD = fit$results$ROC[which.max(fit$results$ROC)] * 0 + fit$results$ROCSD[which.max(fit$results$ROC)],
      Best_tuning = paste(paste0(names(fit$bestTune), "=", as.character(fit$bestTune[1, ])), collapse = "; "))
    oof_rows[[length(oof_rows) + 1L]] <- data.frame(Feature_count = n_features, Model = spec$name, RowIndex = oof$RowIndex, Hospital_ID = train$hospital_id, Outcome = train$outcome, Probability_Leak = oof$Probability_Leak)
    internal_rows[[length(internal_rows) + 1L]] <- data.frame(Feature_count = n_features, Model = spec$name, Hospital_ID = validation$hospital_id, Outcome = validation$outcome, Probability_Leak = p_validation)
    external_rows[[length(external_rows) + 1L]] <- data.frame(Feature_count = n_features, Model = spec$name, Case_ID = external$case_id, Outcome = external$outcome, Probability_Leak = p_external)
    models[[paste(n_features, spec$name, sep = "__")]] <- fit
  }
  for (spec_index in seq_along(native_specs)) {
    spec <- native_specs[[spec_index]]
    message("Running ", scenario, " / ", spec$name)
    native <- tryCatch(run_native_cv(x_train, y_train_factor, fold_index, spec$tune, spec$predict, seed + n_features * 100L + spec_index * 10000L), error = function(e) e)
    if (inherits(native, "error")) {
      failure_rows[[length(failure_rows) + 1L]] <- data.frame(Feature_count = n_features, Model = spec$name, Error = conditionMessage(native)); next
    }
    oof <- aggregate_repeated_oof(native$result$row_index, native$result$outcome, native$result$probability, nrow(x_train))
    threshold <- youden_threshold(train$outcome, oof$Probability_Leak)
    fit <- fit_native_full(spec$name, x_train, y_train_factor, native$best, seed + n_features * 1000L + spec_index)
    p_validation <- fit$predict(x_validation); p_external <- fit$predict(x_external)
    performance_rows[[length(performance_rows) + 1L]] <- rbind(
      cbind(Feature_count = n_features, Stomach_features = sum(startsWith(features, "stomach__")), Tumor_features = sum(startsWith(features, "tumor__")), Model = spec$name, Dataset = "Training repeated 5x5 OOF", auc_metrics(train$outcome, oof$Probability_Leak, threshold)),
      cbind(Feature_count = n_features, Stomach_features = sum(startsWith(features, "stomach__")), Tumor_features = sum(startsWith(features, "tumor__")), Model = spec$name, Dataset = "Internal validation", auc_metrics(validation$outcome, p_validation, threshold)),
      cbind(Feature_count = n_features, Stomach_features = sum(startsWith(features, "stomach__")), Tumor_features = sum(startsWith(features, "tumor__")), Model = spec$name, Dataset = "Ningbo external current", auc_metrics(external$outcome, p_external, threshold))
    )
    tuning_rows[[length(tuning_rows) + 1L]] <- data.frame(Feature_count = n_features, Model = spec$name, Method = spec$method,
      CV_ROC = native$result$ROC, CV_ROC_SD = native$result$ROCSD, Best_tuning = paste(paste0(names(native$best), "=", as.character(native$best[1, ])), collapse = "; "))
    oof_rows[[length(oof_rows) + 1L]] <- data.frame(Feature_count = n_features, Model = spec$name, RowIndex = oof$RowIndex, Hospital_ID = train$hospital_id, Outcome = train$outcome, Probability_Leak = oof$Probability_Leak)
    internal_rows[[length(internal_rows) + 1L]] <- data.frame(Feature_count = n_features, Model = spec$name, Hospital_ID = validation$hospital_id, Outcome = validation$outcome, Probability_Leak = p_validation)
    external_rows[[length(external_rows) + 1L]] <- data.frame(Feature_count = n_features, Model = spec$name, Case_ID = external$case_id, Outcome = external$outcome, Probability_Leak = p_external)
    models[[paste(n_features, spec$name, sep = "__")]] <- fit
  }
}

performance <- do.call(rbind, performance_rows)
performance$Model <- factor(performance$Model, levels = model_order)
performance <- performance[order(performance$Feature_count, performance$Dataset, performance$Model), , drop = FALSE]
performance$Model <- as.character(performance$Model)
tuning <- do.call(rbind, tuning_rows); oof_predictions <- do.call(rbind, oof_rows); internal_predictions <- do.call(rbind, internal_rows); external_predictions <- do.call(rbind, external_rows)
selected_features <- do.call(rbind, feature_rows)
failures <- if (length(failure_rows)) do.call(rbind, failure_rows) else data.frame(Feature_count = integer(), Model = character(), Error = character())

# Clinical LightGBM probabilities are pre-existing locked predictions. Pair them
# with each radiomics-model probability through a logistic meta-model fitted only
# on OOF predictions. Ningbo lacks the complete clinical inputs, so no external
# fusion AUC is computed here.
clinical_oof <- utils::read.csv(copy_to_temp(clinical_oof_csv, "clinical_oof.csv"), check.names = FALSE, fileEncoding = "UTF-8")
clinical_validation <- utils::read.csv(copy_to_temp(clinical_validation_csv, "clinical_validation.csv"), check.names = FALSE, fileEncoding = "UTF-8")
clinical_oof <- clinical_oof[clinical_oof$Model == "LightGBM", , drop = FALSE]
clinical_validation <- clinical_validation[clinical_validation$Model == "LightGBM", , drop = FALSE]
clinical_oof <- clinical_oof[order(clinical_oof$RowIndex), , drop = FALSE]
if (nrow(clinical_oof) != 574L || nrow(clinical_validation) != 246L || !identical(as.integer(clinical_oof$Outcome == "Leak"), train$outcome) || !identical(as.integer(clinical_validation$Outcome == "Leak"), validation$outcome)) stop("Clinical prediction alignment failed")
clinical_train_probability <- clinical_oof$Probability_Leak; clinical_validation_probability <- clinical_validation$Probability_Leak
meta_foldid <- make_foldid(train$outcome, 10L, seed + 8000L)
fusion_rows <- list(); fusion_prediction_rows <- list()
for (n_features in feature_counts) {
  for (model in model_order) {
    r_oof <- oof_predictions[oof_predictions$Feature_count == n_features & oof_predictions$Model == model, , drop = FALSE]
    r_valid <- internal_predictions[internal_predictions$Feature_count == n_features & internal_predictions$Model == model, , drop = FALSE]
    if (nrow(r_oof) != 574L || nrow(r_valid) != 246L) next
    fusion_oof <- rep(NA_real_, 574L)
    for (fold in seq_len(10L)) {
      test <- which(meta_foldid == fold); development <- setdiff(seq_len(574L), test)
      fit <- stats::glm(outcome ~ clinical_logit + radiomics_logit, family = stats::binomial(),
        data = data.frame(outcome = train$outcome[development], clinical_logit = logit_probability(clinical_train_probability[development]), radiomics_logit = logit_probability(r_oof$Probability_Leak[development])))
      fusion_oof[test] <- stats::predict(fit, newdata = data.frame(clinical_logit = logit_probability(clinical_train_probability[test]), radiomics_logit = logit_probability(r_oof$Probability_Leak[test])), type = "response")
    }
    final_fit <- stats::glm(outcome ~ clinical_logit + radiomics_logit, family = stats::binomial(),
      data = data.frame(outcome = train$outcome, clinical_logit = logit_probability(clinical_train_probability), radiomics_logit = logit_probability(r_oof$Probability_Leak)))
    fusion_validation <- stats::predict(final_fit, newdata = data.frame(clinical_logit = logit_probability(clinical_validation_probability), radiomics_logit = logit_probability(r_valid$Probability_Leak)), type = "response")
    threshold <- youden_threshold(train$outcome, fusion_oof)
    fusion_name <- paste0("Clinical LightGBM + ", model)
    fusion_rows[[length(fusion_rows) + 1L]] <- rbind(
      cbind(Feature_count = n_features, Radiomics_model = model, Model = fusion_name, Dataset = "Training 10-fold meta OOF", auc_metrics(train$outcome, fusion_oof, threshold)),
      cbind(Feature_count = n_features, Radiomics_model = model, Model = fusion_name, Dataset = "Internal validation", auc_metrics(validation$outcome, fusion_validation, threshold))
    )
    fusion_prediction_rows[[length(fusion_prediction_rows) + 1L]] <- rbind(
      data.frame(Feature_count = n_features, Radiomics_model = model, Dataset = "Training meta OOF", Hospital_ID = train$hospital_id, Outcome = train$outcome, Clinical_probability = clinical_train_probability, Radiomics_probability = r_oof$Probability_Leak, Fusion_probability = fusion_oof),
      data.frame(Feature_count = n_features, Radiomics_model = model, Dataset = "Internal validation", Hospital_ID = validation$hospital_id, Outcome = validation$outcome, Clinical_probability = clinical_validation_probability, Radiomics_probability = r_valid$Probability_Leak, Fusion_probability = fusion_validation)
    )
  }
}
fusion_performance <- do.call(rbind, fusion_rows)
fusion_performance <- fusion_performance[order(fusion_performance$Feature_count, fusion_performance$Dataset, -fusion_performance$AUC), , drop = FALSE]
fusion_predictions <- do.call(rbind, fusion_prediction_rows)

internal_radiomics <- performance[performance$Dataset == "Internal validation", c("Feature_count", "Stomach_features", "Tumor_features", "Model", "AUC", "AUC_95CI_low", "AUC_95CI_high"), drop = FALSE]
internal_fusion <- fusion_performance[fusion_performance$Dataset == "Internal validation", c("Feature_count", "Radiomics_model", "Model", "AUC", "AUC_95CI_low", "AUC_95CI_high"), drop = FALSE]
best_radiomics <- do.call(rbind, lapply(feature_counts, function(n) internal_radiomics[internal_radiomics$Feature_count == n, , drop = FALSE][which.max(internal_radiomics$AUC[internal_radiomics$Feature_count == n]), , drop = FALSE]))
best_fusion <- do.call(rbind, lapply(feature_counts, function(n) internal_fusion[internal_fusion$Feature_count == n, , drop = FALSE][which.max(internal_fusion$AUC[internal_fusion$Feature_count == n]), , drop = FALSE]))

protocol <- data.frame(
  Item = c("Question", "Feature sets", "Feature-ranking boundary", "Radiomics algorithms", "Training estimate", "Internal validation", "Ningbo imaging-only evaluation", "Fusion", "External fusion boundary", "Model-selection caution"),
  Detail = c(
    "Does the post-LASSO radiomics classifier family change performance for 7/10/15/20/30 fixed ranked features?",
    "Top 7, 10, 15, 20 and 30 from the pre-existing training-only 100-resample LASSO stability ranking.",
    "Validation and Ningbo outcomes were excluded from ranking, feature count selection, tuning and final radiomics model fitting.",
    "Same 12 families as the clinical benchmark: BPNN, Boruta+Logistic, DT, GBDT, GNB, KNN, LightGBM, Logistic, PLS-DA, RF, SVM, XGBoost.",
    "Repeated stratified 5-fold CV with 5 repeats; OOF probabilities are averaged per patient. These estimates are conditional on the already locked feature sets.",
    "246 patients, 19 events. Used only after all feature sets and model hyperparameters were fitted in the 574-patient training cohort.",
    "353 patients, 10 events. Imaging-only model application only; confidence intervals are wide.",
    "Training-only logistic late fusion of the locked Clinical LightGBM probability plus each radiomics-model OOF probability; performance uses 10-fold meta-OOF and internal validation.",
    "No valid Ningbo clinical score is available for all 353 cases, so no external clinical-radiomics fusion AUC is reported.",
    "Do not choose a final model by the highest internal AUC across 60 exploratory combinations. Prioritize parsimony, repeated-OOF stability, internal validation, and future prespecified external evaluation."
  ), stringsAsFactors = FALSE)

workbook <- openxlsx::createWorkbook()
add_sheet <- function(name, data) {
  openxlsx::addWorksheet(workbook, name); openxlsx::writeData(workbook, name, data, withFilter = TRUE)
  header <- openxlsx::createStyle(fgFill = "#1F4E78", fontColour = "#FFFFFF", textDecoration = "bold", halign = "center")
  openxlsx::addStyle(workbook, name, header, rows = 1, cols = seq_len(ncol(data)), gridExpand = TRUE)
  openxlsx::freezePane(workbook, name, firstRow = TRUE); openxlsx::setColWidths(workbook, name, cols = seq_len(ncol(data)), widths = "auto")
}
add_sheet("Radiomics_performance", performance)
add_sheet("Best_internal_radiomics", best_radiomics)
add_sheet("Fusion_performance", fusion_performance)
add_sheet("Best_internal_fusion", best_fusion)
add_sheet("Feature_sets", selected_features)
add_sheet("Tuning", tuning)
add_sheet("Radiomics_OOF_predictions", oof_predictions)
add_sheet("Radiomics_internal_predictions", internal_predictions)
add_sheet("Radiomics_Ningbo_predictions", external_predictions)
add_sheet("Fusion_predictions", fusion_predictions)
add_sheet("Failures", failures)
add_sheet("Methods_and_limits", protocol)
xlsx_temp <- file.path(tempdir(), "radiomics_ml_benchmark.xlsx")
openxlsx::saveWorkbook(workbook, xlsx_temp, overwrite = TRUE)
file.copy(xlsx_temp, file.path(output_dir, "双ROI影像组学_机器学习算法比较及临床融合_7_10_15_20_30项.xlsx"), overwrite = TRUE)

write.csv(performance, file.path(output_dir, "影像模型_全部性能.csv"), row.names = FALSE, fileEncoding = "UTF-8")
write.csv(fusion_performance, file.path(output_dir, "临床影像融合模型_全部性能.csv"), row.names = FALSE, fileEncoding = "UTF-8")
write.csv(best_radiomics, file.path(output_dir, "每个特征数_内部验证最佳影像模型.csv"), row.names = FALSE, fileEncoding = "UTF-8")
write.csv(best_fusion, file.path(output_dir, "每个特征数_内部验证最佳融合模型.csv"), row.names = FALSE, fileEncoding = "UTF-8")

fmt <- function(x) ifelse(is.finite(x), sprintf("%.3f", x), "NA")
report <- c("# 双ROI影像组学机器学习算法比较与临床融合报告", "", "## 设计", "- 固定训练集 LASSO 稳定性排序前 7、10、15、20、30 项特征。", "- 每套特征比较与临床阶段一致的 12 种算法；训练使用重复5折交叉验证（5次重复）。", "- 融合使用已锁定 Clinical LightGBM 概率与影像模型概率的训练集内 logistic late fusion。", "- 宁波仅报告影像模型；当前无全队列锁定临床评分，不能报告外部融合 AUC。", "", "## 每个特征数的内部验证最佳模型")
for (n in feature_counts) {
  radiomics_row <- best_radiomics[best_radiomics$Feature_count == n, , drop = FALSE]
  fusion_row <- best_fusion[best_fusion$Feature_count == n, , drop = FALSE]
  external_row <- performance[performance$Feature_count == n & performance$Model == radiomics_row$Model & performance$Dataset == "Ningbo external current", , drop = FALSE]
  report <- c(report,
    paste0("### ", n, " 项（胃 ", radiomics_row$Stomach_features, "，肿瘤 ", radiomics_row$Tumor_features, "）"),
    paste0("- 最佳影像算法：", radiomics_row$Model, "；内部验证 AUC ", fmt(radiomics_row$AUC), "（95% CI ", fmt(radiomics_row$AUC_95CI_low), "-", fmt(radiomics_row$AUC_95CI_high), "）。"),
    paste0("- 同一影像算法宁波 AUC：", fmt(external_row$AUC), "（95% CI ", fmt(external_row$AUC_95CI_low), "-", fmt(external_row$AUC_95CI_high), "；10例漏）。"),
    paste0("- 最佳临床+影像融合：", fusion_row$Model, "；内部验证 AUC ", fmt(fusion_row$AUC), "（95% CI ", fmt(fusion_row$AUC_95CI_low), "-", fmt(fusion_row$AUC_95CI_high), "）。"), "")
}
report <- c(report, "## 解释边界", "- 本次为探索性模型家族比较；不要因 60 个组合中某个数值最高就替换已锁定的简约模型。", "- 内部验证仅有19例漏、宁波仅有10例漏，CI宽；待宁波补齐病例后，应一次性评价预先锁定的候选模型，不再依据外部AUC筛选。", "- 明细、特征来源、调参、OOF概率、内部与宁波影像预测和融合预测均保存于同目录 Excel。")
writeLines(report, file.path(output_dir, "机器学习比较与临床融合_中文总结.md"), useBytes = TRUE)
writeLines(capture.output(sessionInfo()), file.path(output_dir, "R_sessionInfo.txt"))
current_script <- normalizePath(plain_path("REPO", "scripts/benchmark_v5_radiomics_ml_and_clinical_fusion.R"))
file.copy(current_script, file.path(output_dir, basename(current_script)), overwrite = TRUE)
saveRDS(models, file.path(output_dir, "影像模型_训练集全量拟合对象_R.rds"))
print(best_radiomics); print(best_fusion)
