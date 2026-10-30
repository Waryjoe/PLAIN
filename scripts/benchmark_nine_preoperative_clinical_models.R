#!/usr/bin/env Rscript
source("scripts/plain_paths.R")

# Benchmark clinical prediction models from nine preoperative candidates.
# Model tuning and threshold selection use only the training set. The held-out
# validation set is read only after each model is locked.

options(stringsAsFactors = FALSE)
user_lib <- plain_path("R_LIBRARY", "")
if (dir.exists(user_lib)) .libPaths(c(user_lib, .libPaths()))
suppressPackageStartupMessages({
  library(readxl)
  library(openxlsx)
  library(caret)
  library(pROC)
})

seed <- 20260806
training_input <- plain_path("PROJECT", "汇总/训练验证集/训练集_70%_574例.xlsx")
validation_input <- plain_path("PROJECT", "汇总/训练验证集/验证集_30%_246例.xlsx")
vif_input <- paste0(
  plain_path("PROJECT", "汇总/"),
  "训练集优先候选变量_相关性与VIF_排除肿瘤部位/",
  "训练集_优先候选变量_VIF数值汇总_R.xlsx"
)
output_dir <- plain_path("PROJECT", "汇总/9个术前候选变量_临床模型比较")
output_file <- file.path(output_dir, "9个术前候选变量_内部验证模型性能比较.xlsx")
predictions_file <- file.path(output_dir, "9个术前候选变量_验证集预测概率.csv")
models_file <- file.path(output_dir, "9个术前候选变量_训练后模型对象_R.rds")
run_log_file <- file.path(output_dir, "9个术前候选变量_运行说明.txt")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# The local R installation uses the C locale. Copy Chinese-path workbooks to
# ASCII temporary paths, then read their first/English-named sheets.
temporary_training <- file.path(tempdir(), "training.xlsx")
temporary_validation <- file.path(tempdir(), "validation.xlsx")
temporary_vif <- file.path(tempdir(), "vif.xlsx")
temporary_output <- file.path(tempdir(), "clinical_model_benchmark.xlsx")
stopifnot(file.copy(training_input, temporary_training, overwrite = TRUE))
stopifnot(file.copy(validation_input, temporary_validation, overwrite = TRUE))
stopifnot(file.copy(vif_input, temporary_vif, overwrite = TRUE))

training_clinical <- as.data.frame(readxl::read_excel(temporary_training, sheet = 1))
validation_clinical <- as.data.frame(readxl::read_excel(temporary_validation, sheet = 1))
priority_candidates <- as.data.frame(readxl::read_excel(temporary_vif, sheet = "Priority_candidates"))
if (!all(c("variable", "variable_type", "p_value") %in% names(priority_candidates))) {
  stop("The VIF workbook does not contain the expected Priority_candidates columns")
}
priority_candidates <- priority_candidates[priority_candidates$p_value < 0.05, , drop = FALSE]
if (nrow(priority_candidates) != 9) stop("Expected exactly nine priority candidates")

continuous_tag <- "\u8fde\u7eed"
outcome_name <- "\u543b\u5408\u53e3"
selected_variables <- priority_candidates$variable
continuous_variables <- priority_candidates$variable[priority_candidates$variable_type == continuous_tag]
categorical_variables <- setdiff(selected_variables, continuous_variables)
if (length(continuous_variables) != 7 || length(categorical_variables) != 2) {
  stop("Unexpected continuous/categorical composition of the nine candidates")
}
if (!outcome_name %in% names(training_clinical) || !outcome_name %in% names(validation_clinical)) {
  stop("Outcome column is missing")
}
if (!all(selected_variables %in% names(training_clinical)) || !all(selected_variables %in% names(validation_clinical))) {
  stop("One or more candidate variables are missing from training or validation data")
}

make_model_data <- function(clinical) {
  output <- data.frame(outcome = ifelse(as.numeric(clinical[[outcome_name]]) == 1, "Leak", "NoLeak"))
  output$outcome <- factor(output$outcome, levels = c("Leak", "NoLeak"))
  for (index in seq_along(selected_variables)) {
    variable <- selected_variables[index]
    column <- paste0("x", index)
    if (variable %in% continuous_variables) {
      values <- suppressWarnings(as.numeric(clinical[[variable]]))
      values[values == 0] <- NA_real_
      output[[column]] <- values
    } else {
      output[[column]] <- factor(as.character(clinical[[variable]]))
    }
  }
  output
}

train_data <- make_model_data(training_clinical)
validation_data <- make_model_data(validation_clinical)
for (column in names(train_data)[-1]) {
  if (is.factor(train_data[[column]])) {
    validation_data[[column]] <- factor(validation_data[[column]], levels = levels(train_data[[column]]))
  }
}
if (any(is.na(validation_data$outcome)) || any(is.na(train_data$outcome))) stop("Outcome cannot be missing")

# Some classifiers (PLS-DA, SVM, and XGBoost) require numeric predictors.
# Explicit one-hot coding avoids model.matrix dropping rows with continuous NA.
# Levels are learned from training data only; continuous median imputation and
# scaling still happen separately within each resample.
categorical_columns <- names(train_data)[-1][vapply(train_data[-1], is.factor, logical(1))]
categorical_level_map <- lapply(categorical_columns, function(column) {
  sort(unique(c(as.character(train_data[[column]]), "Missing")))
})
names(categorical_level_map) <- categorical_columns

make_numeric_predictors <- function(model_data) {
  output <- list()
  for (column in names(model_data)[-1]) {
    if (!is.factor(model_data[[column]])) {
      output[[column]] <- as.numeric(model_data[[column]])
      next
    }
    values <- as.character(model_data[[column]])
    values[is.na(values)] <- "Missing"
    for (level in categorical_level_map[[column]]) {
      encoded_name <- paste0(column, "_", make.names(level))
      output[[encoded_name]] <- as.numeric(values == level)
    }
  }
  as.data.frame(output, check.names = FALSE)
}

train_predictors <- make_numeric_predictors(train_data)
validation_predictors <- make_numeric_predictors(validation_data)
if (!identical(names(train_predictors), names(validation_predictors)) || nrow(train_predictors) != nrow(train_data)) {
  stop("Numeric predictor construction failed")
}
nonconstant_columns <- vapply(train_predictors, function(column) {
  length(unique(column[!is.na(column)])) > 1
}, logical(1))
train_predictors <- train_predictors[, nonconstant_columns, drop = FALSE]
validation_predictors <- validation_predictors[, names(train_predictors), drop = FALSE]

variable_dictionary <- data.frame(
  Model_column = paste0("x", seq_along(selected_variables)),
  Original_variable = selected_variables,
  Type = ifelse(selected_variables %in% continuous_variables, "Continuous", "Categorical"),
  Preoperative_status = "Confirmed by investigator as first-visit/preoperative information"
)
missingness <- do.call(rbind, lapply(names(train_data)[-1], function(column) {
  data.frame(
    Model_column = column,
    Train_missing_n = sum(is.na(train_data[[column]])),
    Train_missing_rate = mean(is.na(train_data[[column]])),
    Validation_missing_n = sum(is.na(validation_data[[column]])),
    Validation_missing_rate = mean(is.na(validation_data[[column]]))
  )
}))

# Use identical repeated stratified folds for all models. All preprocessing is
# learned inside each resample by caret's preProcess pipeline.
set.seed(seed)
fold_index <- caret::createMultiFolds(train_data$outcome, k = 5, times = 5)
control <- caret::trainControl(
  method = "repeatedcv", number = 5, repeats = 5, index = fold_index,
  classProbs = TRUE, summaryFunction = caret::twoClassSummary,
  savePredictions = "final", allowParallel = FALSE
)
preprocess_steps <- c("medianImpute", "center", "scale")

lambda_grid <- exp(seq(log(0.0005), log(0.25), length.out = 5))
model_specs <- list(
  list(name = "Logistic", method = "glm", tune_grid = NULL, args = list()),
  list(name = "BPNN", method = "nnet", tune_grid = expand.grid(size = c(1, 3), decay = c(0, 0.05)), args = list(trace = FALSE, MaxNWts = 10000, maxit = 500)),
  list(name = "DT", method = "rpart", tune_grid = data.frame(cp = c(0.001, 0.01, 0.03, 0.05)), args = list()),
  list(name = "GBDT", method = "gbm", tune_grid = expand.grid(n.trees = c(50, 100), interaction.depth = c(1, 2), shrinkage = 0.05, n.minobsinnode = 10), args = list(verbose = FALSE)),
  list(name = "GNB", method = "nb", tune_grid = expand.grid(fL = c(0, 1), usekernel = c(TRUE, FALSE), adjust = 1), args = list()),
  list(name = "KNN", method = "knn", tune_grid = data.frame(k = c(3, 5, 7, 9, 11)), args = list()),
  list(name = "PLS-DA", method = "pls", tune_grid = data.frame(ncomp = 1:5), args = list()),
  list(name = "RF", method = "rf", tune_grid = data.frame(mtry = c(2, 4, 6)), args = list(ntree = 500)),
  list(name = "SVM", method = "svmRadial", tune_grid = expand.grid(sigma = c(0.03, 0.1), C = c(0.5, 2)), args = list())
)

find_youden_threshold <- function(outcome, probability) {
  candidates <- sort(unique(c(0, 1, probability)))
  scores <- vapply(candidates, function(threshold) {
    predicted <- ifelse(probability >= threshold, "Leak", "NoLeak")
    table_value <- table(factor(predicted, levels = c("Leak", "NoLeak")), outcome)
    sensitivity <- table_value["Leak", "Leak"] / sum(table_value[, "Leak"])
    specificity <- table_value["NoLeak", "NoLeak"] / sum(table_value[, "NoLeak"])
    sensitivity + specificity - 1
  }, numeric(1))
  candidates[which.max(scores)]
}

average_precision <- function(outcome, probability) {
  ordered <- order(probability, decreasing = TRUE)
  labels <- as.integer(outcome[ordered] == "Leak")
  positive_count <- sum(labels)
  if (positive_count == 0) return(NA_real_)
  precision <- cumsum(labels) / seq_along(labels)
  sum(precision[labels == 1]) / positive_count
}

calculate_validation_metrics <- function(outcome, probability, threshold) {
  predicted <- factor(ifelse(probability >= threshold, "Leak", "NoLeak"), levels = c("Leak", "NoLeak"))
  table_value <- table(predicted, outcome)
  true_positive <- table_value["Leak", "Leak"]
  false_positive <- table_value["Leak", "NoLeak"]
  false_negative <- table_value["NoLeak", "Leak"]
  true_negative <- table_value["NoLeak", "NoLeak"]
  sensitivity <- true_positive / (true_positive + false_negative)
  specificity <- true_negative / (true_negative + false_positive)
  ppv <- if ((true_positive + false_positive) == 0) NA_real_ else true_positive / (true_positive + false_positive)
  npv <- if ((true_negative + false_negative) == 0) NA_real_ else true_negative / (true_negative + false_negative)
  accuracy <- (true_positive + true_negative) / sum(table_value)
  f1 <- if (is.na(ppv) || (ppv + sensitivity) == 0) NA_real_ else 2 * ppv * sensitivity / (ppv + sensitivity)
  roc_value <- pROC::roc(response = outcome, predictor = probability, levels = c("NoLeak", "Leak"), direction = "<", quiet = TRUE)
  auc_value <- as.numeric(pROC::auc(roc_value))
  auc_ci <- as.numeric(pROC::ci.auc(roc_value, method = "delong"))
  list(
    AUC = auc_value, AUC_CI_low = auc_ci[1], AUC_CI_high = auc_ci[3],
    Accuracy = accuracy, F1 = f1, NPV = npv, PPV = ppv,
    Sensitivity = sensitivity, Specificity = specificity,
    PR_AUC = average_precision(outcome, probability),
    TP = true_positive, FP = false_positive, FN = false_negative, TN = true_negative,
    Predicted_class = predicted
  )
}

impute_and_scale <- function(training_x, new_x) {
  medians <- vapply(training_x, function(column) {
    value <- stats::median(column, na.rm = TRUE)
    if (is.finite(value)) value else 0
  }, numeric(1))
  training_x <- as.data.frame(Map(function(column, median_value) {
    column[is.na(column)] <- median_value
    column
  }, training_x, medians), check.names = FALSE)
  new_x <- as.data.frame(Map(function(column, median_value) {
    column[is.na(column)] <- median_value
    column
  }, new_x, medians), check.names = FALSE)
  names(training_x) <- names(new_x) <- names(medians)
  means <- vapply(training_x, mean, numeric(1))
  sds <- vapply(training_x, stats::sd, numeric(1))
  sds[!is.finite(sds) | sds == 0] <- 1
  list(
    training = sweep(sweep(as.matrix(training_x), 2, means, "-"), 2, sds, "/"),
    new = sweep(sweep(as.matrix(new_x), 2, means, "-"), 2, sds, "/")
  )
}

fold_auc <- function(outcome, probability) {
  as.numeric(pROC::auc(pROC::roc(outcome, probability, levels = c("NoLeak", "Leak"), direction = "<", quiet = TRUE)))
}

run_native_cv <- function(model_name, tuning_grid, fit_predict) {
  tuning_results <- vector("list", nrow(tuning_grid))
  for (grid_index in seq_len(nrow(tuning_grid))) {
    fold_scores <- numeric(length(fold_index))
    oof_outcome <- factor(character(), levels = c("Leak", "NoLeak"))
    oof_probability <- numeric()
    for (fold_position in seq_along(fold_index)) {
      train_index <- fold_index[[fold_position]]
      test_index <- setdiff(seq_len(nrow(train_predictors)), train_index)
      prepared <- impute_and_scale(train_predictors[train_index, , drop = FALSE], train_predictors[test_index, , drop = FALSE])
      probability <- fit_predict(
        prepared$training, train_data$outcome[train_index], prepared$new,
        tuning_grid[grid_index, , drop = FALSE], seed + grid_index * 1000 + fold_position
      )
      fold_scores[fold_position] <- fold_auc(train_data$outcome[test_index], probability)
      oof_outcome <- c(oof_outcome, train_data$outcome[test_index])
      oof_probability <- c(oof_probability, probability)
    }
    tuning_results[[grid_index]] <- list(
      ROC = mean(fold_scores), ROCSD = stats::sd(fold_scores),
      outcome = oof_outcome, probability = oof_probability
    )
  }
  scores <- vapply(tuning_results, `[[`, numeric(1), "ROC")
  best_index <- which.max(scores)
  list(best = tuning_grid[best_index, , drop = FALSE], result = tuning_results[[best_index]])
}

xgb_predict_probability <- function(training_x, training_y, new_x, parameters, run_seed) {
  set.seed(run_seed)
  model <- xgboost::xgb.train(
    params = list(
      objective = "binary:logistic", eval_metric = "auc", max_depth = parameters$max_depth,
      eta = parameters$eta, gamma = parameters$gamma, subsample = parameters$subsample,
      colsample_bytree = parameters$colsample_bytree, min_child_weight = parameters$min_child_weight,
      seed = run_seed, nthread = 1
    ),
    data = xgboost::xgb.DMatrix(data = training_x, label = as.integer(training_y == "Leak")),
    nrounds = parameters$nrounds, verbose = 0
  )
  as.numeric(predict(model, new_x))
}

lgb_predict_probability <- function(training_x, training_y, new_x, parameters, run_seed) {
  set.seed(run_seed)
  model <- lightgbm::lgb.train(
    params = list(
      objective = "binary", metric = "auc", learning_rate = parameters$learning_rate,
      num_leaves = parameters$num_leaves, min_data_in_leaf = parameters$min_data_in_leaf,
      feature_fraction = parameters$feature_fraction, bagging_fraction = parameters$bagging_fraction,
      bagging_freq = 1L, verbosity = -1L, seed = run_seed, num_threads = 1L
    ),
    data = lightgbm::lgb.Dataset(data = training_x, label = as.integer(training_y == "Leak")),
    nrounds = parameters$nrounds, verbose = -1
  )
  as.numeric(predict(model, new_x))
}

boruta_predict_probability <- function(training_x, training_y, new_x, parameters, run_seed) {
  set.seed(run_seed)
  boruta_fit <- Boruta::Boruta(
    x = as.data.frame(training_x), y = training_y, maxRuns = parameters$maxRuns, doTrace = 0
  )
  selected <- Boruta::getSelectedAttributes(boruta_fit, withTentative = FALSE)
  if (length(selected) == 0) {
    boruta_fit <- Boruta::TentativeRoughFix(boruta_fit)
    selected <- Boruta::getSelectedAttributes(boruta_fit, withTentative = FALSE)
  }
  if (length(selected) == 0) selected <- colnames(training_x)
  fitting_data <- data.frame(y = as.integer(training_y == "Leak"), training_x[, selected, drop = FALSE])
  logistic_fit <- stats::glm(y ~ ., data = fitting_data, family = stats::binomial())
  as.numeric(stats::predict(logistic_fit, newdata = as.data.frame(new_x[, selected, drop = FALSE]), type = "response"))
}

add_native_model <- function(model_name, method_name, tuning_grid, fit_predict, full_fit) {
  native_cv <- run_native_cv(model_name, tuning_grid, fit_predict)
  threshold <- find_youden_threshold(native_cv$result$outcome, native_cv$result$probability)
  final_prepared <- impute_and_scale(train_predictors, validation_predictors)
  final_object <- full_fit(final_prepared$training, train_data$outcome, native_cv$best, seed)
  probability <- final_object$predict(final_prepared$new)
  metrics <- calculate_validation_metrics(validation_data$outcome, probability, threshold)
  performance <<- rbind(performance, data.frame(
    Model = model_name, CV_ROC = native_cv$result$ROC, CV_ROC_SD = native_cv$result$ROCSD,
    Threshold_from_training = threshold,
    `AUC (95% CI)` = sprintf("%.3f (%.3f-%.3f)", metrics$AUC, metrics$AUC_CI_low, metrics$AUC_CI_high),
    Accuracy = metrics$Accuracy, `F1 score` = metrics$F1, NPV = metrics$NPV, PPV = metrics$PPV,
    Sensitivity = metrics$Sensitivity, Specificity = metrics$Specificity, PR_AUC = metrics$PR_AUC,
    TP = metrics$TP, FP = metrics$FP, FN = metrics$FN, TN = metrics$TN, check.names = FALSE
  ))
  tuning_rows[[model_name]] <<- data.frame(
    Model = model_name, Method = method_name, CV_ROC = native_cv$result$ROC,
    CV_ROC_SD = native_cv$result$ROCSD,
    Best_tuning_parameters = paste(paste0(names(native_cv$best), "=", as.character(native_cv$best[1, ])), collapse = "; "),
    stringsAsFactors = FALSE
  )
  validation_prediction_rows[[model_name]] <<- data.frame(
    Model = model_name, Outcome = validation_data$outcome, Probability_Leak = probability,
    Threshold_from_training = threshold, Predicted_class = metrics$Predicted_class
  )
  fits[[model_name]] <<- final_object$model
}

fits <- list()
performance <- data.frame(
  Model = character(), CV_ROC = numeric(), CV_ROC_SD = numeric(), Threshold_from_training = numeric(),
  `AUC (95% CI)` = character(), Accuracy = numeric(), `F1 score` = numeric(),
  NPV = numeric(), PPV = numeric(), Sensitivity = numeric(), Specificity = numeric(),
  PR_AUC = numeric(), TP = integer(), FP = integer(), FN = integer(), TN = integer(),
  stringsAsFactors = FALSE, check.names = FALSE
)
tuning_rows <- list()
validation_prediction_rows <- list()
failure_log <- data.frame(Model = character(), Method = character(), Error = character())

for (index in seq_along(model_specs)) {
  spec <- model_specs[[index]]
  cat(sprintf("Running %s (%d/%d)\\n", spec$name, index, length(model_specs)))
  flush.console()
  set.seed(seed + index)
  arguments <- c(
    list(
      x = train_predictors, y = train_data$outcome,
      method = spec$method, metric = "ROC",
      trControl = control, preProcess = preprocess_steps
    ),
    if (is.null(spec$tune_grid)) list(tuneLength = 1) else list(tuneGrid = spec$tune_grid),
    spec$args
  )
  fit <- tryCatch(do.call(caret::train, arguments), error = function(error) error)
  if (inherits(fit, "error")) {
    cat(sprintf("FAILED %s: %s\\n", spec$name, conditionMessage(fit)))
    flush.console()
    failure_log <- rbind(failure_log, data.frame(Model = spec$name, Method = spec$method, Error = conditionMessage(fit)))
    next
  }
  fits[[spec$name]] <- fit
  best_result <- fit$results[which.max(fit$results$ROC), , drop = FALSE]
  best_tune <- fit$bestTune
  tuning_rows[[spec$name]] <- data.frame(
    Model = spec$name,
    Method = spec$method,
    CV_ROC = best_result$ROC,
    CV_ROC_SD = best_result$ROCSD,
    Best_tuning_parameters = paste(
      paste0(names(best_tune), "=", as.character(best_tune[1, ])),
      collapse = "; "
    ),
    stringsAsFactors = FALSE
  )
  oof <- fit$pred
  probability_column <- "Leak"
  if (!probability_column %in% names(oof)) stop(paste0("Missing Leak probability for ", spec$name))
  threshold <- find_youden_threshold(oof$obs, oof[[probability_column]])
  probability <- predict(fit, validation_predictors, type = "prob")[[probability_column]]
  metrics <- calculate_validation_metrics(validation_data$outcome, probability, threshold)
  performance <- rbind(performance, data.frame(
    Model = spec$name,
    CV_ROC = best_result$ROC,
    CV_ROC_SD = best_result$ROCSD,
    Threshold_from_training = threshold,
    `AUC (95% CI)` = sprintf("%.3f (%.3f-%.3f)", metrics$AUC, metrics$AUC_CI_low, metrics$AUC_CI_high),
    Accuracy = metrics$Accuracy, `F1 score` = metrics$F1, NPV = metrics$NPV, PPV = metrics$PPV,
    Sensitivity = metrics$Sensitivity, Specificity = metrics$Specificity, PR_AUC = metrics$PR_AUC,
    TP = metrics$TP, FP = metrics$FP, FN = metrics$FN, TN = metrics$TN,
    check.names = FALSE
  ))
  validation_prediction_rows[[spec$name]] <- data.frame(
    Model = spec$name, Outcome = validation_data$outcome, Probability_Leak = probability,
    Threshold_from_training = threshold, Predicted_class = metrics$Predicted_class
  )
}

add_native_model(
  "Boruta", "Boruta + Logistic",
  data.frame(maxRuns = 50), boruta_predict_probability,
  function(training_x, training_y, parameters, run_seed) {
    set.seed(run_seed)
    boruta_fit <- Boruta::Boruta(x = as.data.frame(training_x), y = training_y, maxRuns = parameters$maxRuns, doTrace = 0)
    selected <- Boruta::getSelectedAttributes(boruta_fit, withTentative = FALSE)
    if (length(selected) == 0) {
      boruta_fit <- Boruta::TentativeRoughFix(boruta_fit)
      selected <- Boruta::getSelectedAttributes(boruta_fit, withTentative = FALSE)
    }
    if (length(selected) == 0) selected <- colnames(training_x)
    logistic_fit <- stats::glm(
      y ~ ., data = data.frame(y = as.integer(training_y == "Leak"), training_x[, selected, drop = FALSE]),
      family = stats::binomial()
    )
    list(
      model = list(boruta = boruta_fit, selected = selected, logistic = logistic_fit),
      predict = function(new_x) as.numeric(stats::predict(
        logistic_fit, newdata = as.data.frame(new_x[, selected, drop = FALSE]), type = "response"
      ))
    )
  }
)

add_native_model(
  "LightGBM", "lightgbm::lgb.train",
  expand.grid(nrounds = c(50, 100), num_leaves = c(3, 7), learning_rate = 0.05,
              min_data_in_leaf = 10, feature_fraction = 0.8, bagging_fraction = 0.8),
  lgb_predict_probability,
  function(training_x, training_y, parameters, run_seed) {
    set.seed(run_seed)
    model <- lightgbm::lgb.train(
      params = list(
        objective = "binary", metric = "auc", learning_rate = parameters$learning_rate,
        num_leaves = parameters$num_leaves, min_data_in_leaf = parameters$min_data_in_leaf,
        feature_fraction = parameters$feature_fraction, bagging_fraction = parameters$bagging_fraction,
        bagging_freq = 1L, verbosity = -1L, seed = run_seed, num_threads = 1L
      ),
      data = lightgbm::lgb.Dataset(data = training_x, label = as.integer(training_y == "Leak")),
      nrounds = parameters$nrounds, verbose = -1
    )
    list(model = model, predict = function(new_x) as.numeric(predict(model, new_x)))
  }
)

add_native_model(
  "XGBoost", "xgboost::xgb.train",
  expand.grid(nrounds = c(50, 100), max_depth = c(1, 2), eta = 0.05, gamma = 0,
              colsample_bytree = 0.8, min_child_weight = 1, subsample = 0.8),
  xgb_predict_probability,
  function(training_x, training_y, parameters, run_seed) {
    set.seed(run_seed)
    model <- xgboost::xgb.train(
      params = list(
        objective = "binary:logistic", eval_metric = "auc", max_depth = parameters$max_depth,
        eta = parameters$eta, gamma = parameters$gamma, subsample = parameters$subsample,
        colsample_bytree = parameters$colsample_bytree, min_child_weight = parameters$min_child_weight,
        seed = run_seed, nthread = 1
      ),
      data = xgboost::xgb.DMatrix(data = training_x, label = as.integer(training_y == "Leak")),
      nrounds = parameters$nrounds, verbose = 0
    )
    list(model = model, predict = function(new_x) as.numeric(predict(model, new_x)))
  }
)

if (nrow(performance) < 2) stop("Fewer than two models ran successfully")
final_model_order <- c("BPNN", "Boruta", "DT", "GBDT", "GNB", "KNN", "LightGBM", "Logistic", "PLS-DA", "RF", "SVM", "XGBoost")
performance <- performance[match(final_model_order, performance$Model, nomatch = 0), , drop = FALSE]
performance <- performance[!is.na(performance$Model), , drop = FALSE]
manuscript_table <- performance[, c("Model", "AUC (95% CI)", "Accuracy", "F1 score", "NPV", "PPV", "Sensitivity", "Specificity")]
names(manuscript_table)[1] <- "ML"
tuning_table <- do.call(rbind, tuning_rows)
validation_predictions <- do.call(rbind, validation_prediction_rows)
write.csv(validation_predictions, predictions_file, row.names = FALSE, fileEncoding = "UTF-8")
saveRDS(fits, models_file)

notes <- data.frame(
  Item = c(
    "Purpose", "Endpoint", "Training set", "Validation set", "Input variables", "Candidate rule", "Classifier set",
    "Timing", "Excluded variables", "Continuous zero rule", "Resampling", "Tuning metric",
    "Missing data", "Threshold", "Validation use", "AUC CI", "Performance boundary", "Boruta"
  ),
  Detail = c(
    "Clinical-model comparison using nine preoperative priority candidate variables",
    "Anastomotic leakage (Leak=1, NoLeak=0)",
    paste0(nrow(train_data), " patients; ", sum(train_data$outcome == "Leak"), " leakage events"),
    paste0(nrow(validation_data), " patients; ", sum(validation_data$outcome == "Leak"), " leakage events"),
    paste(selected_variables, collapse = "; "),
    "Univariable raw P<0.05 priority candidates only",
    "BPNN, Boruta + Logistic, DT, GBDT, GNB, KNN, LightGBM, Logistic, PLS-DA, RF, SVM, and XGBoost",
    "T stage and differentiation confirmed by investigator as first-visit/preoperative information",
    "Tumour location and TG2 excluded; no radiomics features included",
    "Exact zero in continuous measurements treated as missing; categorical zero codes retained",
    "Training-only repeated stratified 5-fold cross-validation (5 repeats), shared folds for all models",
    "Mean cross-validated ROC AUC",
    "Training median imputation for continuous variables within each resample; categorical predictors converted to training-level one-hot columns",
    "Youden threshold estimated from training out-of-fold probabilities; not optimised in validation",
    "Validation data not used for tuning, threshold selection, or choosing the final model",
    "Validation AUC 95% CI by DeLong method",
    "Validation has few leakage events; PPV, NPV, F1, and sensitivity are imprecise and should not alone select a final model",
    "Not included because Boruta is feature selection, not a classifier"
  )
)
package_versions <- data.frame(
  Package = c("R", "caret", "Boruta", "lightgbm", "randomForest", "e1071", "gbm", "klaR", "pROC", "pls", "xgboost", "readxl", "openxlsx"),
  Version = c(
    as.character(getRversion()),
    sapply(c("caret", "Boruta", "lightgbm", "randomForest", "e1071", "gbm", "klaR", "pROC", "pls", "xgboost", "readxl", "openxlsx"), function(package) as.character(packageVersion(package)))
  )
)

workbook <- openxlsx::createWorkbook()
header_style <- openxlsx::createStyle(fgFill = "#1F4E78", fontColour = "#FFFFFF", textDecoration = "bold", halign = "center", valign = "center", wrapText = TRUE)
number_style <- openxlsx::createStyle(numFmt = "0.000")
add_sheet <- function(name, data) {
  openxlsx::addWorksheet(workbook, name)
  openxlsx::writeData(workbook, name, data, withFilter = TRUE)
  openxlsx::addStyle(workbook, name, header_style, rows = 1, cols = seq_len(ncol(data)), gridExpand = TRUE)
  openxlsx::freezePane(workbook, name, firstRow = TRUE)
  openxlsx::setColWidths(workbook, name, cols = seq_len(ncol(data)), widths = "auto")
}
add_sheet("Validation_performance", performance)
openxlsx::addStyle(workbook, "Validation_performance", number_style, rows = 2:(nrow(performance) + 1), cols = c(2, 3, 4, 6:ncol(performance)), gridExpand = TRUE, stack = TRUE)
add_sheet("Manuscript_table", manuscript_table)
openxlsx::addStyle(workbook, "Manuscript_table", number_style, rows = 2:(nrow(manuscript_table) + 1), cols = 3:ncol(manuscript_table), gridExpand = TRUE, stack = TRUE)
add_sheet("CV_tuning_results", tuning_table)
add_sheet("Validation_predictions", validation_predictions)
add_sheet("Model_input_variables", variable_dictionary)
add_sheet("Missingness_QC", missingness)
add_sheet("Analysis_notes", notes)
add_sheet("Model_failures", failure_log)
add_sheet("Software_versions", package_versions)
openxlsx::saveWorkbook(workbook, temporary_output, overwrite = TRUE)
stopifnot(file.copy(temporary_output, output_file, overwrite = TRUE))

writeLines(c(
  paste0("output_file=", output_file),
  paste0("prediction_file=", predictions_file),
  paste0("models_file=", models_file),
  paste0("train_n=", nrow(train_data), "; train_events=", sum(train_data$outcome == "Leak")),
  paste0("validation_n=", nrow(validation_data), "; validation_events=", sum(validation_data$outcome == "Leak")),
  paste0("models_succeeded=", paste(performance$Model, collapse = "; "))
), run_log_file)
print(performance)
