#!/usr/bin/env Rscript
source("scripts/plain_paths.R")

# Final locked evaluation package for the selected 10-feature dual-ROI LightGBM
# radiomics model, Clinical LightGBM, and training-only logistic late fusion.
# Uses saved OOF/validation predictions only; no model is refitted on validation.

options(stringsAsFactors = FALSE)
user_lib <- plain_path("R_LIBRARY", "")
if (dir.exists(user_lib)) .libPaths(c(user_lib, .libPaths()))
suppressPackageStartupMessages({
  library(openxlsx)
  library(ggplot2)
  library(pROC)
  library(lightgbm)
})

project_root <- plain_path("PROJECT", "")
source_dir <- file.path(project_root, "汇总/胃肿瘤双ROI影像组学_v5_LASSO特征数_机器学习比较与临床融合")
output_dir <- file.path(project_root, "汇总/最终模型_10特征双ROI_LightGBM_全流程评估图表代码包_20260812")
tables_dir <- file.path(output_dir, "01_Tables")
figures_dir <- file.path(output_dir, "02_Figures")
code_dir <- file.path(output_dir, "03_R_Code")
notes_dir <- file.path(output_dir, "04_Notes")
for (directory in c(output_dir, tables_dir, figures_dir, code_dir, notes_dir)) dir.create(directory, recursive = TRUE, showWarnings = FALSE)

workbook_path <- file.path(source_dir, "双ROI影像组学_机器学习算法比较及临床融合_7_10_15_20_30项.xlsx")
model_rds <- file.path(source_dir, "影像模型_训练集全量拟合对象_R.rds")
external_predictions_path <- file.path(source_dir, "影像模型_全部性能.csv")
if (!all(file.exists(c(workbook_path, model_rds, external_predictions_path)))) stop("Required final-model inputs are missing")

rad_oof <- read.xlsx(workbook_path, sheet = "Radiomics_OOF_predictions")
rad_internal <- read.xlsx(workbook_path, sheet = "Radiomics_internal_predictions")
rad_external <- read.xlsx(workbook_path, sheet = "Radiomics_Ningbo_predictions")
fusion_predictions <- read.xlsx(workbook_path, sheet = "Fusion_predictions")
feature_sets <- read.xlsx(workbook_path, sheet = "Feature_sets")
all_performance <- read.csv(external_predictions_path, check.names = FALSE)

rad_oof <- rad_oof[rad_oof$Feature_count == 10 & rad_oof$Model == "LightGBM", ]
rad_internal <- rad_internal[rad_internal$Feature_count == 10 & rad_internal$Model == "LightGBM", ]
rad_external <- rad_external[rad_external$Feature_count == 10 & rad_external$Model == "LightGBM", ]
fusion_oof <- fusion_predictions[fusion_predictions$Feature_count == 10 & fusion_predictions$Radiomics_model == "LightGBM" & fusion_predictions$Dataset == "Training meta OOF", ]
fusion_internal <- fusion_predictions[fusion_predictions$Feature_count == 10 & fusion_predictions$Radiomics_model == "LightGBM" & fusion_predictions$Dataset == "Internal validation", ]
if (nrow(rad_oof) != 574L || nrow(rad_internal) != 246L || nrow(rad_external) != 353L || nrow(fusion_oof) != 574L || nrow(fusion_internal) != 246L) stop("Unexpected locked prediction dimensions")

make_prediction_set <- function(dataset, outcome, clinical, radiomics, fusion) {
  data.frame(
    Dataset = dataset,
    Outcome = as.integer(outcome),
    Clinical_LightGBM = as.numeric(clinical),
    Radiomics_10Feature_LightGBM = as.numeric(radiomics),
    Clinical_Radiomics_Fusion = as.numeric(fusion)
  )
}
training <- make_prediction_set("Training OOF", fusion_oof$Outcome, fusion_oof$Clinical_probability, fusion_oof$Radiomics_probability, fusion_oof$Fusion_probability)
internal <- make_prediction_set("Internal validation", fusion_internal$Outcome, fusion_internal$Clinical_probability, fusion_internal$Radiomics_probability, fusion_internal$Fusion_probability)
if (!identical(training$Outcome, as.integer(rad_oof$Outcome)) || !identical(internal$Outcome, as.integer(rad_internal$Outcome))) stop("Prediction outcome alignment failed")
external <- data.frame(Dataset = "Ningbo external imaging validation", Outcome = as.integer(rad_external$Outcome), Radiomics_10Feature_LightGBM = as.numeric(rad_external$Probability_Leak))

model_labels <- c(
  Clinical_LightGBM = "Clinical LightGBM",
  Radiomics_10Feature_LightGBM = "10-feature dual-ROI radiomics LightGBM",
  Clinical_Radiomics_Fusion = "Clinical-radiomics fusion"
)
model_colors <- c(
  "Clinical LightGBM" = "#3366CC",
  "10-feature dual-ROI radiomics LightGBM" = "#D95F02",
  "Clinical-radiomics fusion" = "#009E73"
)

safe_prob <- function(x) pmin(pmax(as.numeric(x), 1e-6), 1 - 1e-6)
youden_threshold <- function(y, p) {
  cutoffs <- sort(unique(c(0, 1, p)))
  statistic <- vapply(cutoffs, function(cut) {
    prediction <- as.integer(p >= cut)
    mean(prediction[y == 1L] == 1L) + mean(prediction[y == 0L] == 0L) - 1
  }, numeric(1))
  cutoffs[which.max(statistic)]
}
get_metrics <- function(y, p, threshold) {
  roc <- pROC::roc(y, p, levels = c(0, 1), direction = "<", quiet = TRUE)
  ci <- as.numeric(pROC::ci.auc(roc, method = "delong"))
  pred <- as.integer(p >= threshold)
  tp <- sum(pred == 1L & y == 1L); fp <- sum(pred == 1L & y == 0L)
  fn <- sum(pred == 0L & y == 1L); tn <- sum(pred == 0L & y == 0L)
  sensitivity <- tp / (tp + fn); specificity <- tn / (tn + fp)
  ppv <- ifelse(tp + fp == 0L, NA_real_, tp / (tp + fp)); npv <- ifelse(tn + fn == 0L, NA_real_, tn / (tn + fn))
  data.frame(AUC = as.numeric(pROC::auc(roc)), AUC_95CI_low = ci[1], AUC_95CI_high = ci[3],
    Threshold = threshold, Accuracy = mean(pred == y), Sensitivity = sensitivity, Specificity = specificity,
    PPV = ppv, NPV = npv, F1 = ifelse(ppv + sensitivity == 0, NA_real_, 2 * ppv * sensitivity / (ppv + sensitivity)),
    Brier_score = mean((p - y)^2), TP = tp, FP = fp, FN = fn, TN = tn)
}
get_calibration <- function(y, p) {
  p <- safe_prob(p)
  fit <- stats::glm(y ~ stats::qlogis(p), family = stats::binomial())
  data.frame(Calibration_intercept = unname(stats::coef(fit)[1]), Calibration_slope = unname(stats::coef(fit)[2]), Brier_score = mean((p - y)^2))
}
calibration_curve <- function(y, p, bins = 10L) {
  group <- cut(p, breaks = unique(quantile(p, probs = seq(0, 1, length.out = bins + 1L), type = 8)), include.lowest = TRUE)
  raw <- data.frame(Predicted = p, Observed = y, Bin = group)
  output <- aggregate(cbind(Predicted, Observed) ~ Bin, raw, mean)
  counts <- aggregate(Observed ~ Bin, raw, length); names(counts)[2] <- "N"
  merge(output, counts, by = "Bin", all.x = TRUE)
}
dca_curve <- function(y, p, thresholds = seq(0.01, 0.80, by = 0.01)) {
  n <- length(y); prevalence <- mean(y)
  data.frame(Threshold = thresholds,
    Model = vapply(thresholds, function(t) sum(p >= t & y == 1L) / n - sum(p >= t & y == 0L) / n * t / (1 - t), numeric(1)),
    Treat_all = prevalence - (1 - prevalence) * thresholds / (1 - thresholds),
    Treat_none = 0)
}

metric_rows <- list(); calibration_rows <- list(); confusion_rows <- list(); curve_rows <- list(); dca_rows <- list()
threshold_map <- list()
for (model in names(model_labels)) {
  train_probability <- training[[model]]; internal_probability <- internal[[model]]
  threshold <- youden_threshold(training$Outcome, train_probability)
  threshold_map[[model]] <- threshold
  for (item in list(list("Training OOF", training$Outcome, train_probability), list("Internal validation", internal$Outcome, internal_probability))) {
    metric_rows[[length(metric_rows) + 1L]] <- cbind(Model = model_labels[[model]], Dataset = item[[1]], get_metrics(item[[2]], item[[3]], threshold))
    calibration_rows[[length(calibration_rows) + 1L]] <- cbind(Model = model_labels[[model]], Dataset = item[[1]], get_calibration(item[[2]], item[[3]]))
    calibration_data <- calibration_curve(item[[2]], item[[3]])
    calibration_data$Model <- model_labels[[model]]; calibration_data$Dataset <- item[[1]]
    curve_rows[[length(curve_rows) + 1L]] <- calibration_data
    predicted <- as.integer(item[[3]] >= threshold)
    confusion_rows[[length(confusion_rows) + 1L]] <- data.frame(Model = model_labels[[model]], Dataset = item[[1]],
      Observed = rep(c("No leakage", "Leakage"), each = 2), Predicted = rep(c("No leakage", "Leakage"), 2),
      Count = c(sum(item[[2]] == 0L & predicted == 0L), sum(item[[2]] == 0L & predicted == 1L), sum(item[[2]] == 1L & predicted == 0L), sum(item[[2]] == 1L & predicted == 1L)))
    dca <- dca_curve(item[[2]], item[[3]])
    dca_rows[[length(dca_rows) + 1L]] <- rbind(
      data.frame(Dataset = item[[1]], Model = model_labels[[model]], Threshold = dca$Threshold, Net_benefit = dca$Model),
      data.frame(Dataset = item[[1]], Model = "Treat all", Threshold = dca$Threshold, Net_benefit = dca$Treat_all),
      data.frame(Dataset = item[[1]], Model = "Treat none", Threshold = dca$Threshold, Net_benefit = dca$Treat_none)
    )
  }
}
metrics <- do.call(rbind, metric_rows)
calibration_summary <- do.call(rbind, calibration_rows)
calibration_data <- do.call(rbind, curve_rows)
confusion_data <- do.call(rbind, confusion_rows)
dca_data <- do.call(rbind, dca_rows)

# External Ningbo analyses are radiomics-only; its 10 events make calibration and
# DCA explicitly exploratory, not a basis for model or feature selection.
external_threshold <- threshold_map[["Radiomics_10Feature_LightGBM"]]
external_metrics <- cbind(Model = "10-feature dual-ROI radiomics LightGBM", Dataset = "Ningbo external imaging validation", get_metrics(external$Outcome, external$Radiomics_10Feature_LightGBM, external_threshold))
external_calibration_summary <- cbind(Model = "10-feature dual-ROI radiomics LightGBM", Dataset = "Ningbo external imaging validation", get_calibration(external$Outcome, external$Radiomics_10Feature_LightGBM))
external_calibration_data <- calibration_curve(external$Outcome, external$Radiomics_10Feature_LightGBM)
external_calibration_data$Model <- "10-feature dual-ROI radiomics LightGBM"; external_calibration_data$Dataset <- "Ningbo external imaging validation"
external_dca_raw <- dca_curve(external$Outcome, external$Radiomics_10Feature_LightGBM)
external_dca_data <- rbind(
  data.frame(Dataset = "Ningbo external imaging validation", Model = "10-feature dual-ROI radiomics LightGBM", Threshold = external_dca_raw$Threshold, Net_benefit = external_dca_raw$Model),
  data.frame(Dataset = "Ningbo external imaging validation", Model = "Treat all", Threshold = external_dca_raw$Threshold, Net_benefit = external_dca_raw$Treat_all),
  data.frame(Dataset = "Ningbo external imaging validation", Model = "Treat none", Threshold = external_dca_raw$Threshold, Net_benefit = external_dca_raw$Treat_none)
)

# External Ningbo radiomics-only summary remains valid for the selected model.
external_summary <- all_performance[all_performance$Feature_count == 10 & all_performance$Model == "LightGBM" & all_performance$Dataset == "Ningbo external current", ]
external_summary$Model <- "10-feature dual-ROI radiomics LightGBM"

# Recover LightGBM feature importance from the locked full-training fit.
models <- readRDS(model_rds)
predict_env <- environment(models[["10__LightGBM"]]$predict)
lightgbm_object <- get("object", envir = predict_env)
importance <- as.data.frame(lightgbm::lgb.importance(lightgbm_object))
importance$Feature_original <- importance$Feature
importance$Feature_display <- sub("^stomach__", "Stomach: ", importance$Feature_original)
importance$Feature_display <- sub("^tumor__", "Tumor: ", importance$Feature_display)
importance$ROI <- ifelse(startsWith(importance$Feature_original, "stomach__"), "Stomach ROI", "Tumor ROI")
importance <- importance[order(-importance$Gain), ]
selected_features <- feature_sets[feature_sets$Feature_count == 10, ]

save_figure <- function(filename, plot, width, height) {
  temporary_file <- file.path(tempdir(), filename)
  ggsave(temporary_file, plot, width = width, height = height, dpi = 320)
  if (!file.copy(temporary_file, file.path(figures_dir, filename), overwrite = TRUE)) stop("Could not copy figure: ", filename)
}
theme_publication <- function() {
  theme_bw(base_size = 11) + theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
    legend.title = element_blank(), plot.title = element_text(face = "bold", hjust = 0), plot.subtitle = element_text(colour = "#444444"),
    axis.title = element_text(face = "bold"), strip.background = element_rect(fill = "#EEF3F8", colour = NA), strip.text = element_text(face = "bold"))
}
get_roc_data <- function(y, p, model, dataset) {
  roc <- pROC::roc(y, p, levels = c(0, 1), direction = "<", quiet = TRUE)
  data.frame(False_positive_rate = 1 - roc$specificities, True_positive_rate = roc$sensitivities, Model = model, Dataset = dataset,
    AUC = as.numeric(pROC::auc(roc)))
}
roc_data <- do.call(rbind, c(
  lapply(names(model_labels), function(model) get_roc_data(training$Outcome, training[[model]], model_labels[[model]], "Training OOF")),
  lapply(names(model_labels), function(model) get_roc_data(internal$Outcome, internal[[model]], model_labels[[model]], "Internal validation"))
))
roc_labels <- aggregate(AUC ~ Dataset + Model, roc_data, unique)
roc_annotations <- do.call(rbind, lapply(unique(roc_labels$Dataset), function(dataset) {
  block <- roc_labels[roc_labels$Dataset == dataset, ]
  display_order <- c("Clinical LightGBM", "10-feature dual-ROI radiomics LightGBM", "Clinical-radiomics fusion")
  block <- block[match(display_order, block$Model), ]
  data.frame(Dataset = dataset, x = 0.53, y = 0.31,
    Label = paste0(c("Clinical", "Radiomics", "Fusion"), ": ", sprintf("%.3f", block$AUC), collapse = "\n"))
}))

fig1 <- ggplot(roc_data, aes(x = False_positive_rate, y = True_positive_rate, colour = Model, group = Model)) +
  geom_abline(linetype = 2, colour = "#888888") + geom_path(linewidth = 0.9) + coord_equal() +
  geom_text(data = roc_annotations, aes(x = x, y = y, label = Label), inherit.aes = FALSE, hjust = 0, size = 3.5, lineheight = 0.95) +
  scale_colour_manual(values = model_colors) + facet_wrap(~Dataset) +
  labs(title = "ROC curves for the final selected models", x = "1 - Specificity", y = "Sensitivity") + theme_publication()
save_figure("Fig1_ROC_training_OOF_and_internal_validation.png", fig1, 11, 5.4)

fig2 <- ggplot(calibration_data, aes(x = Predicted, y = Observed, colour = Model)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "#777777") + geom_point(size = 1.8) + geom_line(linewidth = 0.7) +
  scale_colour_manual(values = model_colors) + coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) + facet_wrap(~Dataset) +
  labs(title = "Calibration curves", subtitle = "Points are deciles of predicted risk", x = "Mean predicted probability", y = "Observed leakage proportion") + theme_publication()
save_figure("Fig2_calibration_training_OOF_and_internal_validation.png", fig2, 11, 5.4)

dca_plot_data <- dca_data[dca_data$Threshold <= 0.50, ]
fig3 <- ggplot(dca_plot_data, aes(x = Threshold, y = Net_benefit, colour = Model, linetype = Model)) +
  geom_hline(yintercept = 0, colour = "#555555", linewidth = 0.35) + geom_line(linewidth = 0.85) +
  scale_colour_manual(values = c(model_colors, "Treat all" = "#7F7F7F", "Treat none" = "#1A1A1A")) +
  scale_linetype_manual(values = c("Clinical LightGBM" = "solid", "10-feature dual-ROI radiomics LightGBM" = "solid", "Clinical-radiomics fusion" = "solid", "Treat all" = "dashed", "Treat none" = "dotted")) +
  facet_wrap(~Dataset) + labs(title = "Decision curve analysis", subtitle = "Net benefit shown for threshold probabilities from 0.01 to 0.50", x = "Threshold probability", y = "Net benefit") + theme_publication()
save_figure("Fig3_DCA_training_OOF_and_internal_validation.png", fig3, 11, 5.4)

internal_metrics <- metrics[metrics$Dataset == "Internal validation", ]
metric_long <- rbind(
  data.frame(Model = internal_metrics$Model, Metric = "AUC", Value = internal_metrics$AUC),
  data.frame(Model = internal_metrics$Model, Metric = "Brier score", Value = internal_metrics$Brier_score),
  data.frame(Model = internal_metrics$Model, Metric = "Sensitivity", Value = internal_metrics$Sensitivity),
  data.frame(Model = internal_metrics$Model, Metric = "Specificity", Value = internal_metrics$Specificity)
)
fig4 <- ggplot(metric_long, aes(x = Model, y = Value, fill = Model)) + geom_col(width = 0.66) + facet_wrap(~Metric, scales = "free_y", nrow = 1) +
  scale_fill_manual(values = model_colors) + labs(title = "Internal validation performance summary", x = NULL, y = "Value") +
  theme_publication() + theme(axis.text.x = element_text(angle = 20, hjust = 1), legend.position = "none")
save_figure("Fig4_internal_validation_performance_summary.png", fig4, 12, 5.2)

importance$Feature_display <- factor(importance$Feature_display, levels = rev(importance$Feature_display))
fig5 <- ggplot(importance, aes(x = Feature_display, y = Gain, fill = ROI)) + geom_col(width = 0.70) + coord_flip() +
  scale_fill_manual(values = c("Stomach ROI" = "#4C78A8", "Tumor ROI" = "#F58518")) +
  labs(title = "Feature importance of the locked 10-feature radiomics LightGBM", x = NULL, y = "Relative gain") + theme_publication()
save_figure("Fig5_locked_10_feature_LightGBM_importance.png", fig5, 10, 6.3)

fig6 <- ggplot(confusion_data[confusion_data$Dataset == "Internal validation", ], aes(x = Predicted, y = Observed, fill = Count)) +
  geom_tile(colour = "white", linewidth = 1) + geom_text(aes(label = Count), size = 6, fontface = "bold") +
  scale_fill_gradient(low = "#EAF2F8", high = "#1F618D") + facet_wrap(~Model, nrow = 1) + coord_equal() +
  labs(title = "Confusion matrices in the internal validation cohort", x = "Predicted class", y = "Observed class", fill = "Patients") + theme_publication()
save_figure("Fig6_internal_validation_confusion_matrices.png", fig6, 12, 4.6)

# Predictive probability distribution for the locked 10-feature radiomics model.
risk_distribution <- rbind(
  data.frame(Dataset = "Training OOF", Outcome = factor(training$Outcome, levels = c(0, 1), labels = c("No leakage", "Leakage")), Probability = training$Radiomics_10Feature_LightGBM),
  data.frame(Dataset = "Internal validation", Outcome = factor(internal$Outcome, levels = c(0, 1), labels = c("No leakage", "Leakage")), Probability = internal$Radiomics_10Feature_LightGBM)
)
fig7 <- ggplot(risk_distribution, aes(x = Outcome, y = Probability, fill = Outcome)) + geom_boxplot(outlier.shape = NA, width = 0.62, alpha = 0.75) + geom_jitter(width = 0.10, alpha = 0.35, size = 1) +
  scale_fill_manual(values = c("No leakage" = "#8DA0CB", "Leakage" = "#FC8D62")) + facet_wrap(~Dataset) +
  labs(title = "Distribution of the 10-feature radiomics LightGBM probability", x = NULL, y = "Predicted leakage probability") + theme_publication() + theme(legend.position = "none")
save_figure("Fig7_radiomics_probability_distribution.png", fig7, 9, 5.2)

# Feature selection profile for the locked 10-feature set.
selected_features$Feature_display <- sub("^stomach__", "Stomach: ", selected_features$Feature)
selected_features$Feature_display <- sub("^tumor__", "Tumor: ", selected_features$Feature_display)
selected_features$ROI_display <- ifelse(selected_features$ROI == "Stomach", "Stomach ROI", "Tumor ROI")
selected_features$Feature_short <- selected_features$Feature
selected_features$Feature_short <- sub("^stomach__", "Stomach | ", selected_features$Feature_short)
selected_features$Feature_short <- sub("^tumor__", "Tumor | ", selected_features$Feature_short)
selected_features$Feature_short <- gsub("_", " | ", selected_features$Feature_short, fixed = TRUE)
selected_features$Feature_short <- factor(selected_features$Feature_short, levels = rev(selected_features$Feature_short[order(selected_features$Selection_frequency)]))
fig8 <- ggplot(selected_features, aes(x = Feature_short, y = Selection_frequency, fill = ROI_display)) + geom_col(width = 0.70) +
  scale_y_continuous(limits = c(0, 1.05), expand = expansion(mult = c(0, 0.02))) + coord_flip() +
  scale_fill_manual(values = c("Stomach ROI" = "#4C78A8", "Tumor ROI" = "#F58518")) +
  labs(title = "Stability selection profile of the locked 10 radiomics features", x = NULL, y = "LASSO selection frequency") + theme_publication()
save_figure("Fig8_stability_selection_locked_10_features.png", fig8, 10, 6.3)

external_roc <- get_roc_data(external$Outcome, external$Radiomics_10Feature_LightGBM, "10-feature dual-ROI radiomics LightGBM", "Ningbo external imaging validation")
external_auc <- unique(external_roc$AUC)
fig9 <- ggplot(external_roc, aes(x = False_positive_rate, y = True_positive_rate)) + geom_abline(linetype = 2, colour = "#888888") +
  geom_path(linewidth = 1.0, colour = "#D95F02") + coord_equal() +
  annotate("text", x = 0.56, y = 0.26, label = paste0("AUC = ", sprintf("%.3f", external_auc), "\n10 leakage events"), hjust = 0, size = 4) +
  labs(title = "Ningbo external imaging validation: radiomics ROC", subtitle = "Exploratory estimate with 10 leakage events", x = "1 - Specificity", y = "Sensitivity") + theme_publication()
save_figure("Fig9_Ningbo_external_radiomics_ROC.png", fig9, 7.2, 6.0)

fig10 <- ggplot(external_calibration_data, aes(x = Predicted, y = Observed)) + geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "#777777") +
  geom_point(colour = "#D95F02", size = 2) + geom_line(colour = "#D95F02", linewidth = 0.75) + coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
  labs(title = "Ningbo external radiomics calibration", subtitle = "Exploratory: 10 leakage events", x = "Mean predicted probability", y = "Observed leakage proportion") + theme_publication()
save_figure("Fig10_Ningbo_external_radiomics_calibration.png", fig10, 7.2, 6.0)

fig11 <- ggplot(external_dca_data[external_dca_data$Threshold <= 0.50, ], aes(x = Threshold, y = Net_benefit, colour = Model, linetype = Model)) +
  geom_hline(yintercept = 0, colour = "#555555", linewidth = 0.35) + geom_line(linewidth = 0.90) +
  scale_colour_manual(values = c("10-feature dual-ROI radiomics LightGBM" = "#D95F02", "Treat all" = "#7F7F7F", "Treat none" = "#1A1A1A")) +
  scale_linetype_manual(values = c("10-feature dual-ROI radiomics LightGBM" = "solid", "Treat all" = "dashed", "Treat none" = "dotted")) +
  labs(title = "Ningbo external radiomics decision curve", subtitle = "Exploratory: 10 leakage events", x = "Threshold probability", y = "Net benefit") + theme_publication()
save_figure("Fig11_Ningbo_external_radiomics_DCA.png", fig11, 8.0, 6.0)

write_sheet <- function(workbook, sheet, data) {
  addWorksheet(workbook, sheet); writeData(workbook, sheet, data, withFilter = TRUE)
  addStyle(workbook, sheet, createStyle(fgFill = "#1F4E78", fontColour = "#FFFFFF", textDecoration = "bold", halign = "center"), rows = 1, cols = seq_len(ncol(data)), gridExpand = TRUE)
  freezePane(workbook, sheet, firstRow = TRUE); setColWidths(workbook, sheet, cols = seq_len(ncol(data)), widths = "auto")
}
design <- data.frame(Item = c("Selected radiomics model", "Feature source", "Training cohort", "Internal validation cohort", "External imaging cohort", "Clinical model", "Fusion model", "Calibration", "DCA", "External fusion limitation"),
  Detail = c(
    "LightGBM using the top 10 training-only LASSO stability-ranked dual-ROI radiomics features (5 stomach and 5 tumor features).",
    "Feature ranking was fixed before this final evaluation; validation and Ningbo outcomes were not used for feature or algorithm selection.",
    "574 patients, 44 leakage events. Radiomics performance uses repeated 5x5-fold averaged OOF probabilities.",
    "246 patients, 19 leakage events. All prediction metrics/plots use predictions from models fitted on the training cohort.",
    "Ningbo imaging-only cohort: 353 patients, 10 leakage events. Only the radiomics model can be evaluated with currently available inputs.",
    "Previously locked Clinical LightGBM score based on the 9 preoperative candidate variables.",
    "Training-only logistic late fusion of Clinical LightGBM and 10-feature radiomics LightGBM probabilities.",
    "Calibration curves are grouped deciles; intercept/slope are apparent descriptive calibration estimates. Limited events make validation calibration uncertain.",
    "Net benefit calculated from saved probabilities. DCA supports comparison but does not establish clinical utility without prospective evaluation.",
    "No valid Ningbo clinical score for the full external cohort; no external clinical-radiomics fusion AUC, calibration or DCA is claimed."
  ))
summary_wb <- createWorkbook()
write_sheet(summary_wb, "Design_and_limits", design)
write_sheet(summary_wb, "Performance", metrics)
write_sheet(summary_wb, "Calibration_summary", calibration_summary)
write_sheet(summary_wb, "Calibration_curve_data", calibration_data)
write_sheet(summary_wb, "DCA_curve_data", dca_data)
write_sheet(summary_wb, "Ningbo_external_metrics", external_metrics)
write_sheet(summary_wb, "Ningbo_external_calibration", external_calibration_summary)
write_sheet(summary_wb, "Ningbo_external_cal_curve", external_calibration_data)
write_sheet(summary_wb, "Ningbo_external_DCA", external_dca_data)
write_sheet(summary_wb, "Internal_confusion_matrices", confusion_data[confusion_data$Dataset == "Internal validation", ])
write_sheet(summary_wb, "Selected_10_features", selected_features)
write_sheet(summary_wb, "LightGBM_importance", importance)
write_sheet(summary_wb, "Ningbo_radiomics_only", external_summary)
write_sheet(summary_wb, "Predictions", rbind(training, internal))
saveWorkbook(summary_wb, file.path(tables_dir, "Final_10Feature_DualROI_LightGBM_Evaluation.xlsx"), overwrite = TRUE)

file.copy(workbook_path, tables_dir, overwrite = TRUE)
file.copy(model_rds, tables_dir, overwrite = TRUE)
file.copy(file.path(source_dir, "影像模型_全部性能.csv"), tables_dir, overwrite = TRUE)
file.copy(file.path(source_dir, "临床影像融合模型_全部性能.csv"), tables_dir, overwrite = TRUE)

readme <- c(
  "# Final 10-feature dual-ROI LightGBM evaluation package", "",
  "## Contents", "- 01_Tables: final performance, calibration, DCA curve data, confusion matrices, selected features, importance, and saved prediction inputs.",
  "- 02_Figures: ROC, calibration, DCA, performance comparison, feature stability/importance, confusion matrices, risk-probability distribution, and separate Ningbo imaging-only evaluation figures.",
  "- 03_R_Code: this reproducible evaluation script plus the preceding model-selection script.",
  "- 04_Notes: environment, logs and a concise methods/limitations statement.", "",
  "## Locked final candidate", "The selected radiomics candidate is LightGBM using 10 dual-ROI features: 5 stomach ROI and 5 tumor ROI features.",
  "The clinical-radiomics model is logistic late fusion of the pre-existing Clinical LightGBM probability and the locked radiomics LightGBM probability.", "",
  "## Interpretation limits", "Training metrics use OOF probabilities; internal validation metrics use an untouched internal validation cohort.",
  "Ningbo currently supports imaging-only evaluation, not external clinical-radiomics fusion. The external cohort has only 10 leakage events; all external estimates are imprecise.",
  "Calibration curves and DCA are descriptive validation outputs and must be reported with the limited-event caveat."
)
writeLines(readme, file.path(notes_dir, "README.md"), useBytes = TRUE)
writeLines(capture.output(sessionInfo()), file.path(notes_dir, "R_sessionInfo.txt"))
this_script <- normalizePath(plain_path("REPO", "scripts/create_final_10feature_radiomics_evaluation_package.R"))
previous_script <- file.path(source_dir, "benchmark_v5_radiomics_ml_and_clinical_fusion.R")
file.copy(this_script, code_dir, overwrite = TRUE); file.copy(previous_script, code_dir, overwrite = TRUE)

message("Final evaluation package created: ", output_dir)
