#!/usr/bin/env Rscript
source("scripts/plain_paths.R")

# Publication figures for the locked clinical LightGBM model.
# Performance uses training OOF and the sealed internal validation cohort.
# SHAP explains the final LightGBM fitted on all training patients; no re-training,
# re-selection, or validation-driven tuning occurs here.

options(stringsAsFactors = FALSE)
user_lib <- plain_path("R_LIBRARY", "")
if (dir.exists(user_lib)) .libPaths(c(user_lib, .libPaths()))
suppressPackageStartupMessages({
  library(readxl)
  library(openxlsx)
  library(pROC)
  library(lightgbm)
  library(ggplot2)
  library(patchwork)
  library(svglite)
  library(ragg)
  library(scales)
})

seed <- 20260806
training_input <- plain_path("PROJECT", "汇总/训练验证集/训练集_70%_574例.xlsx")
vif_input <- plain_path("PROJECT", "汇总/训练集优先候选变量_相关性与VIF_排除肿瘤部位/训练集_优先候选变量_VIF数值汇总_R.xlsx")
model_input <- plain_path("PROJECT", "汇总/9个术前候选变量_临床模型比较--内部验证集/9个术前候选变量_训练后模型对象_R.rds")
training_oof_input <- plain_path("PROJECT", "汇总/9个术前候选变量_训练集OOF模型比较。训练集/9个术前候选变量_训练集重复交叉验证OOF预测概率.csv")
validation_prediction_input <- plain_path("PROJECT", "汇总/9个术前候选变量_临床模型比较--内部验证集/9个术前候选变量_验证集预测概率.csv")
output_dir <- plain_path("PROJECT", "汇总/LightGBM_性能与可解释性")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
copy_to_temp <- function(path, name) {
  target <- file.path(tempdir(), name)
  stopifnot(file.copy(path, target, overwrite = TRUE))
  target
}
training_file <- copy_to_temp(training_input, "lightgbm_training.xlsx")
vif_file <- copy_to_temp(vif_input, "lightgbm_vif.xlsx")
training_oof_file <- copy_to_temp(training_oof_input, "lightgbm_training_oof.csv")
validation_prediction_file <- copy_to_temp(validation_prediction_input, "lightgbm_validation_predictions.csv")

# Figure contract: a quantitative-grid figure pair. Figure 1 establishes
# discrimination, calibration and net benefit; Figure 2 explains the locked model.
palette <- c(
  train = "#0072B2", validation = "#D55E00", signal = "#0072B2",
  accent = "#C64E63", neutral_dark = "#252525", neutral_mid = "#7A7A7A",
  neutral_light = "#D9D9D9", treat_all = "#8C8C8C"
)
# R's optional Cairo module is unavailable on this workstation because its X11 dependency is absent.
# SVG remains the primary editable vector export and the base PDF device is used as a local fallback.
cairo_pdf_available <- FALSE
theme_publication <- function(base_size = 7.0) {
  theme_classic(base_size = base_size, base_family = "Helvetica") +
    theme(
      axis.line = element_line(linewidth = 0.35, colour = "#252525"),
      axis.ticks = element_line(linewidth = 0.35, colour = "#252525"),
      axis.title = element_text(size = base_size, colour = "#252525"),
      axis.text = element_text(size = base_size - 0.5, colour = "#252525"),
      plot.title = element_text(size = base_size + 0.5, face = "bold", colour = "#252525"),
      plot.subtitle = element_text(size = base_size - 0.3, colour = "#555555"),
      plot.caption = element_text(size = base_size - 1, colour = "#666666", hjust = 0),
      legend.title = element_text(size = base_size - 0.3),
      legend.text = element_text(size = base_size - 0.7),
      panel.grid = element_blank(),
      plot.tag = element_text(size = base_size + 2, face = "bold", colour = "#111111")
    )
}

save_figure_bundle <- function(plot, stem, width_mm, height_mm) {
  temporary_base <- file.path(tempdir(), paste0("lightgbm_fig_", stem))
  width_in <- width_mm / 25.4
  height_in <- height_mm / 25.4
  temporary_svg <- paste0(temporary_base, ".svg")
  temporary_pdf <- paste0(temporary_base, ".pdf")
  temporary_png <- paste0(temporary_base, "_600dpi.png")
  temporary_tiff <- paste0(temporary_base, "_600dpi.tiff")
  svglite::svglite(temporary_svg, width = width_in, height = height_in, bg = "white")
  print(plot)
  grDevices::dev.off()
  # Prefer cairo_pdf when supported; this local R installation may lack the optional X11 Cairo module.
  if (cairo_pdf_available && capabilities("cairo")) {
    grDevices::cairo_pdf(temporary_pdf, width = width_in, height = height_in, family = "Helvetica")
  } else {
    grDevices::pdf(temporary_pdf, width = width_in, height = height_in, family = "Helvetica", useDingbats = FALSE)
  }
  print(plot)
  grDevices::dev.off()
  if (!file.exists(temporary_pdf) || file.info(temporary_pdf)$size == 0) {
    grDevices::pdf(temporary_pdf, width = width_in, height = height_in, family = "Helvetica", useDingbats = FALSE)
    print(plot)
    grDevices::dev.off()
  }
  ragg::agg_png(temporary_png, width = width_in, height = height_in, units = "in", res = 600, background = "white")
  print(plot)
  grDevices::dev.off()
  ragg::agg_tiff(temporary_tiff, width = width_in, height = height_in, units = "in", res = 600, background = "white")
  print(plot)
  grDevices::dev.off()
  outputs <- c(temporary_svg, temporary_pdf, temporary_png, temporary_tiff)
  destination <- file.path(output_dir, c(
    paste0(stem, ".svg"), paste0(stem, ".pdf"), paste0(stem, "_600dpi.png"), paste0(stem, "_600dpi.tiff")
  ))
  copied <- mapply(function(source, target) file.copy(source, target, overwrite = TRUE), outputs, destination)
  if (!all(copied)) stop(paste0("Figure export copy failed for: ", paste(outputs[!copied], collapse = "; ")))
  destination
}

priority_candidates <- as.data.frame(readxl::read_excel(vif_file, sheet = "Priority_candidates"))
priority_candidates <- priority_candidates[priority_candidates$p_value < 0.05, , drop = FALSE]
if (nrow(priority_candidates) != 9) stop("Expected nine locked priority variables")
selected_variables <- priority_candidates$variable
continuous_tag <- "\u8fde\u7eed"
continuous_variables <- priority_candidates$variable[priority_candidates$variable_type == continuous_tag]
categorical_variables <- setdiff(selected_variables, continuous_variables)
if (length(continuous_variables) != 7 || length(categorical_variables) != 2) stop("Unexpected variable types")

# The VIF workbook fixes this candidate order; using its order avoids locale-dependent Chinese-name matching.
display_labels <- stats::setNames(c(
  "Prothrombin time (s)", "T stage", "Differentiation", "LDH (U/L)",
  "Thrombin time (s)", "Basophils (%)", "Age (years)", "HDL-C (mmol/L)", "Albumin/globulin ratio"
), selected_variables)

training_clinical <- as.data.frame(readxl::read_excel(training_file, sheet = 1))
outcome_name <- "\u543b\u5408\u53e3"
if (!all(c(outcome_name, selected_variables) %in% names(training_clinical))) stop("Training input lacks required columns")

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
categorical_columns <- names(train_data)[-1][vapply(train_data[-1], is.factor, logical(1))]
categorical_level_map <- lapply(categorical_columns, function(column) sort(unique(c(as.character(train_data[[column]]), "Missing"))))
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
    for (level in categorical_level_map[[column]]) output[[paste0(column, "_", make.names(level))]] <- as.numeric(values == level)
  }
  as.data.frame(output, check.names = FALSE)
}

impute_and_scale_training <- function(x) {
  medians <- vapply(x, function(column) {
    value <- stats::median(column, na.rm = TRUE)
    if (is.finite(value)) value else 0
  }, numeric(1))
  imputed <- as.data.frame(Map(function(column, median_value) {
    column[is.na(column)] <- median_value
    column
  }, x, medians), check.names = FALSE)
  names(imputed) <- names(medians)
  means <- vapply(imputed, mean, numeric(1))
  sds <- vapply(imputed, stats::sd, numeric(1))
  sds[!is.finite(sds) | sds == 0] <- 1
  list(values = sweep(sweep(as.matrix(imputed), 2, means, "-"), 2, sds, "/"), medians = medians, means = means, sds = sds)
}

train_predictors <- make_numeric_predictors(train_data)
nonconstant_columns <- vapply(train_predictors, function(column) length(unique(column[!is.na(column)])) > 1, logical(1))
train_predictors <- train_predictors[, nonconstant_columns, drop = FALSE]
prepared_training <- impute_and_scale_training(train_predictors)

fits <- readRDS(model_input)
if (!"LightGBM" %in% names(fits) || !inherits(fits$LightGBM, "lgb.Booster")) stop("Locked LightGBM model object is unavailable")
lightgbm_model <- fits$LightGBM
contributions <- as.matrix(predict(lightgbm_model, prepared_training$values, type = "contrib"))
if (nrow(contributions) != nrow(train_data) || ncol(contributions) != ncol(prepared_training$values) + 1) stop("Unexpected LightGBM SHAP contribution dimensions")
colnames(contributions) <- c(colnames(prepared_training$values), "Bias")
if (max(abs(rowSums(contributions) - qlogis(pmin(pmax(as.numeric(predict(lightgbm_model, prepared_training$values)), 1e-8), 1 - 1e-8)))) > 1e-5) stop("SHAP additivity check failed")

feature_mapping <- do.call(rbind, lapply(seq_along(selected_variables), function(index) {
  variable <- selected_variables[index]
  model_column <- paste0("x", index)
  encoded_columns <- if (variable %in% continuous_variables) model_column else paste0(model_column, "_", make.names(categorical_level_map[[model_column]]))
  encoded_columns <- intersect(encoded_columns, colnames(prepared_training$values))
  if (length(encoded_columns) == 0) stop(paste0("No active model columns found for ", variable))
  data.frame(Model_column = encoded_columns, Variable = variable, Display = unname(display_labels[variable]), stringsAsFactors = FALSE)
}))
if (!all(colnames(prepared_training$values) %in% feature_mapping$Model_column)) stop("SHAP feature mapping does not cover model columns")

shap_by_variable <- sapply(selected_variables, function(variable) {
  columns <- feature_mapping$Model_column[feature_mapping$Variable == variable]
  rowSums(contributions[, columns, drop = FALSE])
})
colnames(shap_by_variable) <- selected_variables
importance <- data.frame(
  Variable = selected_variables,
  Display = unname(display_labels[selected_variables]),
  Mean_abs_SHAP = colMeans(abs(shap_by_variable)),
  stringsAsFactors = FALSE
)
importance <- importance[order(importance$Mean_abs_SHAP, decreasing = TRUE), , drop = FALSE]
importance$Display <- factor(importance$Display, levels = rev(importance$Display))

feature_values <- lapply(seq_along(selected_variables), function(index) {
  variable <- selected_variables[index]
  values <- if (variable %in% continuous_variables) {
    as.numeric(train_data[[paste0("x", index)]])
  } else {
    category_values <- as.character(train_data[[paste0("x", index)]])
    factor(category_values, levels = categorical_level_map[[paste0("x", index)]], ordered = TRUE)
  }
  data.frame(Variable = variable, Display = unname(display_labels[variable]), Value = as.numeric(values), stringsAsFactors = FALSE)
})
feature_values <- do.call(rbind, feature_values)
feature_values$Patient <- rep(seq_len(nrow(train_data)), times = length(selected_variables))
shap_long <- do.call(rbind, lapply(seq_along(selected_variables), function(index) {
  variable <- selected_variables[index]
  data.frame(Patient = seq_len(nrow(train_data)), Variable = variable, Display = unname(display_labels[variable]), SHAP = shap_by_variable[, variable], stringsAsFactors = FALSE)
}))
shap_long <- merge(shap_long, feature_values, by = c("Patient", "Variable", "Display"), sort = FALSE)
shap_long$Feature_percentile <- ave(shap_long$Value, shap_long$Variable, FUN = function(value) {
  output <- rep(NA_real_, length(value))
  observed <- !is.na(value)
  if (sum(observed) > 1) output[observed] <- rank(value[observed], ties.method = "average") / sum(observed)
  output
})
shap_long$Display <- factor(shap_long$Display, levels = rev(as.character(importance$Display)))
# Deterministic vertical spreading avoids random sampling while keeping dense SHAP points visible.
shap_long$Bee_y <- as.numeric(shap_long$Display) + 0.26 * sin(shap_long$Patient * 2.399 + as.numeric(shap_long$Display) * 0.731)

training_oof <- utils::read.csv(training_oof_file, fileEncoding = "UTF-8")
validation_predictions <- utils::read.csv(validation_prediction_file, fileEncoding = "UTF-8")
training_lightgbm <- training_oof[training_oof$Model == "LightGBM", c("RowIndex", "Outcome", "Probability_Leak"), drop = FALSE]
training_lightgbm <- training_lightgbm[order(training_lightgbm$RowIndex), , drop = FALSE]
validation_lightgbm <- validation_predictions[validation_predictions$Model == "LightGBM", c("Outcome", "Probability_Leak"), drop = FALSE]
if (nrow(training_lightgbm) != 574 || sum(training_lightgbm$Outcome == "Leak") != 44 || nrow(validation_lightgbm) != 246 || sum(validation_lightgbm$Outcome == "Leak") != 19) stop("Unexpected performance cohort counts")

make_roc_data <- function(data, cohort) {
  outcome <- factor(data$Outcome, levels = c("NoLeak", "Leak"))
  roc_object <- pROC::roc(outcome, data$Probability_Leak, levels = c("NoLeak", "Leak"), direction = "<", quiet = TRUE)
  ci <- as.numeric(pROC::ci.auc(roc_object, method = "delong"))
  list(
    plot = data.frame(Cohort = cohort, False_positive_rate = 1 - roc_object$specificities, Sensitivity = roc_object$sensitivities),
    summary = data.frame(Cohort = cohort, N = nrow(data), Events = sum(outcome == "Leak"), AUC = as.numeric(pROC::auc(roc_object)), AUC_CI_low = ci[1], AUC_CI_high = ci[3], stringsAsFactors = FALSE)
  )
}
train_roc <- make_roc_data(training_lightgbm, "Training OOF")
validation_roc <- make_roc_data(validation_lightgbm, "Internal validation")
roc_data <- rbind(train_roc$plot, validation_roc$plot)
roc_summary <- rbind(train_roc$summary, validation_roc$summary)
roc_summary$Legend <- sprintf("%s  AUC %.3f (%.3f-%.3f)", roc_summary$Cohort, roc_summary$AUC, roc_summary$AUC_CI_low, roc_summary$AUC_CI_high)
roc_data$Legend <- roc_summary$Legend[match(roc_data$Cohort, roc_summary$Cohort)]

calibration_statistics <- function(data, cohort) {
  outcome <- as.integer(data$Outcome == "Leak")
  probability <- pmin(pmax(data$Probability_Leak, 1e-6), 1 - 1e-6)
  logit_probability <- qlogis(probability)
  intercept_fit <- stats::glm(outcome ~ offset(logit_probability), family = stats::binomial())
  slope_fit <- stats::glm(outcome ~ logit_probability, family = stats::binomial())
  data.frame(
    Cohort = cohort, N = length(outcome), Events = sum(outcome), Brier_score = mean((probability - outcome)^2),
    Calibration_intercept = unname(stats::coef(intercept_fit)[1]), Calibration_slope = unname(stats::coef(slope_fit)[2]), stringsAsFactors = FALSE
  )
}
calibration_bins <- function(data, cohort, groups = 5) {
  outcome <- as.integer(data$Outcome == "Leak")
  probability <- data$Probability_Leak
  group <- pmin(groups, ceiling(rank(probability, ties.method = "first") / length(probability) * groups))
  split_index <- split(seq_along(probability), group)
  result <- do.call(rbind, lapply(split_index, function(index) {
    observed <- mean(outcome[index])
    se <- sqrt(observed * (1 - observed) / length(index))
    data.frame(Cohort = cohort, Bin = paste0("Q", unique(group[index])), N = length(index), Predicted = mean(probability[index]), Observed = observed, Lower = max(0, observed - 1.96 * se), Upper = min(1, observed + 1.96 * se))
  }))
  result[order(result$Predicted), , drop = FALSE]
}
calibration_data <- rbind(calibration_bins(training_lightgbm, "Training OOF"), calibration_bins(validation_lightgbm, "Internal validation"))
calibration_summary <- rbind(calibration_statistics(training_lightgbm, "Training OOF"), calibration_statistics(validation_lightgbm, "Internal validation"))

dca_curve <- function(data, cohort, thresholds) {
  outcome <- as.integer(data$Outcome == "Leak")
  prevalence <- mean(outcome)
  do.call(rbind, lapply(thresholds, function(threshold) {
    predicted_positive <- data$Probability_Leak >= threshold
    tp <- sum(predicted_positive & outcome == 1)
    fp <- sum(predicted_positive & outcome == 0)
    weight <- threshold / (1 - threshold)
    data.frame(
      Cohort = cohort, Threshold = threshold,
      Model = tp / length(outcome) - fp / length(outcome) * weight,
      Treat_all = prevalence - (1 - prevalence) * weight,
      Treat_none = 0
    )
  }))
}
dca_data <- rbind(dca_curve(training_lightgbm, "Training OOF", seq(0.01, 0.30, by = 0.01)), dca_curve(validation_lightgbm, "Internal validation", seq(0.01, 0.30, by = 0.01)))
dca_long <- rbind(
  data.frame(Cohort = dca_data$Cohort, Threshold = dca_data$Threshold, Strategy = "LightGBM", Net_benefit = dca_data$Model),
  data.frame(Cohort = dca_data$Cohort, Threshold = dca_data$Threshold, Strategy = "Treat all", Net_benefit = dca_data$Treat_all),
  data.frame(Cohort = dca_data$Cohort, Threshold = dca_data$Threshold, Strategy = "Treat none", Net_benefit = dca_data$Treat_none)
)

p_roc <- ggplot(roc_data, aes(False_positive_rate, Sensitivity, colour = Cohort)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", linewidth = 0.35, colour = palette["neutral_mid"]) +
  geom_step(linewidth = 0.85, direction = "vh") +
  annotate("text", x = 0.95, y = 0.36, label = sprintf("OOF AUC: %.3f", train_roc$summary$AUC), colour = palette[["train"]], hjust = 1, size = 2.15) +
  annotate("text", x = 0.95, y = 0.25, label = sprintf("Validation AUC: %.3f", validation_roc$summary$AUC), colour = palette[["validation"]], hjust = 1, size = 2.15) +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
  scale_colour_manual(values = c("Training OOF" = palette[["train"]], "Internal validation" = palette[["validation"]])) +
  labs(title = "Discrimination", x = "1 - Specificity", y = "Sensitivity", colour = NULL) +
  theme_publication() + theme(legend.position = "none")

p_calibration <- ggplot(calibration_data, aes(Predicted, Observed, colour = Cohort, group = Cohort)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", linewidth = 0.35, colour = palette["neutral_mid"]) +
  geom_errorbar(aes(ymin = Lower, ymax = Upper), width = 0.006, linewidth = 0.35) +
  geom_line(linewidth = 0.55) + geom_point(size = 1.7) +
  coord_equal(xlim = c(0, 0.35), ylim = c(0, 0.35), expand = FALSE) +
  scale_colour_manual(values = c("Training OOF" = palette[["train"]], "Internal validation" = palette[["validation"]])) +
  labs(title = "Calibration", subtitle = "Five equal-frequency risk groups", x = "Mean predicted risk", y = "Observed leakage proportion", colour = NULL) +
  theme_publication() + theme(legend.position = "bottom")

dca_models <- dca_long[dca_long$Strategy == "LightGBM", , drop = FALSE]
dca_treat_all <- dca_long[dca_long$Strategy == "Treat all" & dca_long$Cohort == "Internal validation", , drop = FALSE]
p_dca <- ggplot() +
  geom_hline(yintercept = 0, linewidth = 0.35, colour = palette[["neutral_mid"]]) +
  geom_line(data = dca_models, aes(Threshold, Net_benefit, colour = Cohort, linetype = Cohort), linewidth = 0.75) +
  geom_line(data = dca_treat_all, aes(Threshold, Net_benefit), linewidth = 0.55, linetype = "dotted", colour = palette[["treat_all"]]) +
  annotate("text", x = 0.25, y = 0.017, label = "Treat all", size = 2.0, colour = palette[["treat_all"]], hjust = 0) +
  annotate("text", x = 0.25, y = 0.004, label = "Treat none", size = 2.0, colour = palette[["neutral_mid"]], hjust = 0) +
  scale_colour_manual(values = c("Training OOF" = palette[["train"]], "Internal validation" = palette[["validation"]])) +
  scale_linetype_manual(values = c("Training OOF" = "solid", "Internal validation" = "dashed")) +
  coord_cartesian(xlim = c(0.01, 0.30), ylim = c(-0.03, 0.09)) +
  labs(title = "Clinical utility", subtitle = "Decision-curve analysis", x = "Threshold probability", y = "Net benefit", colour = NULL, linetype = NULL) +
  theme_publication() + theme(legend.position = "bottom", legend.text = element_text(size = 5.8))

figure_1 <- (p_roc | p_calibration | p_dca) +
  patchwork::plot_annotation(tag_levels = "A", title = "Clinical LightGBM model performance") &
  theme(plot.tag = element_text(size = 10, face = "bold"))

p_importance <- ggplot(importance, aes(x = Mean_abs_SHAP, y = Display, fill = Mean_abs_SHAP)) +
  geom_col(width = 0.72) +
  geom_text(aes(label = sprintf("%.3f", Mean_abs_SHAP)), hjust = -0.15, size = 2.4, colour = palette["neutral_dark"]) +
  scale_fill_gradient(low = "#9CC6DD", high = palette["accent"], guide = "none") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.16))) +
  labs(title = "Global feature importance", x = "Mean |SHAP value| (log-odds)", y = NULL, tag = "A") +
  theme_publication() + theme(axis.text.y = element_text(size = 6.2))

p_beeswarm <- ggplot(shap_long, aes(x = SHAP, y = Bee_y, colour = Feature_percentile)) +
  geom_vline(xintercept = 0, linewidth = 0.35, colour = palette["neutral_mid"]) +
  geom_point(size = 0.62, alpha = 0.78, na.rm = TRUE) +
  scale_colour_gradientn(colours = c("#2166AC", "#F7F7F7", "#D73027"), limits = c(0, 1), oob = scales::squish, name = "Feature value", breaks = c(0, 1), labels = c("Low", "High")) +
  scale_y_continuous(breaks = seq_along(levels(shap_long$Display)), labels = levels(shap_long$Display)) +
  labs(title = "SHAP summary", x = "SHAP value (impact on model output)", y = NULL, tag = "B") +
  theme_publication() + theme(axis.text.y = element_text(size = 6.2), legend.position = "right")

dependence_plots <- lapply(seq_along(selected_variables), function(index) {
  variable <- selected_variables[index]
  label <- unname(display_labels[variable])
  data <- shap_long[shap_long$Variable == variable, , drop = FALSE]
  categorical <- variable %in% categorical_variables
  plot <- ggplot(data, aes(Value, SHAP)) +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.3, colour = palette["neutral_mid"]) +
    geom_point(colour = palette["signal"], size = 0.72, alpha = 0.78, na.rm = TRUE) +
    labs(title = paste0(letters[index], "  ", label), x = label, y = "SHAP value") +
    theme_publication(6.0) +
    theme(plot.title = element_text(size = 6.4, face = "bold"), axis.title = element_text(size = 5.7), axis.text = element_text(size = 5.5))
  if (categorical) plot <- plot + scale_x_continuous(breaks = sort(unique(data$Value[!is.na(data$Value)]))) + geom_jitter(width = 0.08, height = 0, colour = palette["signal"], size = 0.5, alpha = 0.42, na.rm = TRUE)
  if (index == 1) plot <- plot + labs(tag = "C")
  plot
})
p_dependence <- patchwork::wrap_plots(dependence_plots, ncol = 3) + patchwork::plot_annotation(title = "Feature dependence")

figure_2 <- (p_importance | p_beeswarm) / p_dependence +
  patchwork::plot_layout(heights = c(1.0, 2.25)) +
  patchwork::plot_annotation(title = "LightGBM model interpretation") &
  theme(plot.tag = element_text(size = 10, face = "bold"))

figure_1_paths <- save_figure_bundle(figure_1, "Figure_1_LightGBM_performance", width_mm = 183, height_mm = 70)
figure_2_paths <- save_figure_bundle(figure_2, "Figure_2_LightGBM_SHAP_interpretation", width_mm = 183, height_mm = 235)

calibration_table <- merge(calibration_summary, data.frame(Cohort = roc_summary$Cohort, AUC = roc_summary$AUC, AUC_CI_low = roc_summary$AUC_CI_low, AUC_CI_high = roc_summary$AUC_CI_high), by = "Cohort", all.x = TRUE)
utils::write.csv(importance, file.path(output_dir, "LightGBM_SHAP_global_importance_source_data.csv"), row.names = FALSE, fileEncoding = "UTF-8")
utils::write.csv(shap_long, file.path(output_dir, "LightGBM_SHAP_patient_level_source_data.csv"), row.names = FALSE, fileEncoding = "UTF-8")
utils::write.csv(calibration_data, file.path(output_dir, "LightGBM_calibration_curve_source_data.csv"), row.names = FALSE, fileEncoding = "UTF-8")
utils::write.csv(calibration_table, file.path(output_dir, "LightGBM_performance_and_calibration_summary.csv"), row.names = FALSE, fileEncoding = "UTF-8")
utils::write.csv(dca_long, file.path(output_dir, "LightGBM_DCA_source_data.csv"), row.names = FALSE, fileEncoding = "UTF-8")

workbook <- openxlsx::createWorkbook()
header_style <- openxlsx::createStyle(fgFill = "#24415B", fontColour = "#FFFFFF", textDecoration = "bold", halign = "center", valign = "center", wrapText = TRUE)
number_style <- openxlsx::createStyle(numFmt = "0.000")
add_sheet <- function(name, data) {
  openxlsx::addWorksheet(workbook, name)
  openxlsx::writeData(workbook, name, data, withFilter = TRUE)
  openxlsx::addStyle(workbook, name, header_style, rows = 1, cols = seq_len(ncol(data)), gridExpand = TRUE)
  openxlsx::freezePane(workbook, name, firstRow = TRUE)
  openxlsx::setColWidths(workbook, name, cols = seq_len(ncol(data)), widths = "auto")
}
add_sheet("Performance_calibration", calibration_table)
openxlsx::addStyle(workbook, "Performance_calibration", number_style, rows = 2:(nrow(calibration_table) + 1), cols = 4:ncol(calibration_table), gridExpand = TRUE, stack = TRUE)
add_sheet("Calibration_curve", calibration_data)
openxlsx::addStyle(workbook, "Calibration_curve", number_style, rows = 2:(nrow(calibration_data) + 1), cols = 4:ncol(calibration_data), gridExpand = TRUE, stack = TRUE)
add_sheet("DCA_curve", dca_long)
openxlsx::addStyle(workbook, "DCA_curve", number_style, rows = 2:(nrow(dca_long) + 1), cols = 3:ncol(dca_long), gridExpand = TRUE, stack = TRUE)
add_sheet("SHAP_importance", importance)
openxlsx::addStyle(workbook, "SHAP_importance", number_style, rows = 2:(nrow(importance) + 1), cols = 3, gridExpand = TRUE, stack = TRUE)
add_sheet("SHAP_patient_values", shap_long)
openxlsx::addStyle(workbook, "SHAP_patient_values", number_style, rows = 2:(nrow(shap_long) + 1), cols = c(4:6), gridExpand = TRUE, stack = TRUE)
notes <- data.frame(
  Item = c("Core conclusion", "Figure 1 data", "Figure 2 data", "Training set", "Internal validation", "SHAP method", "Feature handling", "Calibration", "DCA", "Model locking", "Exclusions", "Interpretation limit"),
  Detail = c(
    "The locked nine-variable clinical LightGBM model shows discrimination in training OOF and internal validation; SHAP identifies the variables driving its predictions.",
    "Training repeated-CV patient-level OOF predictions and sealed internal validation predictions for LightGBM only.",
    "Exact SHAP contributions from LightGBM native predcontrib on the final model fit to all training patients.",
    "n=574; leakage events=44. Training performance is repeated 5-fold CV with 5 repeats and patient-level average OOF probabilities.",
    "n=246; leakage events=19. Validation is not used for model selection, tuning, threshold selection, or SHAP ranking.",
    "Native TreeSHAP contributions (log-odds scale); additivity checked against final-model prediction log-odds.",
    "The model retains 9 investigator-confirmed preoperative variables. Continuous exact zero values are handled as missing before model imputation; categorical variables are one-hot encoded and SHAP values are summed back to the original variable.",
    "Five equal-frequency risk groups; vertical bars are approximate 95% binomial intervals. Reported Brier score, calibration intercept and slope are descriptive and imprecise in the small-event validation cohort.",
    "Net benefit computed for thresholds 0.01-0.30; treat-all and treat-none are reference strategies.",
    "LightGBM architecture and hyperparameters were locked before this figure generation. No model was retrained or re-tuned.",
    "No patients, variables or outcomes were excluded from the training or validation performance plots. SHAP uses all 574 training patients.",
    "Internal validation has 19 events; confidence intervals are wide. Figures support internal validation and explanation, not external generalizability or causal inference."
  ),
  stringsAsFactors = FALSE
)
add_sheet("Figure_notes", notes)
openxlsx::saveWorkbook(workbook, file.path(output_dir, "LightGBM_性能与可解释性_统计与源数据.xlsx"), overwrite = TRUE)

writeLines(c(
  paste0("Figure_1_outputs=", paste(figure_1_paths, collapse = "; ")),
  paste0("Figure_2_outputs=", paste(figure_2_paths, collapse = "; ")),
  "figure_1=ROC, five-bin calibration, DCA; Figure_2=mean absolute SHAP, SHAP beeswarm, nine dependence plots",
  "training=574 patients, 44 leakage events; training performance uses averaged repeated-CV OOF probabilities",
  "internal_validation=246 patients, 19 leakage events; no validation-set selection or tuning",
  "SHAP=LightGBM native predcontrib on final locked full-training model; additivity checked",
  "exports=SVG, PDF, 600 dpi PNG, 600 dpi TIFF"
), file.path(output_dir, "LightGBM_性能与可解释性_运行说明.txt"))
writeLines(c(
  "Figure contract: quantitative-grid figure pair; Figure 1 is performance validation and Figure 2 is model explanation.",
  "QA: training OOF n=574/events=44 and internal validation n=246/events=19 were verified.",
  "QA: all nine locked clinical variables and all 574 training patients were retained for SHAP.",
  "QA: SHAP feature columns were mapped back to the original nine variables; one-hot encoded categorical SHAP values were summed.",
  "QA: native SHAP additivity was checked against LightGBM prediction log-odds.",
  "QA: calibration uses five equal-frequency risk groups; DCA thresholds span 0.01 to 0.30.",
  "QA: no validation outcome was used for model selection, hyperparameter tuning, threshold setting or SHAP ranking.",
  "QA: all figures are white-background vector/raster exports with editable SVG/PDF text; SVG is the primary editable vector file because this R installation lacks optional Cairo/X11 support."
), file.path(output_dir, "LightGBM_性能与可解释性_QA说明.txt"))

print(importance[, c("Display", "Mean_abs_SHAP")])
print(calibration_table)
