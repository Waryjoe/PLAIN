#!/usr/bin/env Rscript
source("scripts/plain_paths.R")

# Training-only selection of the locked stomach + tumor CT2 radiomics matrix.
# Validation data are read only after the Rad-score is locked.

options(stringsAsFactors = FALSE)
user_lib <- plain_path("R_LIBRARY", "")
if (dir.exists(user_lib)) .libPaths(c(user_lib, .libPaths()))
suppressPackageStartupMessages({
  library(readxl)
  library(openxlsx)
  library(glmnet)
  library(pROC)
  library(ggplot2)
  library(patchwork)
})

seed <- 20260806L
set.seed(seed)

arguments <- commandArgs(trailingOnly = TRUE)
read_argument <- function(name, default = NULL) {
  value <- sub(paste0("^", name, "="), "", arguments[startsWith(arguments, paste0(name, "="))])
  if (length(value)) value[[1]] else default
}

project_root <- plain_path("PROJECT", "")
default_training_input <- file.path(
  project_root,
  "汇总/影像组学_训练验证集_当前提取版本_待一致化/影像组学训练集_70%_574例_3190特征.xlsx"
)
default_validation_input <- file.path(
  project_root,
  "汇总/影像组学_训练验证集_当前提取版本_待一致化/影像组学验证集_30%_246例_3190特征.xlsx"
)
training_input <- read_argument("--training-input", default_training_input)
validation_input <- read_argument("--validation-input", default_validation_input)
allowed_features_path <- read_argument("--allowed-features", NULL)
output_dir <- read_argument(
  "--output-dir",
  file.path(project_root, "汇总/胃肿瘤双ROI影像组学_训练集LASSO稳定性筛选")
)
code_archive_dir <- file.path(project_root, "代码与复现说明/scripts")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(code_archive_dir, recursive = TRUE, showWarnings = FALSE)

# Prespecified choices for this exploratory-but-reproducible selection analysis.
univariable_screen_p <- 0.05
minimum_screened_features <- 20L
spearman_redundancy_cutoff <- 0.90
resample_count <- 100L
resample_fraction <- 0.80
selection_frequency_cutoff <- 0.70
direction_consistency_cutoff <- 0.90
outer_folds <- as.integer(read_argument("--folds", "5"))
if (!is.finite(outer_folds) || outer_folds < 2L) stop("--folds must be an integer >= 2")
roi_mode <- tolower(read_argument("--roi", "combined"))
if (!roi_mode %in% c("combined", "stomach", "tumor")) stop("--roi must be combined, stomach, or tumor")
analysis_label <- switch(
  roi_mode,
  combined = "Stomach + tumor dual-ROI",
  stomach = "Stomach ROI",
  tumor = "Tumor ROI"
)
file_tag <- switch(roi_mode, combined = "双ROI", stomach = "胃ROI", tumor = "肿瘤ROI")

copy_to_temp <- function(path, name) {
  destination <- file.path(tempdir(), name)
  stopifnot(file.copy(path, destination, overwrite = TRUE))
  destination
}

read_feature_table <- function(path, expected_rows) {
  temporary_path <- copy_to_temp(path, paste0("input_", expected_rows, ".xlsx"))
  data <- as.data.frame(readxl::read_excel(temporary_path, sheet = 1))
  if (nrow(data) != expected_rows || ncol(data) != 3194) {
    stop("Unexpected patient-level radiomics matrix dimensions")
  }
  feature_names <- names(data)[grepl("^(stomach|tumor)__", names(data))]
  if (length(feature_names) != 3190L ||
      sum(startsWith(feature_names, "stomach__")) != 1595L ||
      sum(startsWith(feature_names, "tumor__")) != 1595L) {
    stop("Expected 1595 stomach and 1595 tumor feature columns")
  }
  list(
    ids = as.character(data[[1]]),
    case_ids = as.character(data[[2]]),
    outcome = as.integer(data[[3]]),
    split = as.character(data[[4]]),
    feature_names = feature_names,
    x = as.matrix(data[, feature_names, drop = FALSE])
  )
}

make_stratified_foldid <- function(y, k, local_seed) {
  set.seed(local_seed)
  foldid <- integer(length(y))
  for (class_value in sort(unique(y))) {
    indices <- which(y == class_value)
    foldid[indices] <- sample(rep(seq_len(k), length.out = length(indices)))
  }
  foldid
}

score_test_pvalues <- function(x, y) {
  prevalence <- mean(y)
  centered <- sweep(x, 2, colMeans(x), "-")
  information <- prevalence * (1 - prevalence) * colSums(centered^2)
  score <- as.numeric(crossprod(y - prevalence, centered))
  z_score <- score / sqrt(information)
  p_value <- 2 * stats::pnorm(abs(z_score), lower.tail = FALSE)
  p_value[!is.finite(p_value)] <- 1
  p_value
}

univariable_results <- function(x, y, feature_names) {
  n_positive <- sum(y == 1)
  n_negative <- sum(y == 0)
  output <- vector("list", ncol(x))
  for (index in seq_len(ncol(x))) {
    value <- x[, index]
    standard_deviation <- stats::sd(value)
    if (!is.finite(standard_deviation) || standard_deviation == 0) next
    z_value <- (value - mean(value)) / standard_deviation
    fit <- try(suppressWarnings(stats::glm.fit(
      x = cbind(Intercept = 1, Feature = z_value),
      y = y,
      family = stats::binomial()
    )), silent = TRUE)
    if (inherits(fit, "try-error") || any(!is.finite(fit$coefficients))) next
    if (fit$rank < 2L) next
    # glm.fit returns a low-level list rather than a glm object. Recover the
    # unscaled coefficient covariance directly from its QR decomposition.
    covariance <- try(chol2inv(qr.R(fit$qr)[seq_len(fit$rank), seq_len(fit$rank), drop = FALSE]), silent = TRUE)
    if (inherits(covariance, "try-error") || !is.finite(covariance[2, 2])) next
    beta <- fit$coefficients[2]
    standard_error <- sqrt(covariance[2, 2])
    z_statistic <- beta / standard_error
    p_value <- 2 * stats::pnorm(abs(z_statistic), lower.tail = FALSE)
    ranks <- rank(value, ties.method = "average")
    auc <- (sum(ranks[y == 1]) - n_positive * (n_positive + 1) / 2) /
      (n_positive * n_negative)
    output[[index]] <- data.frame(
      Feature = feature_names[index],
      ROI = ifelse(startsWith(feature_names[index], "stomach__"), "Stomach", "Tumor"),
      Mean_train = mean(value),
      SD_train = standard_deviation,
      OR_per_1SD = exp(beta),
      OR_95CI_low = exp(beta - 1.96 * standard_error),
      OR_95CI_high = exp(beta + 1.96 * standard_error),
      Beta_per_1SD = beta,
      P_value = p_value,
      Single_feature_AUC = auc,
      AUC_direction_free = max(auc, 1 - auc),
      stringsAsFactors = FALSE
    )
  }
  result <- do.call(rbind, output[!vapply(output, is.null, logical(1))])
  result$BH_FDR <- stats::p.adjust(result$P_value, method = "BH")
  result <- result[order(result$P_value, result$Feature), , drop = FALSE]
  rownames(result) <- NULL
  result
}

screen_feature_indices <- function(p_values) {
  indices <- which(p_values < univariable_screen_p)
  fallback_used <- FALSE
  if (length(indices) < minimum_screened_features) {
    indices <- order(p_values, decreasing = FALSE)[seq_len(min(minimum_screened_features, length(p_values)))]
    fallback_used <- TRUE
  }
  list(indices = indices, fallback_used = fallback_used)
}

prune_spearman <- function(x, feature_names, p_values) {
  if (ncol(x) == 1) {
    return(list(keep = 1L, audit = data.frame(
      Feature = feature_names, Decision = "Retained", Reason = "Only screened feature",
      Retained_against = NA_character_, Spearman_rho = NA_real_, P_value = p_values,
      stringsAsFactors = FALSE
    )))
  }
  correlation <- suppressWarnings(stats::cor(x, method = "spearman"))
  ordered <- order(p_values, feature_names)
  retained <- integer()
  audit <- vector("list", length(ordered))
  for (position in seq_along(ordered)) {
    index <- ordered[position]
    if (!length(retained)) {
      retained <- c(retained, index)
      audit[[position]] <- data.frame(
        Feature = feature_names[index], Decision = "Retained", Reason = "Best univariable P value",
        Retained_against = NA_character_, Spearman_rho = NA_real_, P_value = p_values[index],
        stringsAsFactors = FALSE
      )
      next
    }
    correlations <- abs(correlation[index, retained])
    strongest <- which.max(correlations)
    if (max(correlations) > spearman_redundancy_cutoff) {
      reference <- retained[strongest]
      audit[[position]] <- data.frame(
        Feature = feature_names[index], Decision = "Removed", Reason = paste0("|Spearman rho| > ", spearman_redundancy_cutoff),
        Retained_against = feature_names[reference], Spearman_rho = correlation[index, reference],
        P_value = p_values[index], stringsAsFactors = FALSE
      )
    } else {
      retained <- c(retained, index)
      audit[[position]] <- data.frame(
        Feature = feature_names[index], Decision = "Retained", Reason = "No redundant retained feature",
        Retained_against = NA_character_, Spearman_rho = NA_real_, P_value = p_values[index],
        stringsAsFactors = FALSE
      )
    }
  }
  list(keep = retained, audit = do.call(rbind, audit))
}

fit_lasso <- function(x, y, local_seed) {
  means <- colMeans(x)
  standard_deviations <- apply(x, 2, stats::sd)
  standard_deviations[!is.finite(standard_deviations) | standard_deviations == 0] <- 1
  scaled <- sweep(sweep(x, 2, means, "-"), 2, standard_deviations, "/")
  foldid <- make_stratified_foldid(y, outer_folds, local_seed)
  cv_fit <- glmnet::cv.glmnet(
    x = scaled, y = y, family = "binomial", alpha = 1,
    standardize = FALSE, foldid = foldid, type.measure = "deviance"
  )
  coefficient_matrix <- as.matrix(stats::coef(cv_fit, s = "lambda.1se"))
  selected <- rownames(coefficient_matrix)[coefficient_matrix[, 1] != 0]
  selected <- setdiff(selected, "(Intercept)")
  list(
    cv_fit = cv_fit, means = means, standard_deviations = standard_deviations,
    scaled = scaled, selected = selected, coefficients = coefficient_matrix
  )
}

calculate_auc <- function(y, probability) {
  roc_object <- pROC::roc(y, probability, levels = c(0, 1), direction = "<", quiet = TRUE)
  interval <- as.numeric(pROC::ci.auc(roc_object, method = "delong"))
  c(AUC = as.numeric(pROC::auc(roc_object)), AUC_95CI_low = interval[1], AUC_95CI_high = interval[3])
}

training <- read_feature_table(training_input, 574L)
validation <- read_feature_table(validation_input, 246L)
if (!identical(training$feature_names, validation$feature_names)) stop("Training/validation feature schema mismatch")
if (!identical(sort(unique(training$outcome)), c(0L, 1L)) || sum(training$outcome == 1) != 44L) stop("Unexpected training outcome")
if (sum(validation$outcome == 1) != 19L) stop("Unexpected validation outcome")

if (roi_mode != "combined") {
  roi_prefix <- paste0(roi_mode, "__")
  training$feature_names <- training$feature_names[startsWith(training$feature_names, roi_prefix)]
  validation$feature_names <- validation$feature_names[startsWith(validation$feature_names, roi_prefix)]
}
if (!is.null(allowed_features_path)) {
  if (!file.exists(allowed_features_path)) stop("--allowed-features does not exist")
  allowed_table <- utils::read.csv(allowed_features_path, check.names = FALSE, stringsAsFactors = FALSE)
  names(allowed_table)[grep("Feature$", names(allowed_table))] <- "Feature"
  if (!"Feature" %in% names(allowed_table)) stop("--allowed-features must contain a Feature column")
  allowed_features <- unique(as.character(allowed_table$Feature))
  training$feature_names <- training$feature_names[training$feature_names %in% allowed_features]
  validation$feature_names <- validation$feature_names[validation$feature_names %in% allowed_features]
  if (!length(training$feature_names)) stop("No requested analysis features were retained by --allowed-features")
}
if (length(training$feature_names) != length(validation$feature_names) ||
    !identical(training$feature_names, validation$feature_names)) {
  stop("Unexpected ROI feature pool width")
}

x_all <- apply(training$x[, training$feature_names, drop = FALSE], 2, as.numeric)
colnames(x_all) <- training$feature_names
x_validation_all <- apply(validation$x[, validation$feature_names, drop = FALSE], 2, as.numeric)
colnames(x_validation_all) <- validation$feature_names
if (any(!is.finite(x_all)) || any(!is.finite(x_validation_all))) stop("Non-finite radiomics values detected")

feature_sd <- apply(x_all, 2, stats::sd)
feature_iqr <- apply(x_all, 2, stats::IQR)
valid_indices <- which(is.finite(feature_sd) & feature_sd > 1e-10 & is.finite(feature_iqr) & feature_iqr > 1e-10)
x_qc <- x_all[, valid_indices, drop = FALSE]
feature_qc <- data.frame(
  Feature = colnames(x_all),
  ROI = ifelse(startsWith(colnames(x_all), "stomach__"), "Stomach", "Tumor"),
  SD_train = feature_sd,
  IQR_train = feature_iqr,
  Finite_values = colSums(is.finite(x_all)),
  Decision = ifelse(seq_len(ncol(x_all)) %in% valid_indices, "Retained", "Removed"),
  Reason = ifelse(seq_len(ncol(x_all)) %in% valid_indices, "Finite and non-constant", "Constant or near-zero variance"),
  stringsAsFactors = FALSE
)

univariable <- univariable_results(x_qc, training$outcome, colnames(x_qc))
univariable_lookup <- setNames(univariable$P_value, univariable$Feature)
screen <- screen_feature_indices(univariable$P_value)
screened_features <- univariable$Feature[screen$indices]
x_screened <- x_qc[, screened_features, drop = FALSE]
screened_p <- univariable_lookup[screened_features]
pruned <- prune_spearman(x_screened, screened_features, screened_p)
final_candidate_features <- screened_features[pruned$keep]
x_candidates <- x_qc[, final_candidate_features, drop = FALSE]

final_lasso <- fit_lasso(x_candidates, training$outcome, seed + 1000L)
final_lasso_table <- data.frame(
  Feature = rownames(final_lasso$coefficients),
  Coefficient_lambda_1se = final_lasso$coefficients[, 1],
  Selected_lambda_1se = final_lasso$coefficients[, 1] != 0,
  stringsAsFactors = FALSE
)
final_lasso_table <- merge(final_lasso_table, univariable[, c("Feature", "ROI", "P_value", "BH_FDR", "AUC_direction_free")], by = "Feature", all.x = TRUE, sort = FALSE)

# Repeats the selection cascade in stratified 80% training subsamples. The
# quick score test is used only for the preliminary screen; LASSO is fully
# re-tuned by the requested stratified CV within every resample.
selection_count <- setNames(integer(ncol(x_qc)), colnames(x_qc))
positive_count <- setNames(integer(ncol(x_qc)), colnames(x_qc))
negative_count <- setNames(integer(ncol(x_qc)), colnames(x_qc))
resample_log <- vector("list", resample_count)
set.seed(seed + 2000L)
for (resample_index in seq_len(resample_count)) {
  pos_indices <- which(training$outcome == 1)
  neg_indices <- which(training$outcome == 0)
  chosen <- c(
    sample(pos_indices, size = floor(length(pos_indices) * resample_fraction), replace = FALSE),
    sample(neg_indices, size = floor(length(neg_indices) * resample_fraction), replace = FALSE)
  )
  y_subsample <- training$outcome[chosen]
  x_subsample <- x_qc[chosen, , drop = FALSE]
  p_subsample <- score_test_pvalues(x_subsample, y_subsample)
  screen_subsample <- screen_feature_indices(p_subsample)
  selected_screened <- colnames(x_subsample)[screen_subsample$indices]
  pruning_subsample <- prune_spearman(
    x_subsample[, selected_screened, drop = FALSE], selected_screened,
    p_subsample[screen_subsample$indices]
  )
  candidate_subsample <- selected_screened[pruning_subsample$keep]
  lasso_subsample <- fit_lasso(
    x_subsample[, candidate_subsample, drop = FALSE], y_subsample,
    seed + 2000L + resample_index
  )
  selected <- lasso_subsample$selected
  if (length(selected)) {
    coefficient_values <- lasso_subsample$coefficients[selected, 1]
    selection_count[selected] <- selection_count[selected] + 1L
    positive_count[selected[coefficient_values > 0]] <- positive_count[selected[coefficient_values > 0]] + 1L
    negative_count[selected[coefficient_values < 0]] <- negative_count[selected[coefficient_values < 0]] + 1L
  }
  resample_log[[resample_index]] <- data.frame(
    Resample = resample_index,
    Patients = length(chosen), Events = sum(y_subsample),
    Screened_features = length(selected_screened),
    Correlation_pruned_features = length(candidate_subsample),
    LASSO_lambda_1se_features = length(selected),
    Lambda_1se = lasso_subsample$cv_fit$lambda.1se,
    stringsAsFactors = FALSE
  )
  if (resample_index %% 10L == 0L) message("Completed resample ", resample_index, " of ", resample_count)
}
resample_log <- do.call(rbind, resample_log)
direction <- ifelse(positive_count >= negative_count, "Positive", "Negative")
direction_consistency <- ifelse(selection_count > 0, pmax(positive_count, negative_count) / selection_count, NA_real_)
stability <- data.frame(
  Feature = names(selection_count),
  ROI = ifelse(startsWith(names(selection_count), "stomach__"), "Stomach", "Tumor"),
  Selection_count = as.integer(selection_count),
  Selection_frequency = as.numeric(selection_count) / resample_count,
  Positive_coefficient_count = as.integer(positive_count),
  Negative_coefficient_count = as.integer(negative_count),
  Predominant_direction = ifelse(selection_count > 0, direction, NA_character_),
  Direction_consistency = direction_consistency,
  Stable = selection_count / resample_count >= selection_frequency_cutoff & direction_consistency >= direction_consistency_cutoff,
  stringsAsFactors = FALSE
)
stability <- merge(stability, univariable[, c("Feature", "P_value", "BH_FDR", "AUC_direction_free")], by = "Feature", all.x = TRUE, sort = FALSE)
stability <- stability[order(-stability$Selection_frequency, -stability$Direction_consistency, stability$P_value), , drop = FALSE]

stable_features <- stability$Feature[stability$Stable]
locked_result <- NULL
if (length(stable_features) > 0) {
  locked_result <- fit_lasso(x_qc[, stable_features, drop = FALSE], training$outcome, seed + 5000L)
  locked_selected <- locked_result$selected
  if (!length(locked_selected)) {
    locked_result <- NULL
  }
}

final_radscore <- data.frame()
prediction_table <- data.frame(
  Hospital_ID = c(training$ids, validation$ids),
  Case_ID = c(training$case_ids, validation$case_ids),
  Dataset = c(rep("Training", length(training$outcome)), rep("Validation", length(validation$outcome))),
  Outcome = c(training$outcome, validation$outcome),
  Rad_score = NA_real_,
  Predicted_probability = NA_real_,
  stringsAsFactors = FALSE
)
performance_table <- data.frame()
if (!is.null(locked_result)) {
  locked_coefficients <- locked_result$coefficients
  locked_selected <- locked_result$selected
  final_radscore <- data.frame(
    Feature = locked_selected,
    ROI = ifelse(startsWith(locked_selected, "stomach__"), "Stomach", "Tumor"),
    Mean_training = locked_result$means[locked_selected],
    SD_training = locked_result$standard_deviations[locked_selected],
    Coefficient_per_1SD = locked_coefficients[locked_selected, 1],
    stringsAsFactors = FALSE
  )
  intercept <- locked_coefficients["(Intercept)", 1]
  x_train_locked <- sweep(sweep(x_qc[, locked_selected, drop = FALSE], 2, locked_result$means[locked_selected], "-"), 2, locked_result$standard_deviations[locked_selected], "/")
  x_validation_locked <- sweep(sweep(x_validation_all[, locked_selected, drop = FALSE], 2, locked_result$means[locked_selected], "-"), 2, locked_result$standard_deviations[locked_selected], "/")
  train_score <- as.numeric(intercept + x_train_locked %*% locked_coefficients[locked_selected, 1])
  validation_score <- as.numeric(intercept + x_validation_locked %*% locked_coefficients[locked_selected, 1])
  train_probability <- stats::plogis(train_score)
  validation_probability <- stats::plogis(validation_score)
  prediction_table$Rad_score <- c(train_score, validation_score)
  prediction_table$Predicted_probability <- c(train_probability, validation_probability)
  train_auc <- calculate_auc(training$outcome, train_probability)
  validation_auc <- calculate_auc(validation$outcome, validation_probability)
  performance_table <- rbind(
    data.frame(Dataset = "Training apparent (selection derived here)", N = length(training$outcome), Events = sum(training$outcome), t(train_auc), stringsAsFactors = FALSE),
    data.frame(Dataset = "Internal validation (untouched during selection)", N = length(validation$outcome), Events = sum(validation$outcome), t(validation_auc), stringsAsFactors = FALSE)
  )
}

cascade <- data.frame(
  Step = c("Extracted features", "After quality filter", "Univariable screen", "After Spearman pruning", "LASSO lambda.1se", "Stable features", "Locked Rad-score"),
  Feature_count = c(ncol(x_all), ncol(x_qc), length(screened_features), length(final_candidate_features), length(final_lasso$selected), length(stable_features), if (is.null(locked_result)) 0 else nrow(final_radscore)),
  stringsAsFactors = FALSE
)

# Publication-ready selection figures, exported via ASCII temporary paths.
cascade_plot <- ggplot(cascade, aes(x = reorder(Step, Feature_count), y = Feature_count)) +
  geom_col(fill = "#397A9E", width = 0.72) +
  geom_text(aes(label = Feature_count), hjust = -0.18, size = 3.6) +
  coord_flip(clip = "off") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(title = "Feature-selection cascade", x = NULL, y = "Feature count") +
  theme_minimal(base_family = "Helvetica", base_size = 11) +
  theme(plot.title = element_text(face = "bold"), panel.grid.major.y = element_blank(), plot.margin = margin(5.5, 24, 5.5, 5.5))

cv_data <- data.frame(log_lambda = log(final_lasso$cv_fit$lambda), cvm = final_lasso$cv_fit$cvm, cvsd = final_lasso$cv_fit$cvsd)
lasso_plot <- ggplot(cv_data, aes(x = log_lambda, y = cvm)) +
  geom_errorbar(aes(ymin = cvm - cvsd, ymax = cvm + cvsd), width = 0.02, colour = "grey75") +
  geom_line(colour = "#1F4E78", linewidth = 0.7) +
  geom_vline(xintercept = log(final_lasso$cv_fit$lambda.1se), colour = "#D55E00", linewidth = 0.75, linetype = "dashed") +
  geom_vline(xintercept = log(final_lasso$cv_fit$lambda.min), colour = "#6A9A3D", linewidth = 0.75, linetype = "dotted") +
  labs(title = "Training-only LASSO tuning", subtitle = "Dashed: lambda.1se; dotted: lambda.min", x = "log(lambda)", y = "Cross-validated binomial deviance") +
  theme_minimal(base_family = "Helvetica", base_size = 11) +
  theme(plot.title = element_text(face = "bold"))

top_stability <- stability[seq_len(min(25L, nrow(stability))), , drop = FALSE]
stability_plot <- ggplot(top_stability, aes(x = reorder(Feature, Selection_frequency), y = Selection_frequency, fill = ROI)) +
  geom_col(width = 0.72) +
  geom_hline(yintercept = selection_frequency_cutoff, linetype = "dashed", colour = "#D55E00") +
  coord_flip() +
  scale_fill_manual(values = c(Stomach = "#3B7EA1", Tumor = "#B65D46")) +
  scale_y_continuous(limits = c(0, 1), labels = function(x) paste0(round(100 * x), "%")) +
  labs(title = "Selection stability", x = NULL, y = "Selection frequency", fill = "ROI") +
  theme_minimal(base_family = "Helvetica", base_size = 10) +
  theme(plot.title = element_text(face = "bold"), panel.grid.major.y = element_blank())

combined_plot <- (cascade_plot | lasso_plot) / stability_plot + plot_annotation(tag_levels = "A")
temporary_png <- file.path(tempdir(), "dual_roi_lasso_selection.png")
temporary_pdf <- file.path(tempdir(), "dual_roi_lasso_selection.pdf")
ggsave(temporary_png, combined_plot, width = 14, height = 12, dpi = 600, bg = "white")
ggsave(temporary_pdf, combined_plot, width = 14, height = 12, device = grDevices::pdf, bg = "white")
file.copy(temporary_png, file.path(output_dir, paste0("Figure_1_", file_tag, "影像组学_LASSO稳定性筛选_600dpi.png")), overwrite = TRUE)
file.copy(temporary_pdf, file.path(output_dir, paste0("Figure_1_", file_tag, "影像组学_LASSO稳定性筛选.pdf")), overwrite = TRUE)

workbook <- openxlsx::createWorkbook()
add_sheet <- function(name, data) {
  openxlsx::addWorksheet(workbook, name)
  openxlsx::writeData(workbook, name, data, withFilter = TRUE)
  header_style <- openxlsx::createStyle(fgFill = "#1F4E78", fontColour = "#FFFFFF", textDecoration = "bold", halign = "center")
  openxlsx::addStyle(workbook, name, header_style, rows = 1, cols = seq_len(ncol(data)), gridExpand = TRUE)
  openxlsx::freezePane(workbook, name, firstRow = TRUE)
  openxlsx::setColWidths(workbook, name, cols = seq_len(ncol(data)), widths = "auto")
}

protocol <- data.frame(
  Item = c(
    "Analysis unit", "Outcome", "Training cohort", "Validation cohort", "Feature pool", "Quality filter", "Univariable screen", "Correlation pruning", "Primary LASSO", "Stability resampling", "Stable-feature rule", "Locked Rad-score", "Validation boundary", "Interpretation"
  ),
  Value = c(
    "Patient-level", "Anastomotic leak (1=yes, 0=no)",
    "574 patients; 44 leaks and 530 non-leaks", "246 patients; 19 leaks and 227 non-leaks",
    paste0(ncol(x_all), " features available after the prespecified technical-stability manifest; analysis ROI = ", roi_mode),
    "Remove non-finite, constant and near-zero variance features in training data",
    paste0("Univariable logistic regression; raw P < ", univariable_screen_p, "; BH-FDR reported for transparency"),
    paste0("Greedy Spearman pruning; |rho| > ", spearman_redundancy_cutoff, "; retain lower-P feature"),
    paste0("LASSO logistic regression with training-only stratified ", outer_folds, "-fold CV; lambda.1se is primary"),
    paste0(resample_count, " stratified ", round(100 * resample_fraction), "% subsamples; each repeats score-test screen, Spearman pruning and ", outer_folds, "-fold-CV LASSO"),
    paste0("Selection frequency >= ", selection_frequency_cutoff, " and coefficient direction consistency >= ", direction_consistency_cutoff),
    "Refit LASSO only with stable features on all training patients; standardized using training means/SDs",
    "Validation set not used for any feature screen, correlation filter, scaling, lambda choice, stability count or model selection",
    "Selected features are predictive model components, not causal risk factors"
  ), stringsAsFactors = FALSE
)
add_sheet("Protocol", protocol)
add_sheet("Feature_QC", feature_qc)
add_sheet("Univariable_all_features", univariable)
add_sheet("Univariable_screened", univariable[match(screened_features, univariable$Feature), , drop = FALSE])
add_sheet("Spearman_pruning_audit", pruned$audit)
add_sheet("LASSO_full_training", final_lasso_table)
add_sheet("Resample_run_log", resample_log)
add_sheet("Stability_selection", stability)
add_sheet("Locked_RadScore", if (nrow(final_radscore)) final_radscore else data.frame(Status = "No feature met the prespecified stability rule; no Rad-score was forced."))
add_sheet("Predictions", prediction_table)
add_sheet("Validation_metrics", if (nrow(performance_table)) performance_table else data.frame(Status = "No locked Rad-score was created under the prespecified stability rule."))
add_sheet("Selection_cascade", cascade)

temporary_xlsx <- file.path(tempdir(), "dual_roi_lasso_stability.xlsx")
openxlsx::saveWorkbook(workbook, temporary_xlsx, overwrite = TRUE)
output_xlsx <- file.path(output_dir, paste0(file_tag, "影像组学_训练集LASSO稳定性筛选.xlsx"))
file.copy(temporary_xlsx, output_xlsx, overwrite = TRUE)

write.csv(prediction_table, file.path(output_dir, paste0(file_tag, "_RadScore_训练验证预测.csv")), row.names = FALSE, fileEncoding = "UTF-8")
write.csv(final_radscore, file.path(output_dir, paste0(file_tag, "_RadScore_锁定系数.csv")), row.names = FALSE, fileEncoding = "UTF-8")

run_notes <- c(
  paste0(analysis_label, " radiomics selection completed using only the 574-patient training cohort."),
  paste0("Features: technical-stability eligible=", ncol(x_all), "; quality-retained=", ncol(x_qc), "; screened=", length(screened_features), "; Spearman-retained=", length(final_candidate_features), "; full-training LASSO lambda.1se=", length(final_lasso$selected), "; stable=", length(stable_features), "; locked Rad-score=", if (is.null(locked_result)) 0 else nrow(final_radscore)),
  paste0("Selection: raw univariable P<", univariable_screen_p, "; |Spearman rho|>", spearman_redundancy_cutoff, "; ", outer_folds, "-fold CV with lambda.1se primary."),
  paste0("Stability: ", resample_count, " stratified subsamples at ", resample_fraction, "; frequency cutoff=", selection_frequency_cutoff, "; direction cutoff=", direction_consistency_cutoff, "."),
  "Validation was excluded from all selection and preprocessing steps and was used only after the final Rad-score was locked.",
  "Because feature selection occurs in the training cohort, apparent training AUC is not an unbiased performance estimate; the validation AUC is the relevant held-out estimate.",
  if (screen$fallback_used) "The prespecified top-20 fallback was used because raw P<0.05 retained fewer than 20 features." else "The top-20 fallback was not used."
)
writeLines(run_notes, file.path(output_dir, "运行说明与方法边界.txt"))

current_script <- normalizePath("scripts/select_dual_roi_radiomics_lasso_stability.R")
file.copy(current_script, file.path(output_dir, basename(current_script)), overwrite = TRUE)
file.copy(current_script, file.path(code_archive_dir, basename(current_script)), overwrite = TRUE)
temporary_session <- file.path(tempdir(), "radiomics_lasso_session_info.txt")
writeLines(capture.output(sessionInfo()), temporary_session)
file.copy(temporary_session, file.path(output_dir, "R_sessionInfo.txt"), overwrite = TRUE)
file.copy(temporary_session, file.path(project_root, "代码与复现说明/environment/RADIOMICS_LASSO_STABILITY_R_SESSION_INFO.txt"), overwrite = TRUE)

print(list(
  output = output_dir,
  roi = roi_mode,
  extracted = ncol(x_all),
  quality_retained = ncol(x_qc),
  screened = length(screened_features),
  correlation_retained = length(final_candidate_features),
  lasso_lambda_1se = length(final_lasso$selected),
  stable_features = length(stable_features),
  locked_radscore_features = if (is.null(locked_result)) 0L else nrow(final_radscore)
))
