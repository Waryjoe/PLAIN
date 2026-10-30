#!/usr/bin/env Rscript
source("scripts/plain_paths.R")

# Dedicated R VIF/GVIF workbook for training-set priority candidates.
# R 4.6 on this computer runs with the C locale. To avoid its Chinese-sheet
# encoding limitation, source workbooks are temporarily copied to ASCII paths
# and accessed by sheet number. All VIF/GVIF calculations occur in local R.

options(stringsAsFactors = FALSE)
user_lib <- plain_path("R_LIBRARY", "")
if (dir.exists(user_lib)) .libPaths(c(user_lib, .libPaths()))
suppressPackageStartupMessages({ library(readxl); library(openxlsx) })

training_input <- plain_path("PROJECT", "汇总/训练验证集/训练集_70%_574例.xlsx")
candidate_input <- paste0(
  plain_path("PROJECT", "汇总/"),
  "训练集单因素候选变量筛选_排除肿瘤部位/",
  "训练集_单因素候选变量筛选_排除肿瘤部位.xlsx"
)
output_dir <- paste0(
  plain_path("PROJECT", "汇总/"),
  "训练集优先候选变量_相关性与VIF_排除肿瘤部位"
)
output_file <- file.path(output_dir, "训练集_优先候选变量_VIF数值汇总_R.xlsx")

vif_label <- function(value) {
  if (!is.finite(value)) return("Not estimable")
  if (value < 2) return("Low collinearity (VIF < 2)")
  if (value < 5) return("Mild-to-moderate collinearity (2 <= VIF < 5)")
  if (value < 10) return("Needs attention (5 <= VIF < 10)")
  "High collinearity (VIF >= 10)"
}
gvif_label <- function(value) {
  if (!is.finite(value)) return("Not estimable")
  if (value < 1.4) return("Low collinearity (adjusted GVIF < 1.4)")
  if (value < 2) return("Needs attention (1.4 <= adjusted GVIF < 2)")
  "High collinearity (adjusted GVIF >= 2)"
}
safe_det <- function(value) {
  value <- as.matrix(value)
  if (nrow(value) == 0) return(1)
  result <- tryCatch(det(value), error = function(error) NA_real_)
  if (!is.finite(result) || result <= 0) return(NA_real_)
  as.numeric(result)
}
calculate_vif <- function(design, index) {
  outcome <- design[, index]
  predictors <- design[, -index, drop = FALSE]
  fit <- stats::lm.fit(cbind(intercept = 1, predictors), outcome)
  rss <- sum(fit$residuals^2)
  tss <- sum((outcome - mean(outcome))^2)
  if (tss <= 0) return(NA_real_)
  1 / max(1 - (1 - rss / tss), .Machine$double.eps)
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
temporary_training <- file.path(tempdir(), "training.xlsx")
temporary_candidates <- file.path(tempdir(), "candidates.xlsx")
temporary_output <- file.path(tempdir(), "priority_vif.xlsx")
stopifnot(file.copy(training_input, temporary_training, overwrite = TRUE))
stopifnot(file.copy(candidate_input, temporary_candidates, overwrite = TRUE))

# Sheet 1 is clinical data in the training workbook; sheet 2 is the candidate list.
clinical <- as.data.frame(readxl::read_excel(temporary_training, sheet = 1))
raw_candidates <- as.data.frame(readxl::read_excel(temporary_candidates, sheet = 2))
if (ncol(raw_candidates) < 7) stop("Candidate workbook does not have expected columns")

# Column positions are fixed by the preceding screening output. Unicode escapes
# make the continuous-type test safe under R's C locale.
candidates <- data.frame(
  variable = raw_candidates[[1]],
  variable_type = raw_candidates[[2]],
  p_value = suppressWarnings(as.numeric(raw_candidates[[5]])),
  stringsAsFactors = FALSE
)
candidates <- candidates[is.finite(candidates$p_value) & candidates$p_value < 0.05, , drop = FALSE]
if (nrow(candidates) != 9) stop("Expected 9 P<0.05 priority candidates")
continuous_tag <- "\u8fde\u7eed"
continuous <- candidates$variable[candidates$variable_type == continuous_tag]
categorical <- candidates$variable[candidates$variable_type != continuous_tag]
if (length(continuous) != 7 || length(categorical) != 2) stop("Unexpected candidate type composition")
if (!all(candidates$variable %in% names(clinical))) stop("Training workbook is missing candidate variables")

# Apply the locked analysis rule: exact zero in continuous measurements is missing.
missing_qc <- data.frame(
  Variable = character(), Type = character(), Zero_replaced_as_missing = integer(),
  Missing_after_cleaning = integer(), Missing_rate = numeric()
)
for (variable in continuous) {
  values <- suppressWarnings(as.numeric(clinical[[variable]]))
  zero_count <- sum(values == 0, na.rm = TRUE)
  values[values == 0] <- NA_real_
  clinical[[variable]] <- values
  missing_qc <- rbind(missing_qc, data.frame(
    Variable = variable, Type = "Continuous", Zero_replaced_as_missing = zero_count,
    Missing_after_cleaning = sum(is.na(values)), Missing_rate = mean(is.na(values))
  ))
}
for (variable in categorical) {
  missing_qc <- rbind(missing_qc, data.frame(
    Variable = variable, Type = "Categorical", Zero_replaced_as_missing = 0L,
    Missing_after_cleaning = sum(is.na(clinical[[variable]])),
    Missing_rate = mean(is.na(clinical[[variable]]))
  ))
}

# Pairwise Spearman correlation for continuous candidates.
correlation_rows <- data.frame(
  Variable_1 = character(), Variable_2 = character(), Complete_cases = integer(),
  Spearman_rho = numeric(), P_value = numeric(), Flag_abs_rho_ge_0_70 = logical()
)
for (pair in combn(continuous, 2, simplify = FALSE)) {
  values <- clinical[, pair, drop = FALSE]
  values <- values[stats::complete.cases(values), , drop = FALSE]
  result <- suppressWarnings(stats::cor.test(values[[1]], values[[2]], method = "spearman", exact = FALSE))
  correlation_rows <- rbind(correlation_rows, data.frame(
    Variable_1 = pair[1], Variable_2 = pair[2], Complete_cases = nrow(values),
    Spearman_rho = unname(result$estimate), P_value = result$p.value,
    Flag_abs_rho_ge_0_70 = abs(unname(result$estimate)) >= 0.70
  ))
}
correlation_rows <- correlation_rows[order(abs(correlation_rows$Spearman_rho), decreasing = TRUE), , drop = FALSE]

# Create a diagnostic-only, imputed design matrix. No outcome is used here.
design <- data.frame(row.names = seq_len(nrow(clinical)))
term_indices <- list()
imputation <- data.frame(
  Variable = character(), Type = character(), Diagnostic_imputation = character(),
  Fill_value_or_reference = character(), Imputed_n = integer()
)
column_map <- data.frame(
  Design_column = character(), Variable = character(), Type = character(), Level = character()
)
for (variable in continuous) {
  values <- clinical[[variable]]
  fill_value <- stats::median(values, na.rm = TRUE)
  design[[variable]] <- ifelse(is.na(values), fill_value, values)
  term_indices[[variable]] <- ncol(design)
  imputation <- rbind(imputation, data.frame(
    Variable = variable, Type = "Continuous", Diagnostic_imputation = "Training median",
    Fill_value_or_reference = as.character(fill_value), Imputed_n = sum(is.na(values))
  ))
  column_map <- rbind(column_map, data.frame(
    Design_column = variable, Variable = variable, Type = "Continuous", Level = ""
  ))
}
for (variable in categorical) {
  values <- clinical[[variable]]
  observed <- values[!is.na(values)]
  mode_value <- names(sort(table(observed), decreasing = TRUE))[1]
  filled <- as.character(values)
  filled[is.na(filled)] <- mode_value
  levels <- sort(unique(filled))
  reference <- levels[1]
  indices <- integer()
  for (level in levels[-1]) {
    column_name <- paste0(variable, "=", level)
    design[[column_name]] <- as.numeric(filled == level)
    indices <- c(indices, ncol(design))
    column_map <- rbind(column_map, data.frame(
      Design_column = column_name, Variable = variable, Type = "Categorical", Level = level
    ))
  }
  term_indices[[variable]] <- indices
  imputation <- rbind(imputation, data.frame(
    Variable = variable, Type = "Categorical", Diagnostic_imputation = "Training mode",
    Fill_value_or_reference = paste0("Mode=", mode_value, "; reference=", reference),
    Imputed_n = sum(is.na(values))
  ))
}

design_matrix <- as.matrix(design)
if (any(!is.finite(design_matrix))) stop("Non-finite value remained in VIF design matrix")

column_vif <- data.frame(
  Design_column = names(design), VIF = numeric(ncol(design)), Interpretation = character(ncol(design))
)
for (index in seq_len(ncol(design))) {
  value <- calculate_vif(design_matrix, index)
  column_vif$VIF[index] <- value
  column_vif$Interpretation[index] <- vif_label(value)
}
column_vif <- merge(column_vif, column_map, by = "Design_column", all.x = TRUE, sort = FALSE)
column_vif <- column_vif[order(column_vif$VIF, decreasing = TRUE), , drop = FALSE]

correlation_matrix <- stats::cor(design_matrix)
full_det <- safe_det(correlation_matrix)
vif_summary <- data.frame(
  Variable = character(), Type = character(), df = integer(), Design_columns = character(),
  VIF_or_GVIF = numeric(), Adjusted_GVIF_1_over_2df = numeric(),
  Interpretation = character(), Recommendation = character()
)
for (variable in candidates$variable) {
  indices <- term_indices[[variable]]
  others <- setdiff(seq_len(ncol(design)), indices)
  gvif <- safe_det(correlation_matrix[indices, indices, drop = FALSE]) *
    safe_det(correlation_matrix[others, others, drop = FALSE]) / full_det
  degrees <- length(indices)
  adjusted_gvif <- gvif^(1 / (2 * degrees))
  variable_type <- candidates$variable_type[match(variable, candidates$variable)]
  is_continuous <- variable_type == continuous_tag
  value_for_label <- if (is_continuous) gvif else adjusted_gvif
  good <- if (is_continuous) gvif < 5 else adjusted_gvif < 1.4
  vif_summary <- rbind(vif_summary, data.frame(
    Variable = variable,
    Type = if (is_continuous) "Continuous" else "Categorical",
    df = degrees,
    Design_columns = paste(names(design)[indices], collapse = "; "),
    VIF_or_GVIF = gvif,
    Adjusted_GVIF_1_over_2df = adjusted_gvif,
    Interpretation = if (is_continuous) vif_label(value_for_label) else gvif_label(value_for_label),
    Recommendation = if (good) "Retain for next modelling step" else "Review for clinical prioritisation"
  ))
}
vif_summary <- vif_summary[order(vif_summary$Adjusted_GVIF_1_over_2df, decreasing = TRUE), , drop = FALSE]

notes <- data.frame(
  Item = c(
    "Analysis set", "Variables included", "Variables excluded", "Continuous zero rule",
    "Correlation", "VIF calculation", "Categorical variables", "Diagnostic imputation",
    "Validation set", "R software", "Interpretation boundary"
  ),
  Detail = c(
    paste0("Training set only: ", nrow(clinical), " patients; validation set not read"),
    paste(candidates$variable, collapse = "; "),
    "Extended candidates, tumour location, and TG2 were not included",
    "All exact zero values in continuous measurements were treated as missing",
    "Pairwise complete-case Spearman correlations for continuous variables",
    "Classical VIF = 1/(1-R-squared); continuous variables report VIF",
    "Reference-cell dummy coding; categorical variables report GVIF and adjusted GVIF",
    "Training median for continuous variables and training mode for categorical variables, only for VIF design matrix",
    "Not read or used",
    paste0("R ", getRversion(), "; readxl ", packageVersion("readxl"), "; openxlsx ", packageVersion("openxlsx")),
    "VIF/GVIF describes predictor redundancy only; it does not establish independent association with leakage"
  )
)

workbook <- openxlsx::createWorkbook()
header_style <- openxlsx::createStyle(fgFill = "#1F4E78", fontColour = "#FFFFFF", textDecoration = "bold", halign = "center", valign = "center", wrapText = TRUE)
add_sheet <- function(name, data) {
  openxlsx::addWorksheet(workbook, name)
  openxlsx::writeData(workbook, name, data, withFilter = TRUE)
  openxlsx::addStyle(workbook, name, header_style, rows = 1, cols = seq_len(ncol(data)), gridExpand = TRUE)
  openxlsx::freezePane(workbook, name, firstRow = TRUE)
  openxlsx::setColWidths(workbook, name, cols = seq_len(ncol(data)), widths = "auto")
}
add_sheet("VIF_summary", vif_summary)
add_sheet("Dummy_variable_VIF", column_vif)
add_sheet("Spearman_correlation", correlation_rows)
add_sheet("Priority_candidates", candidates)
add_sheet("VIF_imputation", imputation)
add_sheet("Missingness_QC", missing_qc)
add_sheet("R_method_notes", notes)
openxlsx::saveWorkbook(workbook, temporary_output, overwrite = TRUE)
stopifnot(file.copy(temporary_output, output_file, overwrite = TRUE))

writeLines(c(
  paste0("output_file=", output_file),
  paste0("n_priority_candidates=", nrow(candidates)),
  paste0("max_VIF_or_GVIF=", max(vif_summary$VIF_or_GVIF, na.rm = TRUE)),
  paste0("R_version=", getRversion())
))
