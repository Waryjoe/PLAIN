#!/usr/bin/env Rscript
source("scripts/plain_paths.R")
options(stringsAsFactors = FALSE)
.libPaths(c(plain_path("R_LIBRARY", ""), .libPaths()))
suppressPackageStartupMessages({
  library(lightgbm)
  library(readxl)
  library(openxlsx)
  library(ggplot2)
  library(ggbeeswarm)
  library(svglite)
})

root <- plain_path("PROJECT", "汇总")
out <- file.path(root, "figure/figure汇总/Radiomics_SHAP_20260911")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out, "Reproducibility"), showWarnings = FALSE)
model_path <- file.path(root, "最终模型_10特征双ROI_LightGBM_全流程评估图表代码包_20260812/01_Tables/影像模型_训练集全量拟合对象_R.rds")
source_path <- file.path(root, "论文Figure_最终整理_10特征LightGBM_内外部验证_20260828/03_Source_Data/Final_Figure_Source_Data_10Feature_LightGBM_20260828.xlsx")
train_path <- file.path(root, "影像组学_训练验证集_v5_GRAPE对齐/影像组学v5_训练集_70%_574例_3190特征.xlsx")
valid_path <- file.path(root, "影像组学_训练验证集_v5_GRAPE对齐/影像组学v5_验证集_30%_246例_3190特征.xlsx")
input_paths <- c(model_path, source_path, train_path, valid_path)
input_md5 <- tools::md5sum(input_paths)
checks <- data.frame(Check = character(), Value = character(), Status = character())
record <- function(label, value, pass = TRUE) {
  checks <<- rbind(checks, data.frame(Check = label, Value = as.character(value), Status = ifelse(pass, "PASS", "FAIL")))
  if (!pass) stop(label, ": ", value)
}

message("Loading saved locked model; no refitting.")
fits <- readRDS(model_path)
wrapper <- fits[["10__LightGBM"]]
e <- environment(wrapper$predict)
booster <- e$object
x <- as.matrix(e$x)
y <- as.integer(e$y == "Leak")
selected <- read.xlsx(source_path, sheet = "Selected_10_features")
selected <- selected[order(selected$Rank), ]
record("Model feature schema", paste(dim(x), collapse = " x "), identical(dim(x), c(574L, 10L)))
record("Training events", sum(y), sum(y) == 44L)
record("Feature order matches locked figure source", "Exact", identical(colnames(x), make.names(selected$Feature, unique = TRUE)))
selected$Model_column <- colnames(x)
selected$ID <- paste0(ifelse(selected$ROI == "Stomach", "S", "T"), ave(seq_len(10), selected$ROI, FUN = seq_along))
short <- c(S1 = "90th percentile", S2 = "Correlation", S3 = "Long-run low gray-level emphasis",
           S4 = "Mean", S5 = "Small-area low gray-level emphasis", T1 = "Major axis length",
           T2 = "Minimum", T3 = "Skewness", T4 = "Coarseness", T5 = "Large-dependence low gray-level emphasis")
selected$Short_name <- unname(short[selected$ID])

read_selected <- function(path) {
  headers <- names(read_xlsx(path, n_max = 0))
  types <- rep("skip", length(headers))
  types[match(selected$Feature, headers)] <- "numeric"
  types[1:3] <- c("text", "text", "numeric")
  d <- as.data.frame(read_xlsx(path, col_types = types))
  mat <- as.matrix(d[, selected$Feature, drop = FALSE])
  colnames(mat) <- selected$Model_column
  list(x = mat, y = as.integer(d[[3]]), ids = d[[1]])
}
train <- read_selected(train_path)
valid <- read_selected(valid_path)
record("Training matrix equals stored model input", max(abs(train$x - x)), isTRUE(all.equal(unname(train$x), unname(x), tolerance = 1e-12)))
record("Training outcome order", "Exact", identical(train$y, y))
record("Internal validation sample size", nrow(valid$x), nrow(valid$x) == 246L)
record("Internal validation events", sum(valid$y), sum(valid$y) == 19L)
mu <- colMeans(x)
sigma <- apply(x, 2, sd)
sigma[!is.finite(sigma) | sigma < 1e-10] <- 1
standardize <- function(m) sweep(sweep(m, 2, mu, "-"), 2, sigma, "/")
z <- standardize(x)
z_valid <- standardize(valid$x)
record("Saved preprocessing reproduced", max(abs(z - e$prepared$train)), max(abs(z - e$prepared$train)) < 1e-12)

gain <- lgb.importance(booster, percentage = TRUE)
old_gain <- read.xlsx(source_path, "Radiomics_importance")
gain_delta <- max(abs(gain$Gain - old_gain$Gain[match(gain$Feature, old_gain$Feature)]))
record("Gain agrees with original panel D", format(gain_delta, scientific = TRUE), is.finite(gain_delta) && gain_delta < 1e-10)
pred_source <- read.xlsx(source_path, "Internal_predictions")
pv_saved <- pred_source[pred_source$Dataset == "Internal validation", ]
pv <- as.numeric(predict(booster, z_valid, type = "response"))
record("Validation outcome order matches source", "Exact", identical(valid$y, as.integer(pv_saved$Outcome)))
pred_delta <- max(abs(pv - pv_saved$Radiomics_10Feature_LightGBM))
record("Saved internal predictions reproduced", format(pred_delta, scientific = TRUE), pred_delta < 1e-10)
valid_auc <- as.numeric(pROC::auc(pROC::roc(valid$y, pv, levels = c(0, 1), direction = "<", quiet = TRUE)))
record("Internal validation AUC reproduced", sprintf("%.9f", valid_auc), round(valid_auc, 3) == 0.929)

message("Computing native Tree SHAP and verifying raw-score additivity.")
calculate <- function(m, cohort) {
  contributions <- as.matrix(predict(booster, m, type = "contrib"))
  raw <- as.numeric(predict(booster, m, type = "raw"))
  p <- as.numeric(predict(booster, m, type = "response"))
  record(paste(cohort, "SHAP dimensions"), paste(dim(contributions), collapse = " x "), identical(dim(contributions), c(nrow(m), 11L)))
  record(paste(cohort, "finite SHAP"), "All finite", all(is.finite(contributions)))
  error <- max(abs(rowSums(contributions) - raw))
  record(paste(cohort, "SHAP plus base equals raw log-odds"), format(error, scientific = TRUE), error < 1e-8)
  error <- max(abs(plogis(rowSums(contributions)) - p))
  record(paste(cohort, "sigmoid reconstructs prediction"), format(error, scientific = TRUE), error < 1e-10)
  colnames(contributions) <- c(selected$ID, "Base_value")
  list(contributions = contributions, raw = raw, p = p)
}
sh <- calculate(z, "Training")
sv <- calculate(z_valid, "Internal validation")
selected$Mean_abs_SHAP <- colMeans(abs(sh$contributions[, 1:10]))
selected$Mean_abs_SHAP_internal <- colMeans(abs(sv$contributions[, 1:10]))
selected$Gain <- gain$Gain[match(selected$Model_column, gain$Feature)]
selected$SHAP_rank <- rank(-selected$Mean_abs_SHAP, ties.method = "first")
selected$Training_mean <- mu
selected$Training_SD <- sigma
selected$Color_low_q05 <- apply(x, 2, quantile, probs = 0.05)
selected$Color_high_q95 <- apply(x, 2, quantile, probs = 0.95)
long <- do.call(rbind, lapply(seq_len(10), function(j) {
  v <- x[, j]
  lo <- selected$Color_low_q05[j]
  hi <- selected$Color_high_q95[j]
  color <- if (hi > lo) pmin(pmax((v - lo) / (hi - lo), 0), 1) else rep(0.5, length(v))
  data.frame(Sample = sprintf("TR%04d", seq_len(nrow(x))), ID = selected$ID[j], ROI = selected$ROI[j],
    Raw_feature_value = v, Standardized_value = z[, j], SHAP_log_odds = sh$contributions[, j], Color_value_0_1 = color)
}))
feature_order <- selected$ID[order(selected$Mean_abs_SHAP)]
long$ID <- factor(long$ID, levels = feature_order)
set.seed(20260911)
long <- long[sample.int(nrow(long)), ]

make_plot <- function(named = FALSE, panel_letter = TRUE) {
  labels <- if (named) setNames(paste0(selected$ID, "  ", selected$Short_name), selected$ID) else setNames(selected$ID, selected$ID)
  p <- ggplot(long, aes(SHAP_log_odds, ID)) +
    geom_hline(yintercept = seq_len(10), linewidth = 0.22, color = "#E4E7EB", linetype = "dotted") +
    geom_vline(xintercept = 0, linewidth = 0.38, color = "#9BA3AC") +
    geom_quasirandom(aes(color = Color_value_0_1, group = ID), orientation = "y", width = 0.32,
                    bandwidth = 0.6, size = if (named) 0.95 else 0.76, alpha = 0.8, stroke = 0) +
    scale_y_discrete(labels = labels, expand = expansion(add = 0.6)) +
    scale_x_continuous(breaks = scales::breaks_pretty(n = 5), expand = expansion(mult = c(0.04, 0.04))) +
    scale_color_gradientn(colors = c("#1686D2", "#8173BA", "#D72D72"), limits = c(0, 1),
       breaks = c(0, 1), labels = c("Low", "High"), name = "Feature value",
       guide = guide_colorbar(display = "rectangles", nbin = 100, title.position = "right", title.theme = element_text(angle = 90, size = 8.5, hjust = 0.5, vjust = 0.5),
                              barwidth = grid::unit(2.2, "mm"), barheight = grid::unit(37, "mm"), ticks = FALSE)) +
    labs(x = "SHAP value (log-odds)", y = NULL, tag = if (panel_letter) "D" else NULL) +
    theme_classic(base_size = 9.5, base_family = "Arial") +
    theme(axis.line.y = element_blank(), axis.line.x = element_line(linewidth = 0.4),
          axis.ticks.y = element_blank(), axis.ticks.x = element_line(linewidth = 0.3),
          axis.text = element_text(color = "#17212B", size = 9.5), axis.text.y = element_text(margin = margin(r = 4)),
          axis.title.x = element_text(margin = margin(t = 6), size = 10),
          legend.text = element_text(size = 8), legend.spacing.x = grid::unit(1, "mm"),
          legend.margin = margin(0, 0, 0, 0), legend.box.spacing = grid::unit(2, "mm"),
          plot.tag = element_text(size = 15, face = "plain", color = "#111111"),
          plot.tag.position = "topleft", plot.margin = margin(t = 6, r = 4, b = 4, l = 5))
  p
}
message("Writing editable vector figures and numerical source data.")
export_plot <- function(p, stem, width, height) {
  ggsave(file.path(out, paste0(stem, ".svg")), p, device = svglite, width = width, height = height, units = "mm", bg = "white")
  ggsave(file.path(out, paste0(stem, ".pdf")), p + theme(text = element_text(family = "Helvetica")),
         device = grDevices::pdf, width = width, height = height, units = "mm", bg = "white", useDingbats = FALSE)
  ggsave(file.path(out, paste0(stem, ".png")), p, device = ragg::agg_png, width = width, height = height, units = "mm", dpi = 400, bg = "white")
}
export_plot(make_plot(), "D_Radiomics_SHAP_Editable", 94, 78)
export_plot(make_plot(panel_letter = FALSE), "Radiomics_SHAP_WithoutPanelLetter", 94, 78)
export_plot(make_plot(named = TRUE, panel_letter = FALSE), "Radiomics_SHAP_NamedFeatures", 175, 100)

pred_frame <- function(s, prefix, outcome) data.frame(Sample = sprintf(paste0(prefix, "%04d"), seq_along(outcome)), Outcome = outcome,
    Probability_leak = s$p, Raw_log_odds = s$raw, s$contributions, check.names = FALSE)
raw_frame <- function(m, prefix) {
  colnames(m) <- selected$ID
  data.frame(Sample = sprintf(paste0(prefix, "%04d"), seq_len(nrow(m))), m, check.names = FALSE)
}
notes <- data.frame(Item = c("Model", "Main plot cohort", "Evaluation type", "Method", "Scale", "Positive SHAP", "Feature order", "Dot color", "Preprocessing", "Validation check", "Caveat", "Identifiers", "Software documentation"),
  Detail = c("Saved 10__LightGBM; 100 trees; num_leaves=7; learning_rate=0.05; no retraining",
    "Training cohort: 574 patients, 44 leakage events; all patients shown",
    "Explanations of the fixed full-training-fit model; NOT out-of-fold SHAP and NOT independent performance estimates",
    "LightGBM native Tree SHAP: predict(type='contrib'); tree-path-dependent baseline; no new reference cohort sampled",
    "Raw binary-model log-odds, not probability units; base + sum(SHAP) = raw score",
    "Increases the model's predicted leakage log-odds relative to the base value",
    "Descending training mean absolute SHAP; original S1-S5 and T1-T5 identifiers preserved",
    "Raw feature value normalized separately within each feature using training 5th/95th percentiles; color clipped only, SHAP values unaltered",
    "Exact saved training means/SDs; validation uses training parameters; no imputation needed",
    sprintf("246 patients / 19 events; original validation probabilities reproduced; AUC %.9f", valid_auc),
    "Model attribution, not causation or established biological mechanism. Correlated features can share attribution. Native path-dependent and interventional SHAP are different estimands.",
    "TR/IV row indices only; no hospital IDs or patient names in deliverables. Row order matches verified v5 input workbook.",
    "https://lightgbm.readthedocs.io/en/stable/R/reference/predict.lgb.Booster.html"))
record("Source files unchanged", "Input MD5 hashes identical before/after analysis", identical(input_md5, tools::md5sum(input_paths)))
wb <- createWorkbook()
tables <- list(Readme = notes, Feature_mapping = selected[order(selected$ROI, selected$ID), ],
  Global_importance = selected[order(selected$SHAP_rank), c("SHAP_rank", "ID", "ROI", "Short_name", "Feature", "Mean_abs_SHAP", "Mean_abs_SHAP_internal", "Gain")],
  Training_SHAP = pred_frame(sh, "TR", y), Training_feature_values = raw_frame(x, "TR"),
  Plot_source_long = long, Internal_SHAP = pred_frame(sv, "IV", valid$y), Internal_feature_values = raw_frame(valid$x, "IV"),
  Verification = checks, Input_provenance = data.frame(Path = input_paths, MD5 = unname(input_md5)))
for (nm in names(tables)) {
  addWorksheet(wb, nm)
  writeData(wb, nm, tables[[nm]], withFilter = TRUE)
  freezePane(wb, nm, firstRow = TRUE)
  addStyle(wb, nm, createStyle(fontColour = "#FFFFFF", fgFill = "#24647A", textDecoration = "bold"), rows = 1, cols = seq_len(ncol(tables[[nm]])), gridExpand = TRUE)
  setColWidths(wb, nm, cols = seq_len(ncol(tables[[nm]])), widths = 22)
}
setColWidths(wb, "Readme", 2, 115)
setColWidths(wb, "Verification", 1, 65)
setColWidths(wb, "Feature_mapping", match("Feature", names(tables$Feature_mapping)), 78)
saveWorkbook(wb, file.path(out, "Radiomics_SHAP_SourceData_and_QC.xlsx"), overwrite = TRUE)
saveRDS(list(selected = selected, training = sh, internal = sv, plot_source = long, checks = checks), file.path(out, "Reproducibility/Verified_SHAP_results.rds"))
lgb.save(booster, file.path(out, "Reproducibility/Locked_10Feature_LightGBM.txt"))
capture.output(sessionInfo(), file = file.path(out, "Reproducibility/R_sessionInfo.txt"))
legend_en <- paste("D. SHAP summary plot of the locked 10-feature dual-ROI radiomics LightGBM model in the training cohort (n = 574).",
 "Each dot represents one patient for one feature. Features are ordered by descending mean absolute SHAP value.",
 "Horizontal position indicates the SHAP contribution to the model's predicted log-odds of anastomotic leakage relative to its base value; positive values increase and negative values decrease the model output.",
 "Dot color indicates the relative feature value, from low (blue) to high (magenta), scaled separately for each feature using its training-cohort 5th and 95th percentiles and clipped for color display only.",
 "S1-S5 and T1-T5 denote the previously defined stomach-ROI and tumor-ROI features, respectively.",
 "Native Tree SHAP values were computed from the saved full-training-fit model and verified for additivity. These are in-sample model explanations, not out-of-fold estimates or causal effects.")
legend_zh <- paste("D. 最终锁定的十特征双ROI影像组学LightGBM模型在训练队列（n = 574）中的SHAP汇总图。",
 "每个点代表某位患者的一个特征，特征按平均绝对SHAP值由高至低排列。横轴为各特征相对基线值对模型预测吻合口漏对数优势（log-odds）的贡献，正值使模型输出升高，负值使其降低。",
 "颜色表示该特征取值的相对高低，蓝色为低值、洋红色为高值；每个特征分别以训练队列第5和第95百分位缩放，仅颜色进行截断，SHAP值不截断。",
 "S1-S5和T1-T5分别对应原图定义的胃ROI和肿瘤ROI特征。采用保存的训练集全量拟合模型计算原生Tree SHAP，并验证加和性；本图为训练样本内模型解释，并非折外解释、独立验证性能或因果效应。")
writeLines(c("# Radiomics SHAP: editable replacement for panel D", "", "## Files", "",
 "- D_Radiomics_SHAP_Editable.svg: recommended replacement for panel D (94 x 78 mm).",
 "- Radiomics_SHAP_WithoutPanelLetter.svg: identical plot without the D label.",
 "- Radiomics_SHAP_NamedFeatures.svg: enlarged reading version with feature names.",
 "- PDF and PNG previews accompany all SVGs. SVG dots and text are vector elements; no rasterized plot layer.",
 "- Radiomics_SHAP_SourceData_and_QC.xlsx: feature mapping, per-patient SHAP, feature values, global importance, and verification.",
 "- Reproducibility/: code, saved SHAP results, exported locked model, and software versions.", "",
 "## English Legend", "", legend_en, "", "## 中文图注", "", legend_zh, "",
 "## Verification", "", paste0("- ", checks$Check, ": ", checks$Value, " [", checks$Status, "]"), "",
 "## Method Reference", "", "https://lightgbm.readthedocs.io/en/stable/R/reference/predict.lgb.Booster.html",
 "https://shap.readthedocs.io/en/latest/example_notebooks/api_examples/plots/beeswarm.html", "",
 "Original Illustrator files and existing figures were not modified."), file.path(out, "README_and_FigureLegend.md"), useBytes = TRUE)
args <- commandArgs(trailingOnly = FALSE)
script <- sub("^--file=", "", args[startsWith(args, "--file=")][1])
if (!is.na(script)) file.copy(script, file.path(out, "Reproducibility/generate_radiomics_shap.R"), overwrite = TRUE)
helper <- file.path(dirname(script), "finalize_radiomics_shap_workbook.py")
file.copy(helper, file.path(out, "Reproducibility/finalize_radiomics_shap_workbook.py"), overwrite = TRUE)
python <- plain_path("PYTHON", "")
status <- system2(python, c(shQuote(helper), shQuote(out)))
if (status != 0L) stop("Workbook compatibility and readback check failed")
print(selected[order(selected$SHAP_rank), c("ID", "Short_name", "Mean_abs_SHAP", "Mean_abs_SHAP_internal", "Gain")])
print(checks)
message("DONE: ", out)
