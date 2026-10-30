#!/usr/bin/env Rscript
# Release utility: evaluate supplied frozen probabilities without refitting.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("Usage: Rscript scripts/evaluate_predictions.R INPUT.csv NEW_OUTPUT_DIR")
if (!requireNamespace("pROC", quietly = TRUE)) stop("Install the R package pROC first")
source("scripts/evaluation_helpers.R")
d <- read.csv(args[1], check.names = FALSE, colClasses = c(patient_id = "character", cohort = "character"))
validate_predictions(d)
out <- args[2]
if (file.exists(out)) stop("Output exists; choose a new directory")
dir.create(out, recursive = TRUE)
models <- c(Clinical = "clinical_probability", Radiomics = "radiomics_probability", Fusion = "fusion_probability")
colors <- c(Clinical = "#3D7FA3", Radiomics = "#269F96", Fusion = "#86628D")
thresholds <- seq(.01, .35, by = .005)
metrics <- curves <- calibration <- decisions <- list()
for (cohort in unique(d$cohort)) {
  s <- d[d$cohort == cohort, ]; y <- s$outcome
  for (model in names(models)) {
    p <- s[[models[[model]]]]; key <- paste(cohort, model)
    metrics[[key]] <- cbind(cohort, model, discrimination(y, p))
    r <- pROC::roc(y, p, levels = c(0, 1), direction = "<", quiet = TRUE)
    curves[[key]] <- data.frame(cohort, model, fpr = 1 - r$specificities, tpr = r$sensitivities)
    calibration[[key]] <- cbind(cohort, model, calibration_bins(y, p))
    decisions[[key]] <- data.frame(cohort, model, threshold = thresholds, net_benefit = net_benefit(y, p, thresholds))
  }
  decisions[[paste(cohort, "all")]] <- data.frame(cohort, model = "Treat all", threshold = thresholds,
      net_benefit = mean(y) - (1 - mean(y)) * thresholds / (1 - thresholds))
  decisions[[paste(cohort, "none")]] <- data.frame(cohort, model = "Treat none", threshold = thresholds, net_benefit = 0)
}
metrics <- do.call(rbind, metrics); curves <- do.call(rbind, curves)
calibration <- do.call(rbind, calibration); decisions <- do.call(rbind, decisions)
for (nm in c("metrics", "curves", "calibration", "decisions"))
  write.csv(get(nm), file.path(out, paste0(nm, ".csv")), row.names = FALSE)
pdf(file.path(out, "evaluation.pdf"), width = 12, height = 4, family = "Helvetica", useDingbats = FALSE)
for (cohort in unique(d$cohort)) {
  par(mfrow = c(1, 3), mar = c(4, 4, 3, 1))
  plot(NA, xlim = c(0, 1), ylim = c(0, 1), xlab = "1 - Specificity", ylab = "Sensitivity", main = paste(cohort, "ROC"))
  abline(0, 1, lty = 2, col = "gray70")
  for (model in names(models)) {
    z <- curves[curves$cohort == cohort & curves$model == model, ]
    lines(z$fpr, z$tpr, col = colors[model], lwd = 2)
  }
  legend("bottomright", names(models), col = colors, lty = 1, bty = "n", cex = .8)
  plot(NA, xlim = c(0, 1), ylim = c(0, 1), xlab = "Mean predicted probability", ylab = "Observed event rate", main = "Calibration")
  abline(0, 1, lty = 2, col = "gray70")
  for (model in names(models)) {
    z <- calibration[calibration$cohort == cohort & calibration$model == model, ]
    lines(z$mean_predicted, z$observed, type = "b", col = colors[model], pch = 1)
    segments(z$mean_predicted, z$low95, z$mean_predicted, z$high95, col = colors[model])
  }
  z <- decisions[decisions$cohort == cohort, ]
  limits <- range(c(0, z$net_benefit[z$model %in% names(models)]))
  if (diff(limits) == 0) limits <- limits + c(-.01, .01)
  plot(NA, xlim = range(thresholds), ylim = limits, xlab = "Threshold probability", ylab = "Net benefit", main = "Decision curves")
  abline(h = 0, lty = 3)
  for (model in c(names(models), "Treat all")) {
    a <- z[z$model == model, ]
    lines(a$threshold, a$net_benefit, col = if (model == "Treat all") "gray50" else colors[model],
          lty = if (model == "Treat all") 2 else 1, lwd = 1.5)
  }
}
invisible(dev.off())
writeLines(c("Probabilities were evaluated as supplied; no fitting, rescaling, smoothing, or patient deletion.",
 "AUC 95% CI: DeLong. Calibration bars: Wilson intervals for grouped outcomes.",
 "DCA threshold range 0.01-0.35 is illustrative; prespecify a clinically justified range.",
 "The figures do not establish independent external validation if the input cohort was outcome-selected.",
 "Synthetic demo output is not a result from the study."), file.path(out, "METHOD.txt"))
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
print(metrics, row.names = FALSE)
