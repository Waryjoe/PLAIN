source("scripts/evaluation_helpers.R")
if (!requireNamespace("pROC", quietly = TRUE)) stop("pROC is required")
expect_error <- function(expr) {
  failed <- inherits(tryCatch(force(expr), error = function(e) e), "error")
  stopifnot(failed)
}
y <- c(0, 0, 1, 1); p <- c(.1, .4, .35, .8)
stopifnot(abs(discrimination(y, p)$auc - .75) < 1e-12)
stopifnot(abs(discrimination(y, rep(.5, 4))$auc - .5) < 1e-12)
stopifnot(abs(net_benefit(y, p, .3) - (2 / 4 - 1 / 4 * .3 / .7)) < 1e-12)
stopifnot(net_benefit(y, rep(0, 4), .3) == 0)
stopifnot(net_benefit(y, rep(1, 4), .5) == 0)
expect_error(net_benefit(y, p, c(0, 1)))
cal <- calibration_bins(y, rep(.5, 4))
stopifnot(nrow(cal) == 1, cal$n == 4, cal$observed == .5,
          cal$low95 < .5, cal$high95 > .5)
d <- data.frame(patient_id = paste0("SYN", 1:4), cohort = "SYNTHETIC", outcome = y,
  clinical_probability = p, radiomics_probability = p, fusion_probability = p)
validate_predictions(d)
for (field in c("clinical_probability", "radiomics_probability", "fusion_probability")) {
  b <- d; b[[field]][1] <- NA; expect_error(validate_predictions(b))
  b[[field]][1] <- 1.1; expect_error(validate_predictions(b))
  b[[field]][1] <- Inf; expect_error(validate_predictions(b))
}
b <- d; b$patient_id[2] <- b$patient_id[1]; expect_error(validate_predictions(b))
b <- d; b$outcome <- 0; expect_error(validate_predictions(b))
b <- d; b$outcome[1] <- 2; expect_error(validate_predictions(b))
b <- d; b$patient_id[1] <- ""; expect_error(validate_predictions(b))
b <- d; b$clinical_probability <- as.character(b$clinical_probability); expect_error(validate_predictions(b))
stopifnot(identical(d$clinical_probability, p))
cat("PASS: AUC/ties, DCA arithmetic, calibration, invalid inputs, and non-mutation.\n")
