# Clean evaluation helpers added for this code release, not the historical fit.
validate_predictions <- function(d) {
  required <- c("patient_id", "cohort", "outcome", "clinical_probability",
                "radiomics_probability", "fusion_probability")
  if (!all(required %in% names(d))) stop("Missing required columns")
  if (!nrow(d) || anyNA(d[, required])) stop("Empty input or missing values")
  if (any(!nzchar(trimws(as.character(d$patient_id)))) ||
      any(!nzchar(trimws(as.character(d$cohort))))) stop("Empty ID or cohort")
  if (anyDuplicated(d$patient_id)) stop("Patient IDs must be globally unique")
  if (!is.numeric(d$outcome) || !all(d$outcome %in% c(0, 1))) stop("Outcome must be numeric 0/1")
  for (nm in required[4:6]) {
    if (!is.numeric(d[[nm]]) || any(!is.finite(d[[nm]])) ||
        any(d[[nm]] < 0 | d[[nm]] > 1)) stop("Invalid probability: ", nm)
  }
  for (c in unique(d$cohort)) {
    if (length(unique(d$outcome[d$cohort == c])) != 2L) stop("Each cohort must contain both outcomes")
  }
  invisible(d)
}

net_benefit <- function(y, p, thresholds) {
  if (!length(y) || length(y) != length(p) || anyNA(y) || anyNA(p) ||
      !all(y %in% c(0, 1)) || any(!is.finite(p)) || any(p < 0 | p > 1)) stop("Invalid DCA input")
  if (any(!is.finite(thresholds)) || any(thresholds <= 0 | thresholds >= 1)) stop("DCA thresholds must be in (0,1)")
  vapply(thresholds, function(t) {
    decision <- p >= t
    sum(decision & y == 1) / length(y) - sum(decision & y == 0) / length(y) * t / (1 - t)
  }, numeric(1))
}

calibration_bins <- function(y, p, bins = 5L) {
  breaks <- unique(as.numeric(quantile(p, seq(0, 1, length.out = bins + 1L))))
  group <- if (length(breaks) < 2L) rep(1L, length(p)) else
    cut(p, breaks, include.lowest = TRUE, labels = FALSE)
  do.call(rbind, lapply(split(seq_along(y), group), function(idx) {
    n <- length(idx); k <- sum(y[idx]); z <- qnorm(.975)
    centre <- (k + z^2 / 2) / (n + z^2)
    half <- z * sqrt(k * (n - k) / n + z^2 / 4) / (n + z^2)
    data.frame(n = n, events = k, mean_predicted = mean(p[idx]),
               observed = k / n, low95 = max(0, centre - half), high95 = min(1, centre + half))
  }))
}

discrimination <- function(y, p) {
  r <- pROC::roc(y, p, levels = c(0, 1), direction = "<", quiet = TRUE)
  ci <- as.numeric(pROC::ci.auc(r, method = "delong"))
  data.frame(n = length(y), events = sum(y), auc = as.numeric(pROC::auc(r)),
             auc_low95 = ci[1], auc_high95 = ci[3], brier = mean((p - y)^2))
}
