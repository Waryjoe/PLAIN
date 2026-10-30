#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("Usage: generate_table1.R <input.csv> <output.csv>")
}

invisible(Sys.setlocale("LC_ALL", "zh_CN.UTF-8"))
input_path <- args[[1]]
output_path <- args[[2]]
options(stringsAsFactors = FALSE)

data <- read.csv(
  input_path,
  check.names = FALSE,
  fileEncoding = "UTF-8",
  na.strings = c("", "NA", "N/A", "#VALUE!", "#REF!")
)

group <- as.integer(data[["吻合口"]])
if (!identical(sort(unique(group)), c(0L, 1L))) {
  stop("吻合口 must contain exactly 0 and 1")
}

fmt_num <- function(value) {
  if (is.na(value)) return("")
  formatC(value, format = "f", digits = 2)
}

fmt_continuous <- function(values) {
  values <- suppressWarnings(as.numeric(values))
  values <- values[is.finite(values)]
  if (length(values) == 0) return("")
  q <- unname(quantile(values, c(0.25, 0.5, 0.75), na.rm = TRUE))
  paste0(fmt_num(q[[2]]), " (", fmt_num(q[[1]]), ", ", fmt_num(q[[3]]), ")")
}

fmt_count <- function(count, denominator) {
  if (denominator == 0) return(paste0(count, " (NA)"))
  paste0(count, " (", formatC(100 * count / denominator, format = "f", digits = 1), ")")
}

fmt_p <- function(value) {
  if (is.na(value)) return("")
  if (value < 0.001) return("<0.001")
  formatC(value, format = "f", digits = 3)
}

wilcoxon_z <- function(values, groups) {
  keep <- is.finite(values) & !is.na(groups)
  x <- values[keep]
  g <- groups[keep]
  n1 <- sum(g == 1)
  n0 <- sum(g == 0)
  ranks <- rank(x, ties.method = "average")
  u1 <- sum(ranks[g == 1]) - n1 * (n1 + 1) / 2
  tie_sizes <- as.numeric(table(x))
  n <- n1 + n0
  tie_term <- sum(tie_sizes^3 - tie_sizes) / (n * (n - 1))
  variance <- n1 * n0 * ((n + 1) - tie_term) / 12
  if (variance <= 0) return(c(statistic = NA_real_, p = NA_real_))
  z <- (u1 - n1 * n0 / 2) / sqrt(variance)
  c(statistic = z, p = 2 * pnorm(-abs(z)))
}

categorical_test <- function(values, groups) {
  keep <- !is.na(values) & values != "" & !is.na(groups)
  tab <- table(groups[keep], values[keep])
  if (nrow(tab) < 2 || ncol(tab) < 2) {
    return(list(statistic = NA_real_, p = NA_real_, test = "Not estimable"))
  }
  expected <- suppressWarnings(chisq.test(tab, correct = FALSE)$expected)
  if (any(expected < 5)) {
    if (all(dim(tab) == c(2, 2))) {
      result <- fisher.test(tab)
      return(list(statistic = NA_real_, p = result$p.value, test = "Fisher exact"))
    }
    set.seed(20260729)
    result <- tryCatch(
      fisher.test(tab, workspace = 2e8),
      error = function(e) fisher.test(tab, simulate.p.value = TRUE, B = 100000)
    )
    method <- if (grepl("simulated", result$method, ignore.case = TRUE)) {
      "Fisher-Freeman-Halton (Monte Carlo)"
    } else {
      "Fisher-Freeman-Halton exact"
    }
    return(list(statistic = NA_real_, p = result$p.value, test = method))
  }
  result <- chisq.test(tab, correct = FALSE)
  list(statistic = unname(result$statistic), p = result$p.value, test = "Pearson chi-square")
}

recode_sex <- function(values) {
  x <- trimws(tolower(as.character(values)))
  out <- rep(NA_character_, length(x))
  out[x %in% c("男", "male", "m", "1")] <- "男"
  out[x %in% c("女", "female", "f", "2")] <- "女"
  factor(out, levels = c("女", "男"))
}

recode_family <- function(values) {
  x <- trimws(tolower(as.character(values)))
  out <- rep(NA_character_, length(x))
  out[x %in% c("无", "no", "n", "0")] <- "无"
  out[x %in% c("有", "yes", "y", "1")] <- "有"
  factor(out, levels = c("无", "有"))
}

recode_code <- function(values, levels, labels) {
  x <- trimws(tolower(as.character(values)))
  factor(x, levels = levels, labels = labels)
}

continuous_variables <- c(
  "年龄", "BMI",
  "肿瘤最长直径", "肿瘤次长径", "肿瘤短径",
  "CEA", "CA199", "CA125", "AFP", "CA242", "CA724",
  "红细胞比容Hct", "血红蛋白Hb", "红细胞RBC", "淋巴细胞Lymph",
  "血小板计数PLT", "白细胞WBC",
  "嗜碱性粒细胞百分比", "嗜酸性粒细胞百分比",
  "中性粒细胞百分比", "淋巴细胞百分比", "中性粒细胞Neut",
  "钾K", "钠Na", "钙Ca", "氯Cl",
  "总蛋白TP", "白蛋白ALB", "球蛋白GLO", "白蛋白球蛋白比值",
  "C反应蛋白CRP", "前白蛋白PA", "总胆红素TBIL", "ALT",
  "直接胆红素DBIL", "间接胆红素IBIL", "AST", "碱性磷酸酶ALP",
  "乳酸脱氢酶LDH", "甘油三酯TG", "总胆固醇TC",
  "高密度脂蛋白", "低密度脂蛋白",
  "凝血酶原时间PT", "纤维蛋白原Fbg", "凝血酶时间TT"
)

categorical_variables <- list(
  list(name = "性别", values = recode_sex(data[["性别"]])),
  list(name = "家族史", values = recode_family(data[["家族史"]])),
  list(
    name = "肿瘤部位",
    values = recode_code(
      data[["肿瘤部位2"]],
      c("1", "2", "3", "4"),
      c("胃上部", "胃中部", "胃下部", "全胃")
    )
  ),
  list(
    name = "分化程度",
    values = recode_code(
      data[["分化程度2"]],
      c("1", "2", "3"),
      c("低危", "中危", "高危")
    )
  ),
  list(
    name = "T分期",
    values = recode_code(data[["T1"]], c("1", "2", "3", "4"), c("T1", "T2", "T3", "T4"))
  ),
  list(
    name = "N分期",
    values = recode_code(data[["N1"]], c("0", "1", "2", "3"), c("N0", "N1", "N2", "N3"))
  ),
  list(
    name = "TG2（数据库编码）",
    values = recode_code(data[["TG2"]], c("0", "1", "2"), c("0", "1", "2"))
  )
)

output_rows <- list()
append_row <- function(
  variable, total = "", no_leak = "", leak = "", statistic = "",
  p_value = "", test = "", row_type = "variable", p_numeric = NA_real_
) {
  output_rows[[length(output_rows) + 1]] <<- data.frame(
    Variable = variable,
    Total = total,
    No_leak = no_leak,
    Leak = leak,
    Statistic = statistic,
    P = p_value,
    Test = test,
    Row_type = row_type,
    P_numeric = p_numeric,
    check.names = FALSE
  )
}

categorical_names <- vapply(categorical_variables, function(x) x$name, character(1))
ordered_variables <- c(
  "年龄", "性别", "BMI", "家族史", "肿瘤部位", "分化程度",
  "肿瘤最长直径", "肿瘤次长径", "肿瘤短径", "T分期", "N分期",
  "CEA", "CA199", "CA125", "AFP", "CA242", "CA724",
  "红细胞比容Hct", "血红蛋白Hb", "红细胞RBC", "淋巴细胞Lymph",
  "血小板计数PLT", "白细胞WBC",
  "嗜碱性粒细胞百分比", "嗜酸性粒细胞百分比",
  "中性粒细胞百分比", "淋巴细胞百分比", "中性粒细胞Neut",
  "钾K", "钠Na", "钙Ca", "氯Cl",
  "总蛋白TP", "白蛋白ALB", "球蛋白GLO", "白蛋白球蛋白比值",
  "C反应蛋白CRP", "前白蛋白PA", "总胆红素TBIL", "ALT",
  "直接胆红素DBIL", "间接胆红素IBIL", "AST", "碱性磷酸酶ALP",
  "乳酸脱氢酶LDH", "甘油三酯TG", "TG2（数据库编码）",
  "总胆固醇TC", "高密度脂蛋白", "低密度脂蛋白",
  "凝血酶原时间PT", "纤维蛋白原Fbg", "凝血酶时间TT"
)

for (variable in ordered_variables) {
  if (variable %in% categorical_names) {
    config <- categorical_variables[[match(variable, categorical_names)]]
    values <- config$values
    test <- categorical_test(as.character(values), group)
    statistic <- if (is.na(test$statistic)) "Fisher" else fmt_num(test$statistic)
    append_row(
      paste0(variable, ", n (%)"),
      statistic = statistic,
      p_value = fmt_p(test$p),
      test = test$test,
      row_type = "variable",
      p_numeric = test$p
    )
    for (level in levels(values)) {
      nonmissing_total <- sum(!is.na(values))
      nonmissing_no <- sum(!is.na(values) & group == 0)
      nonmissing_leak <- sum(!is.na(values) & group == 1)
      append_row(
        paste0("  ", level),
        total = fmt_count(sum(values == level, na.rm = TRUE), nonmissing_total),
        no_leak = fmt_count(sum(values == level & group == 0, na.rm = TRUE), nonmissing_no),
        leak = fmt_count(sum(values == level & group == 1, na.rm = TRUE), nonmissing_leak),
        row_type = "level"
      )
    }
    missing_total <- sum(is.na(values))
    if (missing_total > 0) {
      append_row(
        "  Missing, n",
        total = as.character(missing_total),
        no_leak = as.character(sum(is.na(values) & group == 0)),
        leak = as.character(sum(is.na(values) & group == 1)),
        row_type = "missing"
      )
    }
  } else {
    values <- suppressWarnings(as.numeric(data[[variable]]))
    result <- wilcoxon_z(values, group)
    append_row(
      paste0(variable, ", median (Q1, Q3)"),
      total = fmt_continuous(values),
      no_leak = fmt_continuous(values[group == 0]),
      leak = fmt_continuous(values[group == 1]),
      statistic = fmt_num(result[["statistic"]]),
      p_value = fmt_p(result[["p"]]),
      test = "Mann-Whitney U",
      row_type = "variable",
      p_numeric = result[["p"]]
    )
    missing_total <- sum(!is.finite(values))
    if (missing_total > 0) {
      append_row(
        "  Missing, n",
        total = as.character(missing_total),
        no_leak = as.character(sum(!is.finite(values) & group == 0)),
        leak = as.character(sum(!is.finite(values) & group == 1)),
        row_type = "missing"
      )
    }
  }
}

output <- do.call(rbind, output_rows)
write.csv(output, output_path, row.names = FALSE, fileEncoding = "UTF-8")
