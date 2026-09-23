# ==============================================================================
# RENDER TABLE 2 FOR PUBLICATION
# ==============================================================================
# Classification performance at different probability thresholds, where the
# threshold is applied to the cumulative predicted probability of a good outcome
# (CPC 1-3) and poor outcome (CPC 4-5) is the positive class.
#
# Reads the results produced by 04_tables_and_aucs.R and holds no numbers
# of its own, so the rendered table always matches the analysis.
#
# Run after 01 and 04, in the same session or the same directory:
#   source("04_tables_and_aucs.R")   # writes Table2_with_CI.csv
#   source("07_render_table2.R")
# ==============================================================================

INPUT_CSV <- "Table2_with_CI.csv"

for (pkg in c("gt", "flextable", "officer", "dplyr")) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Missing package '", pkg, "'.\n>>> install.packages(\"", pkg, "\")")
  }
}
library(gt); library(flextable); library(officer); library(dplyr)

# ------------------------------------------------------------------------------
# Read the results
# ------------------------------------------------------------------------------
if (!file.exists(INPUT_CSV)) {
  stop(INPUT_CSV, " not found in ", getwd(),
       "\n>>> Run 01_loocv_analysis.R then 04_tables_and_aucs.R first.")
}

t2 <- read.csv(INPUT_CSV, stringsAsFactors = FALSE)

cat("Reading:", normalizePath(INPUT_CSV), "\n")
cat("  last modified:", format(file.info(INPUT_CSV)$mtime), "\n")
cat("  criteria:", paste(sprintf("%.2f", t2$Criterion), collapse = ", "), "\n\n")

required_cols <- c("Criterion", "Sensitivity", "Sens_L", "Sens_U",
                   "Specificity", "Spec_L", "Spec_U",
                   "PPV", "PPV_L", "PPV_U", "NPV", "NPV_L", "NPV_U",
                   "TP", "TN", "FP", "FN")
missing_cols <- setdiff(required_cols, names(t2))
if (length(missing_cols))
  stop("Columns missing from ", INPUT_CSV, ": ", paste(missing_cols, collapse = ", "))

# Table 2 reports six criteria. Warn rather than silently render a shorter
# table if the criteria vector upstream no longer covers all of them.
if (!all(c(0.10, 0.20, 0.50, 0.75, 0.95, 0.98) %in% round(t2$Criterion, 2))) {
  warning("The six criteria reported in Table 2 are not all present in ",
          INPUT_CSV, ". Check the `criteria` vector in 04_tables_and_aucs.R.")
}

# ------------------------------------------------------------------------------
# Format as percentages with confidence intervals
# ------------------------------------------------------------------------------
# sprintf rounds half to even, applied once to the full-precision value. Do not
# substitute a round-half-up helper, and do not round in two steps: the
# sensitivity at criterion 0.10 is 2/31 = 6.45%, which rounds to 6, but rounding
# first to 6.5 and then to an integer gives 7.
pct_ci <- function(value, lower, upper) {
  sprintf("%.0f (%.0f–%.0f)", value * 100, lower * 100, upper * 100)
}

table2_data <- data.frame(
  Criterion   = sprintf("%.2f", t2$Criterion),
  Sensitivity = pct_ci(t2$Sensitivity, t2$Sens_L, t2$Sens_U),
  Specificity = pct_ci(t2$Specificity, t2$Spec_L, t2$Spec_U),
  PPV         = pct_ci(t2$PPV, t2$PPV_L, t2$PPV_U),
  NPV         = pct_ci(t2$NPV, t2$NPV_L, t2$NPV_U),
  TP = t2$TP, TN = t2$TN, FP = t2$FP, FN = t2$FN,
  stringsAsFactors = FALSE
)

cat("Table 2 as it will be rendered:\n\n")
print(table2_data, row.names = FALSE)
cat("\n")

# ------------------------------------------------------------------------------
# gt - HTML and PNG
# ------------------------------------------------------------------------------
table2_gt <- table2_data %>%
  gt() %>%
  cols_label(
    Criterion   = md("**Criterion**"),
    Sensitivity = md("**Sensitivity % (95% CI)**"),
    Specificity = md("**Specificity % (95% CI)**"),
    PPV         = md("**PPV % (95% CI)**"),
    NPV         = md("**NPV % (95% CI)**"),
    TP = md("**TP**"), TN = md("**TN**"),
    FP = md("**FP**"), FN = md("**FN**")
  ) %>%
  cols_align(align = "center", columns = everything()) %>%
  tab_options(
    table.font.size = px(12),
    table.font.names = "Helvetica",
    column_labels.font.weight = "bold",
    table.border.top.style = "solid",
    table.border.top.width = px(2),
    table.border.top.color = "black",
    table.border.bottom.style = "solid",
    table.border.bottom.width = px(2),
    table.border.bottom.color = "black",
    column_labels.border.bottom.style = "solid",
    column_labels.border.bottom.width = px(1),
    column_labels.border.bottom.color = "black",
    table_body.hlines.style = "none",
    table.width = pct(100)
  )

gtsave(table2_gt, filename = "Table2_publication.html")
cat("Saved: Table2_publication.html\n")

png_ok <- tryCatch({
  gtsave(table2_gt, filename = "Table2_publication.png", vwidth = 1400)
  TRUE
}, error = function(e) { cat("  [PNG export skipped:", conditionMessage(e), "]\n"); FALSE })
if (png_ok) cat("Saved: Table2_publication.png\n")

# ------------------------------------------------------------------------------
# flextable - Word
# ------------------------------------------------------------------------------
table2_flex <- flextable(table2_data) %>%
  set_header_labels(
    Criterion = "Criterion",
    Sensitivity = "Sensitivity % (95% CI)",
    Specificity = "Specificity % (95% CI)",
    PPV = "PPV % (95% CI)",
    NPV = "NPV % (95% CI)",
    TP = "TP", TN = "TN", FP = "FP", FN = "FN"
  ) %>%
  bold(part = "header") %>%
  align(align = "center", part = "all") %>%
  font(fontname = "Times New Roman", part = "all") %>%
  fontsize(size = 10, part = "all") %>%
  border_remove() %>%
  hline_top(border = fp_border(color = "black", width = 2), part = "header") %>%
  hline_bottom(border = fp_border(color = "black", width = 1), part = "header") %>%
  hline_bottom(border = fp_border(color = "black", width = 2), part = "body") %>%
  autofit() %>%
  width(j = 1, width = 0.75) %>%
  width(j = 2:5, width = 1.4) %>%
  width(j = 6:9, width = 0.45)

save_as_docx(table2_flex, path = "Table2_publication.docx")
cat("Saved: Table2_publication.docx\n")

img_ok <- tryCatch({
  save_as_image(table2_flex, path = "Table2_flextable.png", res = 300)
  TRUE
}, error = function(e) { cat("  [image export skipped:", conditionMessage(e), "]\n"); FALSE })
if (img_ok) cat("Saved: Table2_flextable.png\n")

cat("\n============================================\n")
cat("TABLE 2 RENDERED from ", INPUT_CSV, "\n", sep = "")
cat("The .docx holds a real Word table - preferable to pasting an image,\n")
cat("since journals typeset tables and reviewers cannot search a picture.\n")
cat("============================================\n")
