# ==============================================================================
# RENDER TABLE 1 FOR PUBLICATION
# ==============================================================================
# One-vs-all classification for each CPC category (N = 114).
#
# Reads the results produced by 04_tables_and_aucs.R and holds no numbers
# of its own, so the rendered table always matches the analysis.
#
# Run after 01 and 04, in the same session or the same directory:
#   source("04_tables_and_aucs.R")   # writes Table1_with_CI.csv
#   source("06_render_table1.R")
# ==============================================================================

INPUT_CSV <- "Table1_with_CI.csv"

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

t1 <- read.csv(INPUT_CSV, stringsAsFactors = FALSE)

cat("Reading:", normalizePath(INPUT_CSV), "\n")
cat("  last modified:", format(file.info(INPUT_CSV)$mtime), "\n")
cat("  rows:", nrow(t1), "\n\n")

required_cols <- c("CPC", "Sensitivity", "Sens_L", "Sens_U",
                   "Specificity", "Spec_L", "Spec_U",
                   "PPV", "PPV_L", "PPV_U", "NPV", "NPV_L", "NPV_U",
                   "N_predicted")
missing_cols <- setdiff(required_cols, names(t1))
if (length(missing_cols))
  stop("Columns missing from ", INPUT_CSV, ": ", paste(missing_cols, collapse = ", "))

# ------------------------------------------------------------------------------
# Format as percentages with confidence intervals
# ------------------------------------------------------------------------------
# sprintf rounds half to even, applied once to the full-precision value. Do not
# substitute a round-half-up helper, and do not round in two steps: 2/31 is
# 6.45%, which rounds to 6, but rounding first to 6.5 and then to an integer
# gives 7.
pct_ci <- function(value, lower, upper) {
  sprintf("%.0f (%.0f–%.0f)", value * 100, lower * 100, upper * 100)
}

table1_data <- data.frame(
  CPC         = paste("CPC", t1$CPC),
  Sensitivity = pct_ci(t1$Sensitivity, t1$Sens_L, t1$Sens_U),
  Specificity = pct_ci(t1$Specificity, t1$Spec_L, t1$Spec_U),
  PPV         = pct_ci(t1$PPV, t1$PPV_L, t1$PPV_U),
  NPV         = pct_ci(t1$NPV, t1$NPV_L, t1$NPV_U),
  N_Predicted = t1$N_predicted,
  stringsAsFactors = FALSE
)

cat("Table 1 as it will be rendered:\n\n")
print(table1_data, row.names = FALSE)
cat("\n")

# ------------------------------------------------------------------------------
# gt - HTML and PNG
# ------------------------------------------------------------------------------
table1_gt <- table1_data %>%
  gt() %>%
  cols_label(
    CPC         = md(""),
    Sensitivity = md("**Sensitivity % (95% CI)**"),
    Specificity = md("**Specificity % (95% CI)**"),
    PPV         = md("**PPV % (95% CI)**"),
    NPV         = md("**NPV % (95% CI)**"),
    N_Predicted = md("**N Predicted**")
  ) %>%
  cols_align(align = "center", columns = everything()) %>%
  tab_style(style = cell_text(weight = "bold"),
            locations = cells_body(columns = CPC)) %>%
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

gtsave(table1_gt, filename = "Table1_publication.html")
cat("Saved: Table1_publication.html\n")

# PNG export needs a headless browser (webshot2 + chromote). Not fatal.
png_ok <- tryCatch({
  gtsave(table1_gt, filename = "Table1_publication.png", vwidth = 1200)
  TRUE
}, error = function(e) { cat("  [PNG export skipped:", conditionMessage(e), "]\n"); FALSE })
if (png_ok) cat("Saved: Table1_publication.png\n")

# ------------------------------------------------------------------------------
# flextable - Word
# ------------------------------------------------------------------------------
table1_flex <- flextable(table1_data) %>%
  set_header_labels(
    CPC = "",
    Sensitivity = "Sensitivity % (95% CI)",
    Specificity = "Specificity % (95% CI)",
    PPV = "PPV % (95% CI)",
    NPV = "NPV % (95% CI)",
    N_Predicted = "N Predicted"
  ) %>%
  bold(part = "header") %>%
  bold(j = 1, part = "body") %>%
  align(align = "center", part = "all") %>%
  font(fontname = "Times New Roman", part = "all") %>%
  fontsize(size = 10, part = "all") %>%
  border_remove() %>%
  hline_top(border = fp_border(color = "black", width = 2), part = "header") %>%
  hline_bottom(border = fp_border(color = "black", width = 1), part = "header") %>%
  hline_bottom(border = fp_border(color = "black", width = 2), part = "body") %>%
  autofit() %>%
  width(j = 1, width = 0.6) %>%
  width(j = 2:5, width = 1.5) %>%
  width(j = 6, width = 0.9)

save_as_docx(table1_flex, path = "Table1_publication.docx")
cat("Saved: Table1_publication.docx\n")

img_ok <- tryCatch({
  save_as_image(table1_flex, path = "Table1_flextable.png", res = 300)
  TRUE
}, error = function(e) { cat("  [image export skipped:", conditionMessage(e), "]\n"); FALSE })
if (img_ok) cat("Saved: Table1_flextable.png\n")

cat("\n============================================\n")
cat("TABLE 1 RENDERED from ", INPUT_CSV, "\n", sep = "")
cat("The .docx holds a real Word table - preferable to pasting an image,\n")
cat("since journals typeset tables and reviewers cannot search a picture.\n")
cat("============================================\n")
