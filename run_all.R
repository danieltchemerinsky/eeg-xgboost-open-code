# ==============================================================================
# RUN THE COMPLETE ANALYSIS
# ==============================================================================
# Reproduces every number, table and figure in the paper, in order, in one
# session.
#
#   1. Put the two data files in  data/   (see data/README.md)
#   2. Run it, either way:
#
#        From R or RStudio, with this folder as the working directory:
#          setwd("path/to/eeg-xgboost-open-code")
#          source("run_all.R")
#
#        From a terminal, in this folder:
#          Rscript run_all.R
#
# Everything produced is written to  output/  . Nothing outside that directory
# is modified, so the run can be repeated safely.
#
# Takes roughly three minutes: the leave-one-out cross-validation fits 114
# models in script 01 and a further 4 x 114 in script 02.
# ==============================================================================

main <- function() {

  REPO_DIR <- getwd()

  if (!file.exists("01_loocv_analysis.R")) {
    cat("\nThis must be run with the repository folder as the working directory.\n")
    cat("Currently: ", REPO_DIR, "\n", sep = "")
    cat('Set it with  setwd("path/to/eeg-xgboost-open-code")  and try again.\n\n')
    return(invisible(FALSE))
  }

  DATA_DIR <- file.path(REPO_DIR, "data")
  OUT_DIR  <- file.path(REPO_DIR, "output")
  on.exit(setwd(REPO_DIR), add = TRUE)   # always restore, even on error

  # Print each warning as it happens, so it is attributed to the step that
  # caused it. The default defers them all to the end, where a bare count of
  # warnings says nothing about where they came from.
  old_warn <- getOption("warn")
  options(warn = 1)
  on.exit(options(warn = old_warn), add = TRUE)

  DATA_FILES <- c("outputminus.RData",
                  "AllPatientsEEGandETCAbsentObs5.RData")

  REQUIRED_PACKAGES <- c(
    "xgboost", "pROC", "caret", "dplyr", "tidyr", "Matrix",            # analysis
    "ggplot2", "patchwork", "scales",                                  # figures
    "gt", "flextable", "officer")                                      # tables

  rule <- function(ch = "=") cat(strrep(ch, 78), "\n")

  rule(); cat("EEG REACTIVITY / XGBOOST - FULL ANALYSIS\n"); rule(); cat("\n")

  # ----------------------------------------------------------------------------
  # Pre-flight
  # ----------------------------------------------------------------------------
  problems <- character()

  missing_data <- DATA_FILES[!file.exists(file.path(DATA_DIR, DATA_FILES))]
  if (length(missing_data))
    problems <- c(problems, paste0(
      "Data file not found in data/: ", paste(missing_data, collapse = ", "),
      "\n      See data/README.md. The data are not distributed with this code;",
      "\n      see the data availability statement in the paper."))

  have <- vapply(REQUIRED_PACKAGES, requireNamespace, logical(1), quietly = TRUE)
  if (any(!have))
    problems <- c(problems, paste0(
      "Missing R package(s): ", paste(REQUIRED_PACKAGES[!have], collapse = ", "),
      "\n      install.packages(c(",
      paste(sprintf('"%s"', REQUIRED_PACKAGES[!have]), collapse = ", "), "))"))

  if (length(problems)) {
    cat("Cannot start:\n\n")
    for (p in problems) cat("  - ", p, "\n", sep = "")
    cat("\n")
    return(invisible(FALSE))
  }

  cat("Environment\n")
  cat("  ", R.version.string, "\n", sep = "")
  for (p in c("xgboost", "pROC", "ggplot2"))
    cat(sprintf("   %-10s %s\n", p, as.character(packageVersion(p))))
  cat("\n")

  # ----------------------------------------------------------------------------
  # Work in output/, with the data alongside so relative paths resolve
  # ----------------------------------------------------------------------------
  dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
  for (f in DATA_FILES) {
    dest <- file.path(OUT_DIR, f)
    if (!file.exists(dest)) file.copy(file.path(DATA_DIR, f), dest)
  }
  setwd(OUT_DIR)

  run <- function(label, expr) {
    cat(sprintf("%-46s", paste0("  ", label, " ... "))); flush.console()
    t0 <- Sys.time()
    ok <- tryCatch({ force(expr); TRUE },
                   error = function(e) {
                     cat("FAILED\n      ", conditionMessage(e), "\n", sep = "")
                     FALSE })
    if (ok) cat(sprintf("ok  (%4.1f s)\n",
                        as.numeric(difftime(Sys.time(), t0, units = "secs"))))
    ok
  }

  src <- function(f) source(file.path(REPO_DIR, f), local = FALSE)

  rule("-"); cat("Running\n"); rule("-")

  steps_ok <- c(
    run("01  leave-one-out cross-validation", src("01_loocv_analysis.R")),
    run("02  reactivity vs all predictors",   src("02_roc_reactivity_vs_all_predictors.R")),
    run("04  tables and AUCs",                src("04_tables_and_aucs.R")),
    run("05  reported statistics",            src("05_reported_statistics.R")),
    run("08  partial dependence values",      src("08_partial_dependence_values.R")),
    run("03  figure functions",               src("03_figures.R")),
    run("    Figures 2, 3 and S1",
        generate_all_figures(output_dir = "figures")),
    run("    Figure 4", {
          fig4 <- generate_figure4_clean(roc_114_react, roc_114_all,
                                         roc_103_react, roc_103_all,
                                         sep_sens, sep_spec)
          ggplot2::ggsave(file.path("figures", "Figure4_ROC_TwoPanel.pdf"), fig4,
                          width = 180, height = 90, units = "mm")
          ggplot2::ggsave(file.path("figures", "Figure4_ROC_TwoPanel.png"), fig4,
                          width = 180, height = 90, units = "mm", dpi = 300)
        }),
    run("09  Figure 1 flowchart",             src("09_figure1_flowchart.R")),
    run("06  render Table 1",                 src("06_render_table1.R")),
    run("07  render Table 2",                 src("07_render_table2.R"))
  )

  cat("\n")

  # ----------------------------------------------------------------------------
  # What was produced
  # ----------------------------------------------------------------------------
  rule("-"); cat("Output\n"); rule("-")

  expected <- c(
    "Table1_with_CI.csv", "Table2_with_CI.csv",
    "AUC_with_CI.csv", "Figure4_AUC_with_CI.csv",
    "reported_statistics.csv", "partial_dependence_summary.csv",
    "Table1_publication.docx", "Table2_publication.docx",
    "Figure1_Flowchart.pdf",
    "figures/Figure2_FeatureImportance.pdf",
    "figures/Figure3_ROC_OneVsAll.pdf",
    "figures/Figure4_ROC_TwoPanel.pdf",
    "figures/FigureS1_PartialDependence.pdf")

  present <- file.exists(expected)
  for (i in seq_along(expected))
    cat(sprintf("  %-44s %s\n", expected[i], if (present[i]) "present" else "MISSING"))

  writeLines(capture.output(sessionInfo()), "sessionInfo.txt")
  cat("\n  sessionInfo.txt written alongside the results.\n\n")

  ok <- all(steps_ok) && all(present)
  rule()
  if (ok) cat("COMPLETE. Every script ran and every expected output is present.\n")
  else    cat("FINISHED WITH PROBLEMS. See the log above.\n")
  cat("Results are in: ", OUT_DIR, "\n", sep = "")
  rule()

  invisible(ok)
}

ok <- main()

# Non-zero exit for Rscript, so this works as a check in a shell or on CI.
# Interactive sessions are left alone.
if (!interactive() && !isTRUE(ok)) quit(save = "no", status = 1)
