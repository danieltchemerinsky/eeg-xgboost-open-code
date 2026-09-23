# ==============================================================================
# PARTIAL DEPENDENCE VALUES - the numbers behind Figure S1
# ==============================================================================
# Figure S1 is a set of partial dependence panels, and the supplementary text
# describes what each one shows ("negative feature 1 values are indicative of a
# good outcome", and so on). This exports the underlying curves and summarises
# each panel numerically, so such a statement can be made from values rather
# than from the appearance of the plot.
#
# The computation matches generate_figure_s1() in 03_figures.R exactly: a model
# fitted to all 114 patients, each feature swept over its 1st-99th percentile in
# 50 steps with all other features left at their observed values, predictions
# averaged over patients, then centred on the feature's own mean.
#
# Run after 01, in the same session:
#   source("01_loocv_analysis.R")
#   source("08_partial_dependence_values.R")
# ==============================================================================

library(xgboost)

required <- c("train_data", "train_label", "params", "nrounds")
missing <- required[!vapply(required, exists, logical(1))]
if (length(missing) > 0)
  stop(paste("Missing objects:", paste(missing, collapse = ", "),
             "\n>>> Run 01_loocv_analysis.R first."))

N_GRID <- 50

# ------------------------------------------------------------------------------
# Model fitted to all patients, and the four features Figure S1 shows
# ------------------------------------------------------------------------------
dtrain_full <- xgb.DMatrix(data = as.matrix(train_data), label = train_label,
                           missing = NA)
set.seed(42)
final_model <- xgb.train(params = params, data = dtrain_full,
                         nrounds = nrounds, verbose = 0)

# Rank by mean gain across the leave-one-out folds, matching Figure 2 and the
# feature selection in generate_all_figures(). A single model fitted to all
# patients gives a different order for ranks 3 and 4 - the two are separated by
# about 3%, so the ranking method decides which comes first. Using one basis
# everywhere keeps the figures and these numbers consistent.
if (exists("all_importance") && nrow(all_importance) > 0) {
  suppressPackageStartupMessages({library(dplyr); library(tidyr)})
  fold_rank <- all_importance %>%
    dplyr::select(Feature, Fold, Gain) %>%
    tidyr::complete(Feature, Fold = 1:max(Fold), fill = list(Gain = 0)) %>%
    dplyr::group_by(Feature) %>%
    dplyr::summarise(Mean = mean(Gain), .groups = "drop") %>%
    dplyr::arrange(desc(Mean))
  top_features <- head(fold_rank$Feature, 4)
  top_gains    <- head(fold_rank$Mean, 4)
  rank_basis   <- "mean gain across the 114 leave-one-out folds"
} else {
  warning("all_importance not found: ranking from a single model instead, which ",
          "may not match Figure S1. Run 01_loocv_analysis.R first.")
  imp <- xgb.importance(model = final_model)
  top_features <- imp$Feature[1:4]
  top_gains    <- imp$Gain[1:4]
  rank_basis   <- "gain in a single model fitted to all patients"
}

cat("Four highest-ranked features, by ", rank_basis, ":\n", sep = "")
for (i in 1:4) cat(sprintf("  Rank %d: %-22s Gain %.4f\n",
                           i, top_features[i], top_gains[i]))
cat(sprintf("\n  Margin between ranks 3 and 4: %.4f (%.1f%% of rank 3) - too\n",
            top_gains[3] - top_gains[4],
            100 * (top_gains[3] - top_gains[4]) / top_gains[3]))
cat("  narrow to treat their order as a finding.\n\n")

# ------------------------------------------------------------------------------
# Partial dependence curves
# ------------------------------------------------------------------------------
pd_all <- data.frame()

for (i in 1:4) {
  raw_name <- top_features[i]
  feat_vals <- train_data[[raw_name]]
  xlims <- quantile(feat_vals, probs = c(0.01, 0.99), na.rm = TRUE)
  grid_vals <- seq(xlims[1], xlims[2], length.out = N_GRID)

  curve <- data.frame()
  for (j in seq_along(grid_vals)) {
    temp <- train_data
    temp[[raw_name]] <- grid_vals[j]
    preds <- matrix(predict(final_model, as.matrix(temp)), ncol = 4, byrow = TRUE)
    curve <- rbind(curve, data.frame(
      x            = grid_vals[j],
      Good         = mean(preds[, 1]),                 # CPC 1
      Intermediate = mean(preds[, 2] + preds[, 3]),    # CPC 2 and 3
      Poor         = mean(preds[, 4])                  # CPC 5
    ))
  }

  # Centre each class on its own mean, as the figure does
  curve$Good         <- curve$Good         - mean(curve$Good)
  curve$Intermediate <- curve$Intermediate - mean(curve$Intermediate)
  curve$Poor         <- curve$Poor         - mean(curve$Poor)

  pd_all <- rbind(pd_all, data.frame(Rank = i, Feature = raw_name, curve))
}

write.csv(pd_all, "partial_dependence_curves.csv", row.names = FALSE)

# ------------------------------------------------------------------------------
# Summary: what does each panel actually separate?
# ------------------------------------------------------------------------------
# For each class: the value at the low end of the range, the value at the high
# end, and the total swing. The class with the largest swing is the one the
# feature moves most; a class whose swing is small is barely affected, however
# suggestive the coloured line looks next to the others.

cat(strrep("=", 78), "\n")
cat("WHAT EACH PANEL SEPARATES\n")
cat(strrep("=", 78), "\n\n")

summary_all <- data.frame()

for (i in 1:4) {
  d <- pd_all[pd_all$Rank == i, ]
  lo <- head(d, 5)     # first 10% of the sweep
  hi <- tail(d, 5)     # last 10%

  s <- data.frame(
    Rank = i, Feature = d$Feature[1],
    Class = c("Good (CPC 1)", "Intermediate (CPC 2-3)", "Poor (CPC 5)"),
    At_low_x  = c(mean(lo$Good), mean(lo$Intermediate), mean(lo$Poor)),
    At_high_x = c(mean(hi$Good), mean(hi$Intermediate), mean(hi$Poor)),
    stringsAsFactors = FALSE
  )
  s$Change <- s$At_high_x - s$At_low_x
  s$Swing  <- c(diff(range(d$Good)), diff(range(d$Intermediate)), diff(range(d$Poor)))
  s$Share_of_total_swing <- s$Swing / sum(s$Swing)

  cat(sprintf("Rank %d: %s\n", i, d$Feature[1]))
  print(data.frame(
    Class     = s$Class,
    low_x     = sprintf("%+.4f", s$At_low_x),
    high_x    = sprintf("%+.4f", s$At_high_x),
    change    = sprintf("%+.4f", s$Change),
    swing     = sprintf("%.4f", s$Swing),
    share     = sprintf("%.0f%%", 100 * s$Share_of_total_swing)
  ), row.names = FALSE)

  dominant <- s$Class[which.max(s$Swing)]
  weakest  <- s$Class[which.min(s$Swing)]
  ratio    <- max(s$Swing) / min(s$Swing)
  cat(sprintf("  -> moves %s most, %s least (%.1fx smaller)\n\n",
              dominant, weakest, ratio))

  summary_all <- rbind(summary_all, s)
}

write.csv(summary_all, "partial_dependence_summary.csv", row.names = FALSE)

cat(strrep("=", 78), "\n")
cat("Saved: partial_dependence_curves.csv, partial_dependence_summary.csv\n")
cat("\nRead the `share` column before describing a panel. A class with a small\n")
cat("share is barely moved by the feature, however far its line sits from the\n")
cat("others - the vertical offset between lines is not the effect, the slope is.\n")
cat(strrep("=", 78), "\n")
