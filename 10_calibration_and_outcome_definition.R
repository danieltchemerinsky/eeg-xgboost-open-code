# ==============================================================================
# CALIBRATION, AND THE CPC 1-2 VERSUS 3-5 OUTCOME DEFINITION
# ==============================================================================
# Two checks on the leave-one-out predictions from 01. Neither refits a model.
#
# 1. Calibration. The paper presents the model's output as a probability, so it
#    reports how close those probabilities are to the observed outcomes, not
#    only how well they rank patients (the AUC). Calibration is assessed for the
#    predicted probability of a poor outcome (CPC 4-5), the positive class in
#    Table 2, with:
#      - the Brier score, and the Brier score of a model that predicts the
#        cohort's poor-outcome rate for every patient, for reference
#      - the scaled Brier score, 1 - Brier / reference (0 = no better than the
#        cohort rate, 1 = perfect)
#      - calibration-in-the-large (intercept with the slope fixed at 1) and the
#        calibration slope, from logistic regression of the outcome on the
#        logit of the predicted probability (ideal: intercept 0, slope 1)
#      - a calibration plot: observed against predicted poor-outcome rates in
#        quintiles of predicted risk, with Wilson intervals, and a loess curve
#    The same Brier score is also given for the full four-class prediction.
#
# 2. Outcome definition. The main analysis dichotomises CPC 1-3 versus 4-5.
#    Guidelines and Johnsen et al. (2022) use CPC 1-2 versus 3-5. The model
#    predicts each CPC category separately, so the alternative definition needs
#    no refit: the probability of a good outcome becomes P(CPC 1) + P(CPC 2)
#    instead of P(CPC 1) + P(CPC 2) + P(CPC 3). Reported: the AUC with a DeLong
#    interval, the Brier score, and poor-outcome prediction at the six Table 2
#    criteria.
#
# Writes calibration_results.csv, calibration_quintiles.csv,
# outcome_definition_cpc12.csv, and figures/FigureS2_Calibration.pdf/.png.
#
# Run after 01, in the same session:
#   source("01_loocv_analysis.R")
#   source("10_calibration_and_outcome_definition.R")
# ==============================================================================

library(pROC)
library(ggplot2)

required <- c("loo_predictions", "actual_outcomes")
missing <- required[!vapply(required, exists, logical(1))]
if (length(missing) > 0)
  stop(paste("Missing objects:", paste(missing, collapse = ", "),
             "\n>>> Run 01_loocv_analysis.R first."))

# Everything below runs in its own scope, so that its working variables
# cannot overwrite objects that other scripts in the same session rely on.
local({

  rule <- function(ch = "=") cat(strrep(ch, 78), "\n")

  wilson_ci <- function(x, n, conf.level = 0.95) {
    if (n == 0) return(c(lower = NA, upper = NA))
    z <- qnorm((1 + conf.level) / 2)
    p <- x / n
    denom <- 1 + z^2 / n
    center <- (p + z^2 / (2 * n)) / denom
    margin <- z * sqrt((p * (1 - p) + z^2 / (4 * n)) / n) / denom
    c(lower = max(0, center - margin), upper = min(1, center + margin))
  }

  # Columns of loo_predictions are classes 0, 1, 2, 3 = CPC 1, 2, 3, 5
  cal_cpc   <- c(1, 2, 3, 5)[actual_outcomes + 1]
  cal_ppoor <- loo_predictions[, 4]                      # P(CPC 5) = P(CPC 4-5)
  cal_poor  <- as.numeric(cal_cpc >= 4)
  n_pat     <- length(cal_poor)

  # ==============================================================================
  # 1. CALIBRATION (poor outcome = CPC 4-5)
  # ==============================================================================
  rule(); cat("1. CALIBRATION - predicted probability of poor outcome (CPC 4-5)\n"); rule()

  brier      <- mean((cal_ppoor - cal_poor)^2)
  brier_ref  <- mean(cal_poor) * (1 - mean(cal_poor))
  brier_scal <- 1 - brier / brier_ref

  # Four-class Brier score: squared error summed over the four classes, averaged
  # over patients (range 0-2)
  onehot <- matrix(0, n_pat, 4); onehot[cbind(seq_len(n_pat), actual_outcomes + 1)] <- 1
  brier_multi <- mean(rowSums((loo_predictions - onehot)^2))
  brier_multi_ref <- mean(rowSums(sweep(-onehot, 2, colMeans(onehot), "+")^2))

  # Logistic recalibration. Probabilities are clipped away from 0 and 1 so the
  # logit is finite.
  eps   <- 1e-6
  lp    <- qlogis(pmin(pmax(cal_ppoor, eps), 1 - eps))
  fit_s <- glm(cal_poor ~ lp, family = binomial)
  fit_i <- glm(cal_poor ~ offset(lp), family = binomial)
  slope     <- unname(coef(fit_s)[2])
  slope_ci  <- unname(suppressMessages(confint.default(fit_s)[2, ]))
  citl      <- unname(coef(fit_i)[1])
  citl_ci   <- unname(suppressMessages(confint.default(fit_i)[1, ]))

  cat(sprintf("  Observed poor-outcome rate        : %.3f  (%d / %d)\n",
              mean(cal_poor), sum(cal_poor), n_pat))
  cat(sprintf("  Mean predicted P(poor)            : %.3f\n", mean(cal_ppoor)))
  cat(sprintf("  Brier score                       : %.3f\n", brier))
  cat(sprintf("  Brier score, cohort rate only     : %.3f\n", brier_ref))
  cat(sprintf("  Scaled Brier score                : %.3f\n", brier_scal))
  cat(sprintf("  Calibration-in-the-large          : %.2f (95%% CI %.2f to %.2f)\n",
              citl, citl_ci[1], citl_ci[2]))
  cat(sprintf("  Calibration slope                 : %.2f (95%% CI %.2f to %.2f)\n",
              slope, slope_ci[1], slope_ci[2]))
  cat(sprintf("  Four-class Brier score            : %.3f (class rates only: %.3f)\n\n",
              brier_multi, brier_multi_ref))
  cat("  Slope < 1: predictions too extreme; slope > 1: too cautious.\n")
  cat("  Intercept > 0: poor outcome under-predicted on average.\n\n")

  # Quintiles of predicted risk
  q_breaks <- unique(quantile(cal_ppoor, probs = seq(0, 1, 0.2), type = 7))
  grp <- cut(cal_ppoor, breaks = q_breaks, include.lowest = TRUE, labels = FALSE)
  quint <- do.call(rbind, lapply(sort(unique(grp)), function(g) {
    idx <- grp == g
    ci  <- wilson_ci(sum(cal_poor[idx]), sum(idx))
    data.frame(Group = g, N = sum(idx), Events = sum(cal_poor[idx]),
               Mean_predicted = mean(cal_ppoor[idx]),
               Observed = mean(cal_poor[idx]),
               Obs_L = unname(ci[1]), Obs_U = unname(ci[2]))
  }))
  cat("  Observed vs predicted, quintiles of predicted risk:\n")
  print(data.frame(Group = quint$Group, N = quint$N, Events = quint$Events,
                   Predicted = sprintf("%.3f", quint$Mean_predicted),
                   Observed  = sprintf("%.3f (%.3f-%.3f)", quint$Observed,
                                       quint$Obs_L, quint$Obs_U)),
        row.names = FALSE)
  cat("\n")

  cal_results <- data.frame(
    Quantity = c("N", "Poor outcomes (CPC 4-5)", "Observed poor-outcome rate",
                 "Mean predicted P(poor)", "Brier score",
                 "Brier score, cohort rate only", "Scaled Brier score",
                 "Calibration-in-the-large", "Calibration-in-the-large, CI lower",
                 "Calibration-in-the-large, CI upper", "Calibration slope",
                 "Calibration slope, CI lower", "Calibration slope, CI upper",
                 "Four-class Brier score", "Four-class Brier score, class rates only"),
    Value = c(n_pat, sum(cal_poor), mean(cal_poor), mean(cal_ppoor), brier,
              brier_ref, brier_scal, citl, citl_ci[1], citl_ci[2], slope,
              slope_ci[1], slope_ci[2], brier_multi, brier_multi_ref))
  cal_results$Value <- round(cal_results$Value, 4)
  write.csv(cal_results, "calibration_results.csv", row.names = FALSE)
  write.csv(transform(quint, Mean_predicted = round(Mean_predicted, 4),
                      Observed = round(Observed, 4), Obs_L = round(Obs_L, 4),
                      Obs_U = round(Obs_U, 4)),
            "calibration_quintiles.csv", row.names = FALSE)

  # Calibration plot
  dir.create("figures", showWarnings = FALSE)
  cal_df <- data.frame(pred = cal_ppoor, obs = cal_poor)
  fig_s2 <- ggplot() +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", colour = "grey55") +
    geom_smooth(data = cal_df, aes(pred, obs), method = "loess", formula = y ~ x,
                span = 0.9, se = FALSE, colour = "grey25", linewidth = 0.6) +
    geom_errorbar(data = quint, aes(x = Mean_predicted, ymin = Obs_L, ymax = Obs_U),
                  width = 0.015, colour = "black") +
    geom_point(data = quint, aes(Mean_predicted, Observed), size = 2.4) +
    geom_rug(data = cal_df[cal_df$obs == 1, ], aes(x = pred), sides = "t",
             length = unit(0.02, "npc"), colour = "grey30") +
    geom_rug(data = cal_df[cal_df$obs == 0, ], aes(x = pred), sides = "b",
             length = unit(0.02, "npc"), colour = "grey30") +
    coord_equal(xlim = c(0, 1), ylim = c(0, 1), expand = TRUE) +
    labs(x = "Predicted probability of poor outcome (CPC 4-5)",
         y = "Observed proportion with poor outcome") +
    annotate("text", x = 0.02, y = 0.97, hjust = 0, vjust = 1, size = 3.2,
             label = sprintf("Brier score %.3f (cohort rate only %.3f)\nCalibration slope %.2f, intercept %.2f",
                             brier, brier_ref, slope, citl)) +
    theme_classic(base_size = 10)
  ggsave(file.path("figures", "FigureS2_Calibration.pdf"), fig_s2,
         width = 110, height = 110, units = "mm")
  ggsave(file.path("figures", "FigureS2_Calibration.png"), fig_s2,
         width = 110, height = 110, units = "mm", dpi = 300)
  cat("  Saved: calibration_results.csv, calibration_quintiles.csv,\n")
  cat("         figures/FigureS2_Calibration.pdf/.png\n\n")

  # ==============================================================================
  # 2. OUTCOME DEFINITION: CPC 1-2 (good) VERSUS CPC 3-5 (poor)
  # ==============================================================================
  rule(); cat("2. OUTCOME DEFINITION - CPC 1-2 versus CPC 3-5\n"); rule()

  prob_good12 <- loo_predictions[, 1] + loo_predictions[, 2]
  true_good12 <- as.numeric(cal_cpc <= 2)
  prob_good13 <- loo_predictions[, 1] + loo_predictions[, 2] + loo_predictions[, 3]
  true_good13 <- as.numeric(cal_cpc <= 3)

  roc12 <- roc(true_good12, prob_good12, quiet = TRUE)
  roc13 <- roc(true_good13, prob_good13, quiet = TRUE)
  ci12  <- ci.auc(roc12, method = "delong")
  ci13  <- ci.auc(roc13, method = "delong")
  brier12 <- mean(((1 - prob_good12) - (1 - true_good12))^2)

  cat(sprintf("  Poor outcomes: CPC 3-5 %d, CPC 4-5 %d (the two definitions differ by the %d CPC 3 patients)\n",
              sum(1 - true_good12), sum(1 - true_good13), sum(cal_cpc == 3)))
  cat(sprintf("  AUC, CPC 1-2 vs 3-5 : %.3f (95%% CI %.3f-%.3f)\n", ci12[2], ci12[1], ci12[3]))
  cat(sprintf("  AUC, CPC 1-3 vs 4-5 : %.3f (95%% CI %.3f-%.3f)  [main analysis]\n",
              ci13[2], ci13[1], ci13[3]))
  cat(sprintf("  Brier score, CPC 1-2 vs 3-5 : %.3f\n\n", brier12))

  criteria <- c(0.10, 0.20, 0.50, 0.75, 0.95, 0.98)
  ops12 <- do.call(rbind, lapply(criteria, function(crit) {
    pred_poor <- as.numeric(prob_good12 < crit)
    truth     <- 1 - true_good12
    TP <- sum(pred_poor == 1 & truth == 1); FN <- sum(pred_poor == 0 & truth == 1)
    TN <- sum(pred_poor == 0 & truth == 0); FP <- sum(pred_poor == 1 & truth == 0)
    s <- wilson_ci(TP, TP + FN); sp <- wilson_ci(TN, TN + FP)
    n <- wilson_ci(TN, TN + FN)
    data.frame(Criterion = crit,
               Sensitivity = round(TP / (TP + FN), 3), Sens_L = round(s[1], 3), Sens_U = round(s[2], 3),
               Specificity = round(TN / (TN + FP), 3), Spec_L = round(sp[1], 3), Spec_U = round(sp[2], 3),
               NPV = round(ifelse(TN + FN > 0, TN / (TN + FN), NA), 3),
               NPV_L = round(n[1], 3), NPV_U = round(n[2], 3),
               TP = TP, TN = TN, FP = FP, FN = FN, row.names = NULL)
  }))
  cat("  Poor-outcome (CPC 3-5) prediction at the Table 2 criteria:\n")
  print(ops12, row.names = FALSE)

  out12 <- rbind(
    data.frame(Quantity = "AUC, CPC 1-2 vs 3-5", Value = round(ci12[2], 3),
               CI_L = round(ci12[1], 3), CI_U = round(ci12[3], 3)),
    data.frame(Quantity = "AUC, CPC 1-3 vs 4-5 (main analysis)", Value = round(ci13[2], 3),
               CI_L = round(ci13[1], 3), CI_U = round(ci13[3], 3)),
    data.frame(Quantity = "Brier score, CPC 1-2 vs 3-5", Value = round(brier12, 4),
               CI_L = NA, CI_U = NA))
  write.csv(out12, "outcome_definition_cpc12.csv", row.names = FALSE)
  write.csv(ops12, "outcome_definition_cpc12_criteria.csv", row.names = FALSE)
  cat("\n  Saved: outcome_definition_cpc12.csv, outcome_definition_cpc12_criteria.csv\n")
  rule()
})
