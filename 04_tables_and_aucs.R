# ==============================================================================
# TABLES 1 AND 2, AND THE AUCs
# ==============================================================================
# Builds Table 1 (one-vs-all classification at a 50% threshold) and Table 2
# (poor-outcome prediction at six probability criteria) from the leave-one-out
# predictions, with Wilson score intervals, and computes the good-versus-poor
# AUC and the one-vs-all AUC for each CPC category with DeLong intervals.
#
# Exports Table1_with_CI.csv, Table2_with_CI.csv and AUC_with_CI.csv. These are
# the source for the rendered tables: 06_render_table1.R and 07_render_table2.R
# read them and hold no numbers of their own.
#
# Run after 01, in the same session:
#   source("01_loocv_analysis.R")
#   source("04_tables_and_aucs.R")
# ==============================================================================

library(pROC)

# ------------------------------------------------------------------------------
# CHECK: Do the required objects exist from 01_loocv_analysis.R?
# ------------------------------------------------------------------------------
required_objects <- c("loo_predictions", "actual_outcomes", "full", "train_label")
missing <- required_objects[!sapply(required_objects, exists)]
if (length(missing) > 0) {
  stop(paste("Missing objects:", paste(missing, collapse = ", "),
             "\n>>> Run 01_loocv_analysis.R first."))
}

cat("============================================================\n")
cat("TABLES 1 AND 2, AND THE AUCs\n")
cat("============================================================\n")
cat("N patients:", nrow(full), "\n\n")

# ------------------------------------------------------------------------------
# SETUP: Derive key variables
# ------------------------------------------------------------------------------
actual_cpc <- c(1, 2, 3, 5)[actual_outcomes + 1]
prob_good  <- loo_predictions[, 1] + loo_predictions[, 2] + loo_predictions[, 3]
true_good  <- as.numeric(actual_cpc <= 3)

# Wilson score CI (self-contained, no external package needed)
wilson_ci <- function(x, n, conf.level = 0.95) {
  if (n == 0) return(c(lower = NA, upper = NA))
  z <- qnorm((1 + conf.level) / 2)
  p <- x / n
  denom <- 1 + z^2 / n
  center <- (p + z^2 / (2 * n)) / denom
  margin <- z * sqrt((p * (1 - p) + z^2 / (4 * n)) / n) / denom
  c(lower = max(0, center - margin), upper = min(1, center + margin))
}

# ==============================================================================
# TABLE 1: One-vs-All at 50% threshold
# ==============================================================================
cat("============================================================\n")
cat("TABLE 1: One-vs-All Classification (50% threshold)\n")
cat("============================================================\n\n")

class_indices <- c(1, 2, 4)   # columns in loo_predictions
class_labels  <- c(1, 2, 5)   # CPC labels

table1_fresh <- do.call(rbind, lapply(seq_along(class_indices), function(i) {
  idx <- class_indices[i]
  lbl <- class_labels[i]

  pred <- as.numeric(loo_predictions[, idx] >= 0.5)
  true <- as.numeric(actual_outcomes == (idx - 1))

  TP <- sum(pred == 1 & true == 1)
  TN <- sum(pred == 0 & true == 0)
  FP <- sum(pred == 1 & true == 0)
  FN <- sum(pred == 0 & true == 1)

  sens <- TP / (TP + FN)
  spec <- TN / (TN + FP)
  ppv  <- TP / (TP + FP)
  npv  <- TN / (TN + FN)

  data.frame(
    CPC          = lbl,
    Sensitivity  = round(sens, 3),
    Sens_L       = round(wilson_ci(TP, TP + FN)[1], 3),
    Sens_U       = round(wilson_ci(TP, TP + FN)[2], 3),
    Specificity  = round(spec, 3),
    Spec_L       = round(wilson_ci(TN, TN + FP)[1], 3),
    Spec_U       = round(wilson_ci(TN, TN + FP)[2], 3),
    PPV          = round(ppv, 3),
    PPV_L        = round(wilson_ci(TP, TP + FP)[1], 3),
    PPV_U        = round(wilson_ci(TP, TP + FP)[2], 3),
    NPV          = round(npv, 3),
    NPV_L        = round(wilson_ci(TN, TN + FN)[1], 3),
    NPV_U        = round(wilson_ci(TN, TN + FN)[2], 3),
    N_predicted  = sum(pred),
    TP = TP, TN = TN, FP = FP, FN = FN
  )
}))

# Print raw values
cat("Raw values:\n")
print(table1_fresh, row.names = FALSE)
cat("\n")

# Print formatted for paper
cat("Formatted for paper:\n")
for (i in 1:nrow(table1_fresh)) {
  r <- table1_fresh[i, ]
  cat(sprintf("  CPC %d: Sens = %d (%d-%d), Spec = %d (%d-%d), PPV = %d (%d-%d), NPV = %d (%d-%d), N = %d\n",
              r$CPC,
              round(r$Sensitivity*100), round(r$Sens_L*100), round(r$Sens_U*100),
              round(r$Specificity*100), round(r$Spec_L*100), round(r$Spec_U*100),
              round(r$PPV*100), round(r$PPV_L*100), round(r$PPV_U*100),
              round(r$NPV*100), round(r$NPV_L*100), round(r$NPV_U*100),
              r$N_predicted))
}

# ==============================================================================
# TABLE 2: Binary (good vs poor) at different criteria
# ==============================================================================
cat("\n============================================================\n")
cat("TABLE 2: Poor outcome prediction at different criteria\n")
cat("============================================================\n\n")

# The six criteria reported in Table 2.
criteria <- c(0.1, 0.2, 0.5, 0.75, 0.95, 0.98)

table2_fresh <- do.call(rbind, lapply(criteria, function(cut) {
  # Predict POOR when prob_good < criterion
  pred_poor <- as.numeric(prob_good < cut)
  true_poor <- as.numeric(true_good == 0)

  TP <- sum(pred_poor == 1 & true_poor == 1)
  TN <- sum(pred_poor == 0 & true_poor == 0)
  FP <- sum(pred_poor == 1 & true_poor == 0)
  FN <- sum(pred_poor == 0 & true_poor == 1)

  sens <- TP / (TP + FN)
  spec <- TN / (TN + FP)
  ppv  <- ifelse(TP + FP > 0, TP / (TP + FP), NA)
  npv  <- ifelse(TN + FN > 0, TN / (TN + FN), NA)

  data.frame(
    Criterion    = cut,
    Sensitivity  = round(sens, 3),
    Sens_L       = round(wilson_ci(TP, TP + FN)[1], 3),
    Sens_U       = round(wilson_ci(TP, TP + FN)[2], 3),
    Specificity  = round(spec, 3),
    Spec_L       = round(wilson_ci(TN, TN + FP)[1], 3),
    Spec_U       = round(wilson_ci(TN, TN + FP)[2], 3),
    PPV          = round(ppv, 3),
    PPV_L        = round(wilson_ci(TP, TP + FP)[1], 3),
    PPV_U        = round(wilson_ci(TP, TP + FP)[2], 3),
    NPV          = round(npv, 3),
    NPV_L        = round(wilson_ci(TN, TN + FN)[1], 3),
    NPV_U        = round(wilson_ci(TN, TN + FN)[2], 3),
    TP = TP, TN = TN, FP = FP, FN = FN
  )
}))

cat("Raw values:\n")
print(table2_fresh, row.names = FALSE)
cat("\n")

cat("Formatted for paper:\n")
for (i in 1:nrow(table2_fresh)) {
  r <- table2_fresh[i, ]
  cat(sprintf("  Crit %.2f: Sens = %d (%d-%d), Spec = %d (%d-%d), PPV = %d (%d-%d), NPV = %d (%d-%d) | TP=%d TN=%d FP=%d FN=%d\n",
              r$Criterion,
              round(r$Sensitivity*100), round(r$Sens_L*100), round(r$Sens_U*100),
              round(r$Specificity*100), round(r$Spec_L*100), round(r$Spec_U*100),
              round(r$PPV*100), round(r$PPV_L*100), round(r$PPV_U*100),
              round(r$NPV*100), round(r$NPV_L*100), round(r$NPV_U*100),
              r$TP, r$TN, r$FP, r$FN))
}

# ==============================================================================
# AUC with DeLong CI
# ==============================================================================
cat("\n============================================================\n")
cat("AUC: Good vs Poor Outcome\n")
cat("============================================================\n\n")

roc_obj   <- roc(true_good, prob_good, quiet = TRUE)
auc_val   <- auc(roc_obj)
auc_ci    <- ci.auc(roc_obj, method = "delong")

cat(sprintf("AUC = %.3f (95%% CI: %.3f - %.3f)  [DeLong method]\n\n",
            auc_val, auc_ci[1], auc_ci[3]))

# The primary AUC is quoted in the abstract, in Results 3.4 and in the
# Discussion, so it needs a saved artifact rather than living only in memory.
auc_export <- data.frame(
  Model = "Good (CPC 1-3) vs poor (CPC 4-5), reactivity only, N = 114",
  AUC   = round(as.numeric(auc_val), 3),
  CI_L  = round(auc_ci[1], 3),
  CI_U  = round(auc_ci[3], 3),
  CI_method = "DeLong",
  stringsAsFactors = FALSE
)

# ==============================================================================
# One-vs-All AUC per CPC category
# ==============================================================================
cat("One-vs-All AUC per CPC category:\n")
ova_export <- data.frame()
for (class_idx in 1:4) {
  true_class <- as.numeric(actual_outcomes == (class_idx - 1))
  pred_prob  <- loo_predictions[, class_idx]
  if (sum(true_class) > 0 && sum(true_class) < length(true_class)) {
    a <- auc(true_class, pred_prob, quiet = TRUE)
    ci_a <- ci.auc(roc(true_class, pred_prob, quiet = TRUE), method = "delong")
    cpc_label <- c(1, 2, 3, 5)[class_idx]
    cat(sprintf("  CPC %d vs All: %.2f (%.2f-%.2f)\n", cpc_label, a, ci_a[1], ci_a[3]))

    # These are the values printed in the Figure 3 legend and quoted in the
    # Results as "an AUC of 0.59 to 0.79".
    ova_export <- rbind(ova_export, data.frame(
      Model = sprintf("CPC %d vs all", cpc_label),
      AUC   = round(as.numeric(a), 3),
      CI_L  = round(ci_a[1], 3),
      CI_U  = round(ci_a[3], 3),
      CI_method = "DeLong",
      stringsAsFactors = FALSE
    ))
  }
}

write.csv(rbind(auc_export, ova_export), "AUC_with_CI.csv", row.names = FALSE)
cat("\nSaved: AUC_with_CI.csv\n")

# Multiclass AUC (Hand & Till)
# The Hand & Till multiclass AUC is computed in 05_reported_statistics.R,
# together with the other values quoted in the Results text.

# ==============================================================================
# EXPORT
# ==============================================================================
cat("\n============================================================\n")
cat("EXPORTING CSVs\n")
cat("============================================================\n\n")

write.csv(table1_fresh[, c("CPC","Sensitivity","Sens_L","Sens_U",
                            "Specificity","Spec_L","Spec_U",
                            "PPV","PPV_L","PPV_U",
                            "NPV","NPV_L","NPV_U","N_predicted")],
          "Table1_with_CI.csv", row.names = FALSE)

write.csv(table2_fresh[, c("Criterion","Sensitivity","Sens_L","Sens_U",
                            "Specificity","Spec_L","Spec_U",
                            "PPV","PPV_L","PPV_U",
                            "NPV","NPV_L","NPV_U","TP","TN","FP","FN")],
          "Table2_with_CI.csv", row.names = FALSE)

cat("Saved: Table1_with_CI.csv\n")
cat("Saved: Table2_with_CI.csv\n")
cat("Saved: AUC_with_CI.csv\n")
cat("\nThese CSVs are the authoritative source for Tables 1 and 2.\n")
cat("06_render_table1.R and 07_render_table2.R read them directly, so the\n")
cat("rendered tables always match the analysis.\n")
cat("============================================================\n")
