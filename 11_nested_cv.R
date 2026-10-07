# ==============================================================================
# NESTED CROSS-VALIDATION
# ==============================================================================
# In 01, the hyperparameters are fixed across the leave-one-out folds and were
# chosen beforehand by cross-validation on all 114 patients. The held-out
# patient in each fold therefore contributed to choosing the settings used to
# predict it, which can make the leave-one-out estimate optimistic.
#
# This script repeats the leave-one-out analysis with the tuning moved inside
# each fold. For every held-out patient, the other 113 are split into five
# stratified inner folds, and inner cross-validation chooses:
#   - max_depth from {2, 3, 4}, and
#   - the number of boosting rounds (up to MAX_ROUNDS, by the lowest inner
#     multiclass log loss, with early stopping),
# with every other setting as in 01. A model with the chosen settings is then
# fitted to the 113 patients and predicts the held-out patient. That patient is
# never seen during tuning.
#
# Compared with 01: the good-versus-poor AUC (paired DeLong test against the
# fixed-hyperparameter predictions), the Brier score, and poor-outcome
# prediction at the six Table 2 criteria.
#
# Writes nested_cv_summary.csv, nested_cv_criteria.csv, nested_cv_folds.csv and
# nested_cv_predictions.csv.
#
# Run after 01, in the same session:
#   source("01_loocv_analysis.R")
#   source("11_nested_cv.R")
#
# Slow: 114 outer folds x 3 depths x 5 inner folds. Roughly 10-20 minutes,
# depending on the machine.
# ==============================================================================

library(xgboost)
library(pROC)

required <- c("train_data", "train_label", "params", "nrounds", "loo_predictions")
missing <- required[!vapply(required, exists, logical(1))]
if (length(missing) > 0)
  stop(paste("Missing objects:", paste(missing, collapse = ", "),
             "\n>>> Run 01_loocv_analysis.R first."))

# Everything below runs in its own scope, so that its working variables
# cannot overwrite objects that other scripts in the same session rely on.
local({

  rule <- function(ch = "=") cat(strrep(ch, 78), "\n")

  DEPTH_GRID   <- c(2, 3, 4)
  MAX_ROUNDS   <- 2000
  EARLY_STOP   <- 100
  N_INNER      <- 5
  NESTED_SEED_BASE <- 1042    # different from 01's base, so the streams differ

  wilson_ci <- function(x, n, conf.level = 0.95) {
    if (n == 0) return(c(lower = NA, upper = NA))
    z <- qnorm((1 + conf.level) / 2)
    p <- x / n
    denom <- 1 + z^2 / n
    center <- (p + z^2 / (2 * n)) / denom
    margin <- z * sqrt((p * (1 - p) + z^2 / (4 * n)) / n) / denom
    c(lower = max(0, center - margin), upper = min(1, center + margin))
  }

  # Stratified fold assignment: patients of each class dealt round-robin into the
  # folds after shuffling, so every inner fold has a similar class mix.
  stratified_folds <- function(y, k) {
    fold <- integer(length(y))
    for (cl in unique(y)) {
      idx <- which(y == cl)
      idx <- idx[sample.int(length(idx))]
      fold[idx] <- ((seq_along(idx) - 1 + sample.int(k, 1)) %% k) + 1
    }
    lapply(seq_len(k), function(f) which(fold == f))
  }

  # Lowest mean inner-fold log loss, whatever the xgboost version calls the column
  best_iteration <- function(cv) {
    log <- as.data.frame(cv$evaluation_log)
    col <- grep("^test.*mlogloss.*mean$", names(log), value = TRUE)[1]
    it  <- which.min(log[[col]])
    c(iter = log$iter[it], loss = log[[col]][it])
  }

  X_all <- as.matrix(train_data)
  y_all <- train_label
  n     <- nrow(X_all)

  nested_pred  <- matrix(NA_real_, n, 4, dimnames = list(NULL, paste0("P_class_", 0:3)))
  nested_folds <- data.frame(Fold = seq_len(n), max_depth = NA_integer_,
                             nrounds = NA_integer_, inner_logloss = NA_real_)

  rule(); cat("NESTED CROSS-VALIDATION\n"); rule()
  cat(sprintf("  %d outer folds; inner %d-fold CV over max_depth {%s} and up to %d rounds\n",
              n, N_INNER, paste(DEPTH_GRID, collapse = ", "), MAX_ROUNDS))
  cat("  Progress: ")
  t0 <- Sys.time()

  for (i in seq_len(n)) {
    set.seed(NESTED_SEED_BASE + i)
    tr <- setdiff(seq_len(n), i)
    dtr <- xgb.DMatrix(X_all[tr, , drop = FALSE], label = y_all[tr], missing = NA)
    inner <- stratified_folds(y_all[tr], N_INNER)

    best <- NULL
    for (depth in DEPTH_GRID) {
      p <- params; p$max_depth <- depth
      cv <- xgb.cv(params = p, data = dtr, nrounds = MAX_ROUNDS, folds = inner,
                   early_stopping_rounds = EARLY_STOP, verbose = 0)
      b <- best_iteration(cv)
      if (is.null(best) || b["loss"] < best$loss)
        best <- list(depth = depth, nrounds = unname(b["iter"]), loss = unname(b["loss"]))
    }

    p <- params; p$max_depth <- best$depth
    model <- xgb.train(params = p, data = dtr, nrounds = best$nrounds, verbose = 0)
    dte <- xgb.DMatrix(X_all[i, , drop = FALSE], missing = NA)
    nested_pred[i, ] <- matrix(predict(model, dte), nrow = 1, ncol = 4, byrow = TRUE)
    nested_folds[i, 2:4] <- c(best$depth, best$nrounds, best$loss)

    if (i %% 10 == 0)
      cat(i, "(", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), "min) ")
  }
  cat("\n  Done in", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), "minutes\n\n")

  # ------------------------------------------------------------------------------
  # Compare with the fixed-hyperparameter predictions from 01
  # ------------------------------------------------------------------------------
  cpc        <- c(1, 2, 3, 5)[y_all + 1]
  true_good  <- as.numeric(cpc <= 3)
  pg_fixed   <- rowSums(loo_predictions[, 1:3])
  pg_nested  <- rowSums(nested_pred[, 1:3])

  roc_fixed  <- roc(true_good, pg_fixed,  quiet = TRUE)
  roc_nested <- roc(true_good, pg_nested, quiet = TRUE)
  ci_fixed   <- ci.auc(roc_fixed,  method = "delong")
  ci_nested  <- ci.auc(roc_nested, method = "delong")
  test       <- roc.test(roc_fixed, roc_nested, method = "delong", paired = TRUE)
  brier_fixed  <- mean(((1 - pg_fixed)  - (1 - true_good))^2)
  brier_nested <- mean(((1 - pg_nested) - (1 - true_good))^2)

  cat("  Chosen settings across the outer folds:\n")
  cat("    max_depth:", paste(names(table(nested_folds$max_depth)),
                              table(nested_folds$max_depth), sep = ": ", collapse = ", "), "\n")
  cat(sprintf("    nrounds: median %d (range %d-%d); fixed analysis uses %d\n\n",
              as.integer(median(nested_folds$nrounds)), min(nested_folds$nrounds),
              max(nested_folds$nrounds), nrounds))
  cat(sprintf("  AUC, fixed hyperparameters (01) : %.3f (95%% CI %.3f-%.3f)\n",
              ci_fixed[2], ci_fixed[1], ci_fixed[3]))
  cat(sprintf("  AUC, nested                     : %.3f (95%% CI %.3f-%.3f)\n",
              ci_nested[2], ci_nested[1], ci_nested[3]))
  cat(sprintf("  Paired DeLong test              : p = %.3f\n", test$p.value))
  cat(sprintf("  Brier score (poor outcome)      : fixed %.3f, nested %.3f\n\n",
              brier_fixed, brier_nested))

  criteria <- c(0.10, 0.20, 0.50, 0.75, 0.95, 0.98)
  ops <- do.call(rbind, lapply(criteria, function(crit) {
    pred_poor <- as.numeric(pg_nested < crit)
    truth     <- 1 - true_good
    TP <- sum(pred_poor == 1 & truth == 1); FN <- sum(pred_poor == 0 & truth == 1)
    TN <- sum(pred_poor == 0 & truth == 0); FP <- sum(pred_poor == 1 & truth == 0)
    s <- wilson_ci(TP, TP + FN); sp <- wilson_ci(TN, TN + FP); nv <- wilson_ci(TN, TN + FN)
    data.frame(Criterion = crit,
               Sensitivity = round(TP / (TP + FN), 3), Sens_L = round(s[1], 3), Sens_U = round(s[2], 3),
               Specificity = round(TN / (TN + FP), 3), Spec_L = round(sp[1], 3), Spec_U = round(sp[2], 3),
               NPV = round(ifelse(TN + FN > 0, TN / (TN + FN), NA), 3),
               NPV_L = round(nv[1], 3), NPV_U = round(nv[2], 3),
               TP = TP, TN = TN, FP = FP, FN = FN, row.names = NULL)
  }))
  cat("  Poor-outcome prediction at the Table 2 criteria, nested predictions:\n")
  print(ops, row.names = FALSE)

  summary_tab <- data.frame(
    Quantity = c("AUC fixed", "AUC fixed CI lower", "AUC fixed CI upper",
                 "AUC nested", "AUC nested CI lower", "AUC nested CI upper",
                 "Paired DeLong p", "Brier fixed", "Brier nested",
                 "nrounds median", "nrounds min", "nrounds max"),
    Value = c(ci_fixed[2], ci_fixed[1], ci_fixed[3],
              ci_nested[2], ci_nested[1], ci_nested[3],
              test$p.value, brier_fixed, brier_nested,
              median(nested_folds$nrounds), min(nested_folds$nrounds),
              max(nested_folds$nrounds)))
  summary_tab$Value <- round(summary_tab$Value, 4)
  for (d in DEPTH_GRID)
    summary_tab <- rbind(summary_tab, data.frame(
      Quantity = sprintf("Folds choosing max_depth %d", d),
      Value = sum(nested_folds$max_depth == d)))

  write.csv(summary_tab, "nested_cv_summary.csv", row.names = FALSE)
  write.csv(ops, "nested_cv_criteria.csv", row.names = FALSE)
  write.csv(nested_folds, "nested_cv_folds.csv", row.names = FALSE)
  write.csv(nested_pred, "nested_cv_predictions.csv", row.names = FALSE)
  nested_cv_predictions <<- nested_pred
  cat("\n  Saved: nested_cv_summary.csv, nested_cv_criteria.csv, nested_cv_folds.csv,\n")
  cat("         nested_cv_predictions.csv\n")
  rule()
})
