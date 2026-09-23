# ==============================================================================
# LEAVE-ONE-OUT CROSS-VALIDATION
# ==============================================================================
# Fits one model per patient, each trained on the other 113, giving an
# out-of-sample predicted class probability for every patient in the cohort.
# Hyperparameters are fixed across folds rather than tuned within each fold.
#
# Leaves behind the objects the later scripts use: loo_predictions,
# actual_outcomes, train_data, train_label, params, nrounds, full and
# all_importance.
# ==============================================================================

library(xgboost)
library(Matrix)
library(dplyr)
library(pROC)
library(caret)

# If data not loaded, load it
if(!exists("full")) {
  # Load into a private environment and take only `combined`. These .RData files
  # are saved workspaces and carry other objects; loading them into the global
  # environment lets those overwrite objects created by other scripts.
  .data_env <- new.env()
  load("outputminus.RData", envir = .data_env)
  combined <- .data_env$combined
  rm(.data_env)

  absentobs <- rowSums(is.na(combined))
  full <- combined[absentobs <= 5, ]
  train_data <- subset(full, select = -c(Outcome, Patient))
  
  # Proper label mapping
  outcome_to_class <- c("1" = 0, "2" = 1, "3" = 2, "5" = 3)
  train_label <- outcome_to_class[as.character(full$Outcome)]
}

cat("==============================================================================\n")
cat("LEAVE-ONE-OUT CROSS-VALIDATION\n")
cat("==============================================================================\n")
cat("Sample size:", nrow(full), "\n")
cat("This will train", nrow(full), "models (one per patient)\n")
cat("Estimated time: 2-5 minutes depending on your machine\n")
cat("==============================================================================\n\n")

# ==============================================================================
# 1. DEFINE FIXED HYPERPARAMETERS
# ==============================================================================
# These are reasonable defaults that approximate the original paper's tuning
# We use FIXED parameters across all folds for fair comparison

params <- list(
  objective = "multi:softprob",
  eval_metric = "mlogloss",
  num_class = 4,           # CPC 1, 2, 3, 5 → classes 0, 1, 2, 3
  tree_method = "hist",    # Modern method
  eta = 0.01,              # Low learning rate
  max_depth = 3,           # Conservative depth
  min_child_weight = 3,    
  gamma = 0.1,             
  subsample = 0.8,         
  colsample_bytree = 0.8,
  lambda = 1,              
  alpha = 0
)

# Number of boosting rounds (determined by CV on full data)
# You can adjust this, but 200-300 is reasonable for eta=0.01
nrounds <- 750

cat("Hyperparameters:\n")
print(params)
cat("\nNumber of rounds:", nrounds, "\n\n")

# ==============================================================================
# 2. PREPARE STORAGE FOR LOO PREDICTIONS
# ==============================================================================

# Store predictions: each row is one patient, columns are class probabilities
loo_predictions <- matrix(NA, nrow = nrow(full), ncol = 4)
colnames(loo_predictions) <- paste0("P_class_", 0:3)

# Store predicted class (argmax) for each patient
loo_predicted_class <- numeric(nrow(full))

# Store actual outcomes
actual_outcomes <- train_label

# Store per-fold feature importance. Figure 2 shows the spread of each feature's
# Gain across folds, which requires the importance from every fit rather than
# from a single model trained on all patients.
all_importance <- data.frame()

# Progress tracking
start_time <- Sys.time()

# ==============================================================================
# 3. RUN LEAVE-ONE-OUT CROSS-VALIDATION
# ==============================================================================

cat("Starting LOO-CV...\n")
cat("Progress: ")

# Seeding is per fold, not once before the loop.
#
# xgboost's subsample and colsample_bytree draw from R's RNG stream. With a
# single seed before the loop, fold i depends on every draw made before it, so
# adding any line inside the loop that touches the RNG silently changes all 114
# folds and every number downstream. Seeding each fold from its own index makes
# each fit reproducible on its own and independent of everything around it.
LOO_SEED_BASE <- 42

for(i in 1:nrow(full)) {

  set.seed(LOO_SEED_BASE + i)
  
  # Progress indicator
  if(i %% 10 == 0) {
    elapsed <- difftime(Sys.time(), start_time, units = "mins")
    cat(i, "(", round(elapsed, 1), "min) ")
  }
  
  # Create train/test split (leave out patient i)
  train_idx <- setdiff(1:nrow(full), i)
  test_idx <- i
  
  # Prepare training data
  X_train <- as.matrix(train_data[train_idx, ])
  y_train <- train_label[train_idx]
  
  # Prepare test data (single patient)
  X_test <- as.matrix(train_data[test_idx, , drop = FALSE])
  y_test <- train_label[test_idx]
  
  # Create DMatrix objects
  dtrain <- xgb.DMatrix(data = X_train, label = y_train, missing = NA)
  dtest <- xgb.DMatrix(data = X_test, label = y_test, missing = NA)
  
  # Train model on this fold
  model <- xgb.train(
    params = params,
    data = dtrain,
    nrounds = nrounds,
    verbose = 0  # Suppress output for clean progress bar
  )
  
  # Predict on held-out patient
  pred_probs_vec <- predict(model, dtest)
  
  # Reshape predictions (modern XGBoost returns flat vector)
  pred_probs <- matrix(pred_probs_vec, nrow = 1, ncol = 4, byrow = TRUE)
  
  # Store predictions
  loo_predictions[i, ] <- pred_probs
  loo_predicted_class[i] <- which.max(pred_probs) - 1  # Convert to 0-indexed class

  # Feature importance for this fold. xgb.importance() draws no random numbers,
  # so collecting it leaves the predictions above bit-identical.
  imp_i <- xgb.importance(model = model)
  if (nrow(imp_i) > 0) {
    all_importance <- rbind(all_importance,
                            data.frame(as.data.frame(imp_i), Fold = i))
  }
}

cat("\n\nLOO-CV completed!\n")
end_time <- Sys.time()
total_time <- difftime(end_time, start_time, units = "mins")
cat("Total time:", round(total_time, 2), "minutes\n\n")

# ==============================================================================
# 4. CALCULATE PERFORMANCE METRICS (LOO-CV)
# ==============================================================================

library(pROC)

# ------------------------------------------------------------------------------
# Prepare common variables
# ------------------------------------------------------------------------------
actual_cpc <- c(1, 2, 3, 5)[actual_outcomes + 1]
predicted_cpc <- c(1, 2, 3, 5)[loo_predicted_class + 1]

true_good <- as.numeric(actual_cpc <= 3)
prob_good <- loo_predictions[, 1] +
  loo_predictions[, 2] +
  loo_predictions[, 3]

# ==============================================================================
# 4A. PRIMARY OUTCOME — AUC (Good vs Poor) WITH DELONG CI
# ==============================================================================

roc_obj  <- roc(true_good, prob_good, quiet = TRUE)
auc_val  <- as.numeric(auc(roc_obj))
auc_ci   <- ci.auc(roc_obj, method = "delong")

auc_results <- data.frame(
  AUC  = round(auc_val, 3),
  CI_L = round(auc_ci[1], 3),
  CI_U = round(auc_ci[3], 3)
)


# ==============================================================================
# END — METRICS COMPLETE
# ==============================================================================

# Overview
auc_results
