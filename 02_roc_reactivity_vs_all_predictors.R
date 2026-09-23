library(xgboost)
library(pROC)
library(caret)
library(dplyr)

# ==============================================================================
# 1. DATA PREPARATION
# ==============================================================================

# Load Data
#
# Into a private environment, taking only `combined`. This file is a saved 2019
# workspace and contains other objects besides the data; loading it into the
# global environment silently overwrites anything of the same name that
# 01_loocv_analysis.R created, which breaks any script sourced afterwards.
.data_env <- new.env()
load("AllPatientsEEGandETCAbsentObs5.RData", envir = .data_env)
combined <- .data_env$combined
rm(.data_env)

# Clean Data (N = 114)
absentobs <- rowSums(is.na(combined))
data_114 <- combined[absentobs <= 5, ]

# Create Subset (N = 103) - Only rows where SEP is NOT missing
data_103 <- data_114[!is.na(data_114$SEP), ]

# Define Feature Sets
# Clinical variables based on your str() output
cols_clinical <- c("shockable_rhythm", "basic_cpr", "SEP", "age")
cols_meta     <- c("Outcome", "Patient")

# Reactivity variables are everything else
cols_reactivity <- setdiff(names(data_114), c(cols_clinical, cols_meta))

# ==============================================================================
# 2. LOO-CV FUNCTION
# ==============================================================================

run_loo <- function(df, features) {
  n <- nrow(df)
  preds <- numeric(n) # Store probability of Good Outcome (CPC 1-3)
  
  # Map Outcome: 1,2,3 -> Good (0); 5 -> Poor (1)
  # We train XGBoost to predict the specific class (0,1,2,3), then sum probs
  y_true_class <- c("1"=0, "2"=1, "3"=2, "5"=3)[as.character(df$Outcome)]
  
  # Hyperparameters, identical to those in 01_loocv_analysis.R. Every one is set
  # explicitly, including those that happen to match an xgboost default: an
  # omitted setting here would silently fit a different model from the main
  # analysis while both are described in the paper as the reactivity model.
  params <- list(
    objective = "multi:softprob",
    eval_metric = "mlogloss",
    num_class = 4,
    tree_method = "hist",
    eta = 0.01,
    max_depth = 3,
    min_child_weight = 3,
    gamma = 0.1,
    subsample = 0.8,
    colsample_bytree = 0.8,
    lambda = 1,
    alpha = 0
  )

  cat("Running LOO-CV on N =", n, "with", length(features), "features...\n")
  pb <- txtProgressBar(min = 0, max = n, style = 3)

  # Seeded per fold, as in 01_loocv_analysis.R. Seeding from the fold index
  # rather than once before the loop makes each fit reproducible on its own and
  # independent of anything else executed inside the loop.
  LOO_SEED_BASE <- 42

  for(i in 1:n) {

    set.seed(LOO_SEED_BASE + i)
    # Leave one out
    train_idx <- setdiff(1:n, i)
    
    dtrain <- xgb.DMatrix(data = as.matrix(df[train_idx, features]),
                          label = y_true_class[train_idx], missing = NA)
    dtest  <- xgb.DMatrix(data = as.matrix(df[i, features, drop=FALSE]), missing = NA)
    
    # Train
    model <- xgb.train(params = params, data = dtrain, nrounds = 750, verbose = 0)
    
    # Predict
    p_vec <- predict(model, dtest) # Returns vector of length 4
    
    # Sum prob of CPC 1, 2, 3 (Indices 1, 2, 3 in R vector)
    # Class 0=CPC1, Class 1=CPC2, Class 2=CPC3, Class 3=CPC5
    preds[i] <- p_vec[1] + p_vec[2] + p_vec[3]
    
    setTxtProgressBar(pb, i)
  }
  close(pb)
  return(preds)
}

# ==============================================================================
# 3. RUN MODELS
# ==============================================================================

# --- LEFT PANEL (N=114) ---
# Model A: Reactivity Only
pred_114_react <- run_loo(data_114, cols_reactivity)
# Model B: All Predictors
pred_114_all   <- run_loo(data_114, c(cols_reactivity, cols_clinical))

# --- RIGHT PANEL (N=103) ---
# Model C: Reactivity Only (on subset)
pred_103_react <- run_loo(data_103, cols_reactivity)
# Model D: All Predictors (on subset)
pred_103_all   <- run_loo(data_103, c(cols_reactivity, cols_clinical))

# ==============================================================================
# 4. CALCULATE ROCs
# ==============================================================================

# Helper to get ROC (Target: Poor Outcome)
# Predictor is Prob_Good. Low Prob_Good => High Prob_Poor.
get_roc_obj <- function(df, preds) {
  # Binary Truth: Poor (CPC 5) = 1, Good (CPC 1-3) = 0
  truth_poor <- ifelse(df$Outcome == 5, 1, 0)
  roc(truth_poor, preds, quiet = TRUE)
}

roc_114_react <- get_roc_obj(data_114, pred_114_react)
roc_114_all   <- get_roc_obj(data_114, pred_114_all)
roc_103_react <- get_roc_obj(data_103, pred_103_react)
roc_103_all   <- get_roc_obj(data_103, pred_103_all)

# ==============================================================================
# 5. CALCULATE SEP POINT (For N=103)
# ==============================================================================

# Logic: SEP=1 (Absent/Abnormal) predicts Poor Outcome
# SEP=0 means present/normal, SEP=1 means absent/abnormal
# Truth: Outcome=5 is Poor (1)
sep_pred_poor <- ifelse(data_103$SEP == 1, 1, 0)
truth_poor_103 <- ifelse(data_103$Outcome == 5, 1, 0)

cm_sep <- confusionMatrix(factor(sep_pred_poor, levels=c(0,1)), 
                          factor(truth_poor_103, levels=c(0,1)), 
                          positive="1")

sep_sens <- cm_sep$byClass["Sensitivity"]
sep_spec <- cm_sep$byClass["Specificity"]

# ==============================================================================
# 5b. EXPORT THE FIGURE 4 AUCs WITH 95% CI
# ==============================================================================
# The AUC of the N = 114 reactivity model is the primary result quoted in the
# paper. It is exported here at full precision with its DeLong confidence
# interval, so the value in the text and the value in the figure legend come
# from one source rather than from two roundings of it.

auc_row <- function(label, roc_obj) {
  ci <- ci.auc(roc_obj, method = "delong")
  data.frame(
    Model = label,
    N     = length(roc_obj$response),
    AUC   = round(as.numeric(auc(roc_obj)), 3),
    CI_L  = round(ci[1], 3),
    CI_U  = round(ci[3], 3),
    CI_method = "DeLong",
    stringsAsFactors = FALSE
  )
}

figure4_auc <- rbind(
  auc_row("Reactivity only, N = 114",      roc_114_react),
  auc_row("All predictors, N = 114",       roc_114_all),
  auc_row("Reactivity only, N = 103",      roc_103_react),
  auc_row("All predictors, N = 103",       roc_103_all)
)

print(figure4_auc, row.names = FALSE)

cat(sprintf("\nSSEP alone (N = 103): sensitivity = %.3f, specificity = %.3f\n",
            sep_sens, sep_spec))

write.csv(figure4_auc, "Figure4_AUC_with_CI.csv", row.names = FALSE)

write.csv(
  data.frame(Measure = c("SSEP sensitivity", "SSEP specificity"),
             Value   = round(c(sep_sens, sep_spec), 3)),
  "Figure4_SSEP_point.csv", row.names = FALSE)

cat("Saved: Figure4_AUC_with_CI.csv, Figure4_SSEP_point.csv\n")

library(pROC)
library(caret)

# ==============================================================================
# PLOT FIGURE 4 (With Manual Legend Placement)
# ==============================================================================

# Setup the canvas
par(mfrow = c(1, 2), pty = "s", mar = c(4, 4, 3, 1))

# ------------------------------------------------------------------------------
# LEFT PANEL (N=114)
# ------------------------------------------------------------------------------
plot(roc_114_react, 
     col = "black", lty = 2, lwd = 3, 
     xlim = c(1, 0), ylim = c(0, 1),
     xlab = paste0("Specificity (N=", nrow(data_114), ")"), 
     ylab = "Sensitivity", 
     main = "", print.auc = FALSE)

lines(roc_114_all, col = "blue", lty = 4, lwd = 3)

# Manual Legend Placement
# x=0.6 starts the box at Specificity 0.6 (towards the right side)
# y=0.25 starts the box at Sensitivity 0.25
legend(x = 0.85, y = 0.35, 
       legend = c(sprintf("Reactivity Alone | AUC: %.3f", auc(roc_114_react)),
                  sprintf("All Predictors | AUC: %.3f", auc(roc_114_all))),
       col = c("black", "blue"),
       lty = c(2, 4), lwd = 3, 
       bty = "n",       # No box border
       cex = 0.6,       # Text size
       adj = 0)         # Left-align text

# ------------------------------------------------------------------------------
# RIGHT PANEL (N=103)
# ------------------------------------------------------------------------------
plot(roc_103_react, 
     col = "black", lty = 2, lwd = 3, 
     xlim = c(1, 0), ylim = c(0, 1),
     xlab = paste0("Specificity (N=", nrow(data_103), ")"), 
     ylab = "Sensitivity", 
     main = "", print.auc = FALSE)

lines(roc_103_all, col = "blue", lty = 4, lwd = 3)

# Add SSEP Point (High Specificity ~1.0, Low Sensitivity)
# This point will appear in the BOTTOM-LEFT area
points(sep_spec, sep_sens, pch = 19, col = "green", cex = 1.5)

# Manual Legend Placement
legend(x = 0.85, y = 0.35, 
       legend = c(sprintf("Reactivity Alone | AUC: %.3f", auc(roc_103_react)),
                  sprintf("All Predictors | AUC: %.3f", auc(roc_103_all)),
                  "SEP | AUC: N/A"),
       col = c("black", "blue", "green"),
       lty = c(2, 4, NA), 
       pch = c(NA, NA, 19),
       lwd = c(3, 3, NA), 
       bty = "n", 
       cex = 0.6,
       adj = 0)

# ------------------------------------------------------------------------------
# CAPTION
# ------------------------------------------------------------------------------
cat("\nFigure 4 generated with corrected legend placement.\n")