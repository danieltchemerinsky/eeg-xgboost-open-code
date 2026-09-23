# ==============================================================================
# REPORTED STATISTICS - the numbers quoted in the Results text
# ==============================================================================
# Everything in Results 3.2, 3.3, 3.5 and the closing sentence of 3.6 that is
# not in Table 1, Table 2 or a figure, so that every value quoted in the text
# has a script behind it.
#
# Run after 01 and 02, in the same session:
#
#   source("01_loocv_analysis.R")
#   source("02_roc_reactivity_vs_all_predictors.R")
#   source("05_reported_statistics.R")
#
# Every value is printed as it is computed and collected into
# reported_statistics.csv.
# ==============================================================================

library(xgboost)
library(pROC)

required <- c("loo_predictions", "actual_outcomes", "train_data", "train_label",
              "params", "nrounds", "full",            # from 01
              "data_114", "data_103", "pred_103_react",
              "cols_reactivity", "cols_clinical")     # from 02
missing <- required[!vapply(required, exists, logical(1))]
if (length(missing) > 0) {
  stop(paste("Missing objects:", paste(missing, collapse = ", "),
             "\n>>> Run 01_loocv_analysis.R and",
             "02_roc_reactivity_vs_all_predictors.R first."))
}

# Consistency checks. If a .RData file has overwritten one of these objects the
# failure otherwise surfaces deep inside xgboost as an unhelpful length error.
if (nrow(train_data) != length(train_label))
  stop(sprintf(paste("train_data has %d rows but train_label has %d values.",
                     "\n>>> Something has overwritten one of them; re-run",
                     "01_loocv_analysis.R in a clean session."),
               nrow(train_data), length(train_label)))
if (nrow(loo_predictions) != length(actual_outcomes))
  stop("loo_predictions and actual_outcomes disagree on the number of patients.")
if (!all(c(cols_reactivity, cols_clinical) %in% names(data_114)))
  stop("data_114 does not contain all of cols_reactivity and cols_clinical.")

reported <- data.frame(Section = character(), Quantity = character(),
                       Value = numeric(), stringsAsFactors = FALSE)

note <- function(section, quantity, value) {
  reported <<- rbind(reported, data.frame(
    Section = section, Quantity = quantity, Value = round(value, 3),
    stringsAsFactors = FALSE))
}

rule <- function(ch = "=") cat(strrep(ch, 78), "\n")


# ==============================================================================
# 3.3  MULTICLASS AUC (Hand & Till, 2001)
# ==============================================================================
rule(); cat("3.3  Multiclass AUC (Hand & Till)\n"); rule()

# A(i|j): the probability that a randomly drawn member of class j has a lower
# estimated probability of belonging to class i than a randomly drawn member of
# class i. This is the Mann-Whitney statistic, so it can be read off the ranks.
A_given <- function(prob_i, is_i, is_j) {
  xi <- prob_i[is_i]; xj <- prob_i[is_j]
  ni <- length(xi);   nj <- length(xj)
  if (ni == 0 || nj == 0) return(NA_real_)
  r <- rank(c(xi, xj))
  (sum(r[seq_len(ni)]) - ni * (ni + 1) / 2) / (ni * nj)
}

# Hand & Till M: the mean of A_hat(i,j) = (A(i|j) + A(j|i)) / 2 over all class
# pairs. Unlike a prevalence-weighted average this is insensitive to class size.
hand_till <- function(pred, truth, classes) {
  pairs <- combn(classes, 2)
  vals <- apply(pairs, 2, function(p) {
    i <- p[1]; j <- p[2]
    is_i <- truth == i; is_j <- truth == j
    # Index by column NAME. The columns are labelled by CPC (1, 2, 3, 5), so
    # numeric indexing would ask for column 5 of a 4-column matrix.
    mean(c(A_given(pred[, as.character(i)], is_i, is_j),
           A_given(pred[, as.character(j)], is_j, is_i)))
  })
  list(M = mean(vals), values = vals)
}

# Columns of loo_predictions are classes 0,1,2,3 = CPC 1,2,3,5. Work on a copy so
# that relabelling the columns cannot affect anything sourced after this script.
cpc_of_col <- c(1, 2, 3, 5)
pred_cpc <- loo_predictions
colnames(pred_cpc) <- as.character(cpc_of_col)
actual_cpc <- cpc_of_col[actual_outcomes + 1]

ht_3 <- hand_till(pred_cpc, actual_cpc, c(1, 2, 5))
ht_4 <- hand_till(pred_cpc, actual_cpc, c(1, 2, 3, 5))

cat("\nPairwise AUCs (CPC 1, 2, 5):\n")
prs <- combn(c(1, 2, 5), 2)
for (k in seq_len(ncol(prs)))
  cat(sprintf("  CPC %d vs CPC %d : %.3f\n", prs[1, k], prs[2, k], ht_3$values[k]))

cat(sprintf("\n  Hand & Till M over CPC 1, 2, 5    : %.3f\n", ht_3$M))
cat(sprintf("  Hand & Till M over all four CPC   : %.3f\n", ht_4$M))

# pROC's own implementation, as an independent check. Its calling convention has
# varied between versions, so a failure here must not stop the script.
mc_auc <- tryCatch(
  as.numeric(auc(multiclass.roc(response  = actual_cpc,
                                predictor = as.data.frame(pred_cpc),
                                quiet = TRUE))),
  error = function(e) NA_real_)
if (is.na(mc_auc)) {
  cat("  pROC::multiclass.roc: unavailable in this pROC version\n\n")
} else {
  cat(sprintf("  pROC::multiclass.roc (all classes): %.3f\n\n", mc_auc))
}

note("3.3", "Multiclass AUC (Hand & Till, CPC 1/2/5)", ht_3$M)


# ==============================================================================
# 3.2  FEATURE IMPORTANCE - 18 reactivity predictors
# ==============================================================================
rule(); cat("3.2  Feature importance, 18 reactivity predictors\n"); rule()

# Single model fitted to all patients, matching the model behind Figure 2.
dtrain_18 <- xgb.DMatrix(data = as.matrix(train_data), label = train_label,
                         missing = NA)
set.seed(42)
model_18 <- xgb.train(params = params, data = dtrain_18, nrounds = nrounds,
                      verbose = 0)
imp_18 <- xgb.importance(model = model_18)

if (!all(c("Feature", "Gain") %in% names(imp_18))) {
  stop("xgb.importance() returned columns: ", paste(names(imp_18), collapse = ", "),
       "\n>>> This script, and 03_figures.R, expect 'Feature' and 'Gain'.")
}

cat("\nTop 8 features by Gain:\n")
print(head(as.data.frame(imp_18)[, c("Feature", "Gain")], 8), row.names = FALSE)

top1 <- imp_18$Gain[1]; top2 <- imp_18$Gain[2]
top6 <- sum(imp_18$Gain[1:6])
cat(sprintf("\n  Most important (%s) Gain : %.3f\n",
            imp_18$Feature[1], top1))
cat(sprintf("  Second      (%s) Gain    : %.3f\n",
            imp_18$Feature[2], top2))
cat(sprintf("  Top 6 share of total Gain      : %.3f\n", top6))

note("3.2", paste0("Gain, ", imp_18$Feature[1]), top1)
note("3.2", paste0("Gain, ", imp_18$Feature[2]), top2)
note("3.2", "Top 6 predictors, share of Gain", top6)

# Stimulation site -> sensory modality. Four of the six sites are painful.
modality_of <- function(f) {
  ifelse(grepl("LeftArm|RightArm|LeftLeg|Stern", f), "pain",
  ifelse(grepl("Eye", f), "visual",
  ifelse(grepl("Sound", f), "auditory", NA)))
}
band_of <- function(f) sub("\\..*$", "", f)

mod_sum  <- tapply(imp_18$Gain, modality_of(imp_18$Feature), sum)
band_sum <- tapply(imp_18$Gain, band_of(imp_18$Feature), sum)

cat("\nSummed Gain by sensory modality:\n")
for (m in c("pain", "visual", "auditory"))
  cat(sprintf("  %-9s : %.3f\n", m, ifelse(is.na(mod_sum[m]), 0, mod_sum[m])))

cat("\nSummed Gain by frequency band:\n")
for (b in c("Theta", "Alpha", "Delta"))
  cat(sprintf("  %-9s : %.3f\n", b, ifelse(is.na(band_sum[b]), 0, band_sum[b])))
cat("\n")

note("3.2", "Summed Gain, pain",     as.numeric(mod_sum["pain"]))
note("3.2", "Summed Gain, visual",   as.numeric(mod_sum["visual"]))
note("3.2", "Summed Gain, auditory", as.numeric(mod_sum["auditory"]))
note("3.2", "Summed Gain, Theta",    as.numeric(band_sum["Theta"]))
note("3.2", "Summed Gain, Alpha",    as.numeric(band_sum["Alpha"]))
note("3.2", "Summed Gain, Delta",    as.numeric(band_sum["Delta"]))

# Modality composition of the top 6, quoted as "2 pain, 2 visual, 2 auditory"
cat("Modality of the top 6 features: ",
    paste(sprintf("%s (%s)", imp_18$Feature[1:6],
                  modality_of(imp_18$Feature[1:6])), collapse = ", "), "\n\n")


# ==============================================================================
# 3.5  FEATURE IMPORTANCE - 18 reactivity + 4 clinical predictors
# ==============================================================================
rule(); cat("3.5  Feature importance, reactivity plus clinical predictors\n"); rule()

features_all <- c(cols_reactivity, cols_clinical)
dtrain_22 <- xgb.DMatrix(
  data   = as.matrix(data_114[, features_all]),
  label  = c("1" = 0, "2" = 1, "3" = 2, "5" = 3)[as.character(data_114$Outcome)],
  missing = NA
)
set.seed(42)
model_22 <- xgb.train(params = params, data = dtrain_22, nrounds = nrounds,
                      verbose = 0)
imp_22 <- xgb.importance(model = model_22)

cat("\nFull importance ranking:\n")
print(as.data.frame(imp_22)[, c("Feature", "Gain")], row.names = TRUE)

gain_of <- function(f) {
  g <- imp_22$Gain[imp_22$Feature == f]
  if (length(g) == 0) 0 else g          # absent = never selected = zero Gain
}
rank_of <- function(f) {
  r <- which(imp_22$Feature == f)
  if (length(r) == 0) NA_integer_ else r
}

top6_22 <- sum(imp_22$Gain[1:6])
cat(sprintf("\n  Top 6 share of total Gain : %.3f\n", top6_22))
cat(sprintf("  Age: rank %s, Gain %.3f\n", rank_of("age"), gain_of("age")))
cat(sprintf("  Shockable rhythm: Gain %.3f\n", gain_of("shockable_rhythm")))
cat(sprintf("  SSEP: Gain %.3f\n", gain_of("SEP")))
cat(sprintf("  Basic CPR: Gain %.3f\n\n", gain_of("basic_cpr")))

note("3.5", "Top 6 predictors, share of Gain", top6_22)
note("3.5", "Age, rank",  rank_of("age"))
note("3.5", "Age, Gain",  gain_of("age"))
note("3.5", "Shockable rhythm, Gain", gain_of("shockable_rhythm"))
note("3.5", "SSEP, Gain", gain_of("SEP"))
note("3.5", "Basic CPR, Gain", gain_of("basic_cpr"))


# ==============================================================================
# 3.6  ALGORITHM VERSUS SSEP AT MATCHED CRITERIA (N = 103)
# ==============================================================================
rule(); cat("3.6  Poor-outcome predictions vs SSEP, N = 103\n"); rule()

truth_poor_103 <- ifelse(data_103$Outcome == 5, 1, 0)

sweep <- do.call(rbind, lapply(c(0.1, 0.2, 0.3, 0.4, 0.5), function(crit) {
  pred_poor <- as.numeric(pred_103_react < crit)
  TP <- sum(pred_poor == 1 & truth_poor_103 == 1)
  FP <- sum(pred_poor == 1 & truth_poor_103 == 0)
  TN <- sum(pred_poor == 0 & truth_poor_103 == 0)
  data.frame(Criterion = crit, N_predicted_poor = TP + FP, TP = TP, FP = FP,
             Specificity = round(TN / sum(truth_poor_103 == 0), 3))
}))

cat("\nAlgorithm (reactivity only, N = 103):\n")
print(sweep, row.names = FALSE)

ssep_n   <- sum(data_103$SEP == 1)
ssep_tp  <- sum(data_103$SEP == 1 & data_103$Outcome == 5)
ssep_spec <- sum(data_103$SEP != 1 & truth_poor_103 == 0) / sum(truth_poor_103 == 0)

cat(sprintf("\nSSEP: %d patients predicted poor, %d correct, specificity %.3f\n",
            ssep_n, ssep_tp, ssep_spec))

row_04 <- sweep[sweep$Criterion == 0.4, ]
cat(sprintf("\n  At criterion 0.40: algorithm %d vs SSEP %d predicted poor\n",
            row_04$N_predicted_poor, ssep_n))
cat(sprintf("  Specificity: %.3f vs %.3f\n\n",
            row_04$Specificity, ssep_spec))

note("3.6", "Algorithm predicted poor at criterion 0.40", row_04$N_predicted_poor)
note("3.6", "SSEP predicted poor", ssep_n)
note("3.6", "Algorithm specificity at criterion 0.40", row_04$Specificity)
note("3.6", "SSEP specificity", ssep_spec)


# ==============================================================================
# SUMMARY
# ==============================================================================
rule(); cat("SUMMARY - every value quoted in the Results text\n"); rule(); cat("\n")

print(reported, row.names = FALSE)

write.csv(reported, "reported_statistics.csv", row.names = FALSE)
write.csv(as.data.frame(imp_18), "feature_importance_18_predictors.csv", row.names = FALSE)
write.csv(as.data.frame(imp_22), "feature_importance_22_predictors.csv", row.names = FALSE)
write.csv(sweep, "ssep_comparison_N103.csv", row.names = FALSE)

cat("\nSaved: reported_statistics.csv,",
    "feature_importance_18_predictors.csv,\n",
    "       feature_importance_22_predictors.csv, ssep_comparison_N103.csv\n")
rule()
