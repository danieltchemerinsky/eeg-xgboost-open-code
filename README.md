# Predicting neurological outcome for comatose cardiac arrest patients from early EEG reactivity data by use of machine learning

Analysis code for Del Pin, Duez, Tchemerinsky Konieczny, Johnsen & Sandberg.

An XGBoost model is trained on EEG reactivity features recorded within 24 hours of cardiac
arrest and used to predict Cerebral Performance Category (CPC) at 6 months, both as a
four-class problem (CPC 1, 2, 3, 5) and dichotomised into good (CPC 1–3) versus poor
(CPC 4–5) outcome. Performance is estimated by leave-one-out cross-validation over all 114
patients.

## Reproducibility

Every model is seeded from its own fold index, so results do not depend on the seed's
position in the code, on anything else executed inside the loop, or on the number of cores
available. On R 4.4.1 / xgboost 3.1.3.1 (aarch64-darwin), a repeat run, a run with the RNG
stream advanced beforehand and a single-threaded run all produce bit-identical predictions.

## Contents

`run_all.R` executes these in the right order. Individually, they must be run in one
session, in the order listed, because each uses objects the previous ones leave behind.

| File | Purpose |
|---|---|
| **`run_all.R`** | **Entry point.** Runs every script below, in order, in one session, into `output/`. |
| `01_loocv_analysis.R` | Main analysis. Leave-one-out cross-validation, 114 models. Produces the per-patient out-of-sample class probabilities, the per-fold feature importance and the good-vs-poor AUC. Everything downstream is derived from its output. |
| `02_roc_reactivity_vs_all_predictors.R` | Control analysis for Figure 4. Repeats the LOO-CV with reactivity features alone and with reactivity plus four clinical predictors (shockable rhythm, bystander CPR, SSEP, age), for all 114 patients and for the 103 with complete SSEP data. Computes the SSEP operating point. |
| `03_figures.R` | Figure functions. `generate_figure2()` feature importance; `generate_figure3()` one-vs-all ROC; `generate_figure4_clean()` good-vs-poor ROC, two panels; `generate_figure_s1()` partial dependence. `generate_all_figures()` writes Figures 2, 3 and S1. |
| `04_tables_and_aucs.R` | Builds Table 1 (one-vs-all at a 50% threshold) and Table 2 (poor-outcome prediction at six criteria) from the cross-validated predictions, with Wilson score intervals, and computes the good-vs-poor and one-vs-all AUCs with DeLong intervals. Exports all three as CSVs. |
| `05_reported_statistics.R` | Every number quoted in the Results text that is not in a table or figure: the Hand & Till multiclass AUC, the Gain sums by sensory modality and frequency band, the importance ranking with clinical predictors included, and the algorithm-versus-SSEP comparison at matched criteria. Prints them and exports them as a CSV. |
| `06_render_table1.R` | Formats Table 1 for publication (`gt` / `flextable`). |
| `07_render_table2.R` | Formats Table 2 for publication (`gt` / `flextable`). |
| `08_partial_dependence_values.R` | Exports the partial dependence curves behind Figure S1, and summarises for each panel how much of the variation each outcome class accounts for — so a panel can be described from values rather than read off the plot. |
| `09_figure1_flowchart.R` | Patient inclusion flowchart (Figure 1). Standalone — derives the EEG, excluded, included and outcome counts from the data, and takes only the two upstream recruitment counts as declared constants. |

## Running it

Put the two data files in `data/` (see `data/README.md`), then from the
repository root:

```
Rscript run_all.R
```

That runs every script in order in one session and writes everything to
`output/`, including `sessionInfo.txt`. It checks for the data and the required
packages first, reports each step, and lists the expected outputs at the end.
About three minutes.

To run the scripts individually instead, work with `output/` as the working
directory (the scripts read and write relative paths) and keep the order above —
each uses objects left in the environment by the ones before it. Scripts `06`,
`07` and `09` are the exceptions: `06` and `07` need only the CSVs written by
`04`, and `09` needs only the data.

## Model

Multiclass XGBoost (`multi:softprob`, 4 classes), `tree_method = "hist"`, `eta = 0.01`,
`max_depth = 3`, `min_child_weight = 3`, `gamma = 0.1`, `subsample = 0.8`,
`colsample_bytree = 0.8`, `lambda = 1`, `alpha = 0`, 750 boosting rounds. Hyperparameters
are fixed across folds rather than tuned within each fold. The probability of a good outcome
is the sum of the predicted probabilities for CPC 1, 2 and 3. Scripts 01 and 02 use
identical settings, and their reactivity-only models agree exactly.

## Requirements

`xgboost`, `pROC`, `caret`, `dplyr`, `tidyr`, `ggplot2`, `patchwork`, `scales` and `Matrix`,
and for the table-rendering scripts `gt`, `flextable` and `officer`.

`webshot2` and `chromote` are optional: without them the tables are still written as `.docx`
and `.html`, only the `.png` previews are skipped.

A full run takes about three minutes.

## Data

**This repository contains no patient data, and none will be added to it.**

The scripts require two data files, placed in `data/` by whoever runs them. See
**`data/README.md`** for what each contains and which scripts use them.

The data come from the TTH48 multicentre trial (Kirkegaard et al., 2017). They are not
distributed here and cannot be deposited in any public repository: data-sharing agreements
between the participating centres and European data-protection law, including the GDPR, do
not permit open publication. `.gitignore` excludes `*.RData` and `output/`, so the
repository holds only code.

De-identified data may be requested from the corresponding author, subject to approval by
the TTH48 Trial Steering Committee and a data-access agreement under the host institutions'
governance framework. See the data availability statement in the paper.

## Licence

MIT — see `LICENSE`.

## Citation

See `CITATION.cff`. If you use this code, please cite both the software archive and the
paper.
