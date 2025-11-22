# Nested Cross-Validation with VSURF and Elastic Net

This repository contains a set of R utilities for running nested cross-validation pipelines that combine VSURF-based variable selection with elastic-net modelling.  Workflows are provided for both binary classification and continuous outcomes, together with helper functions for summarising resampling performance and extracting final coefficient sets.

## Repository Contents

- `Functions_nested_CV(VSURF+elastic_net).R` – core implementation.  Includes:
  - `nested_elastic_binary_outcome()` / `nested_elastic_continuous_outcome()` for repeated nested resampling.
  - `inner_perf_nested_*()` and `outer_perf_nested_*()` helpers that turn resampling objects into tidy performance summaries.
  - `final_model_with_coefs()` to refit a final elastic-net model (after VSURF or self-selection) on the full dataset and return coefficient estimates.
- `demo_usage.R` – minimal synthetic example that sources the core functions and runs the binary workflow end-to-end.
- `proj_transp.R` – end-to-end analysis script for a transplant cytokine dataset (data import, cleaning, imputation, exploratory analysis, and VSURF/elastic-net modelling).
- `requirements.R` – installs and loads the full package stack required by the above scripts.

## Installation

Install the required packages once per R environment:

```r
source("requirements.R")
```

`requirements.R` installs CRAN dependencies (VSURF, caret, glmnet, tidyverse, recipes, ModelMetrics, yardstick, MLeval, CalibrationCurves, dcurves, readxl, statip, pheatmap, janitor, doParallel, devtools, ranger, kknn, foreach, iterators, VIM) and the `prg` package from GitHub before loading everything into the session.

## Usage

### 1. Binary outcome workflow

```r
source("Functions_nested_CV(VSURF+elastic_net).R")

# Data frame 'df' must contain the outcome column plus predictors
# All numeric predictors should be non-negative because log1p preprocessing is applied

results <- nested_elastic_binary_outcome(
    df,
    outcome_var = "outcome",
    positive_class = "Yes",
    negative_class = "No",
    cv_outer_folds = 5,
    cv_outer_repeats = 20,
    seed = 1,
    alpha_grid  = seq(0, 1, by = 0.1),
    lambda_grid = 10^seq(-4, 1, length.out = 50),
    ntree = 1000,
    nforests = 20,
    inner_cv_method = c("LOOCV", "repeatedcv"),
    inner_cv_repeats = 20,
    inner_cv_folds = 5,
    selection_rule = c("best", "oneSE"),
    data_transformation = c("YeoJohnson", "log")
)

# Summaries
inner_perf  <- inner_perf_nested_binary(results)
outer_perf  <- outer_perf_nested_binary(results, positive_class = "Yes")
```

`nested_elastic_binary_outcome()` returns nested CV predictions and metadata, including:

- `avg_final_outer_preds_vs` / `avg_final_outer_preds_elas` – averaged outer-fold probabilities from VSURF-preselected and elastic-net-self-selected models alongside `outer_y` labels.
- `avg_inner_biased_preds_vs` / `avg_inner_biased_preds_elas` – averaged inner-loop predictions used for hyper-parameter selection.
- `avg_final_outer_preds_single_elas` – averaged outer-fold probabilities from elastic net without prior VSURF selection.
- `avg_final_outer_preds_prev` – averaged outer-fold probabilities from prevalence (baseline) model.
- `avg_final_outer_preds_vsurf_inter` – averaged outer-fold probabilities from VSURF using variables from the interpretation step.
- `avg_final_outer_preds_vsurf_pred` – averaged outer-fold probabilities from VSURF using variables from the prediction step.
- `sel_vars_df` – for each outer resample, the variables chosen by VSURF(interpretation and prediction steps) and elastic net plus the coefficients of the final refit.
- `inner_perf` – fold-level metrics (ROC AUC, PR AUC, AUPRG, LogLoss, Brier score, MCC, calibration slope/intercept, etc.).

The helper `outer_perf_nested_binary()` also prints stability tables showing how frequently each predictor is selected across outer folds and returns a tidy tibble of outer vs. inner performance metrics.

### 2. Continuous outcome workflow

```r
reg_results <- nested_elastic_continuous_outcome(
    df,
    outcome_var = "outcome",
    cv_outer_folds = 5,
    cv_outer_repeats = 20,
    seed = 1,
    alpha_grid  = seq(0, 1, by = 0.1),
    lambda_grid = 10^seq(-4, 1, length.out = 50),
    ntree = 1000,
    nforests = 20,
    inner_cv_method = c("LOOCV", "repeatedcv"),
    inner_cv_repeats = 20,
    inner_cv_folds = 5,
    selection_rule = c("best", "oneSE"),
    optim_metric = c("RMSE", "MAE"),
    data_transformation = c("YeoJohnson", "log")
)

inner_reg  <- inner_perf_nested_continuous(reg_results)
outer_reg  <- outer_perf_nested_continuous(reg_results)
```

Continuous workflows mirror the binary setup but report RMSE, MAE, R² (traditional and modern), concordance correlation, MAPE/SMAPE, and calibration diagnostics.

### 3. Fit final models on the full dataset

Use `final_model_with_coefs()` when you are ready to refit a model on all available observations after the resampling study:

```r
final_bin <- final_model_with_coefs(
df,
outcome_var = "outcome",
positive_class = "Yes",
negative_class = "No",
family = c("binomial", "gaussian"),
cv_method = c("repeatedcv", "LOOCV"),
cv_folds = 5,
cv_repeats = 20,
alpha_grid  = seq(0, 1, by = 0.1),
lambda_grid = 10^seq(-4, 1, length.out = 50),
ntree = 1000,
nforests = 20,
selection_rule = c("best", "oneSE"),
cont_optim_metric = c("RMSE", "MAE"),
data_transformation = c("log", "YeoJohnson")
)

# Access VSURF- and elastic-net-based coefficient sets
final_bin$sel_vars_df$final_coefs
```

Set `family = "gaussian"` to obtain analogous continuous-outcome fits.

### 4. Demo script

`demo_usage.R` generates a toy binary dataset, runs `nested_elastic_binary_outcome()`, and produces summaries/plots illustrating how to interrogate the returned object.  Use it as a quick sanity check that all dependencies are installed.

## Customisation Highlights

- **Variable selection** – toggle VSURF tree counts (`ntree`, `nforests`) and let elastic net perform its own selection by adjusting `alpha_grid`/`lambda_grid`.
- **Hyperparameters selection rule** – choose which rule is applied for the hyperparameters selection during the cross-validation (`best` or `oneSE`).
- **Inner resampling** – choose between leave-one-out (`inner_cv_method = "LOOCV"`) or repeated K-fold cross-validation.  The `oneSE` selection rule is available for repeated CV but not LOOCV.
- **Optimisation metric** – binary models optimise LogLoss; continuous models can optimise RMSE or MAE via `optim_metric`.
- **Parallelism** – parallelism incorporated by `allowParallel = TRUE` in trainControl and initializing parallel clusters. 

## Data Preparation Notes

- Numeric predictors must be ≥ 0 to avoid failures in the `log` preprocessing step used by the recipes. If `YeoJohnson` is selected then numeric values can have any sign.
- Character predictors are converted to factors internally; ensure categorical values are coded consistently across rows.
- Inspect `proj_transp.R` for a full cleaning pipeline that demonstrates missing-value imputation with `VIM::kNN`, low-variance filtering, correlation checks, and VSURF exploration prior to the nested CV workflow.

## Troubleshooting & Performance Tips

1. Ensure the outcome column contains both class labels specified via `positive_class`/`negative_class` before fitting binary models.
2. Large nested CV runs (many repeats or high-dimensional predictors) can be computationally intensive—reduce outer repeats or shrink the hyper-parameter grid when prototyping.
3. If VSURF fails because of limited variability, adjust or remove near-zero-variance predictors (see `proj_transp.R` for an example workflow).
4. When using the GitHub-hosted `prg` package, confirm that `devtools` is available and that you have internet access during installation.

## License

This code is provided for educational and research purposes.  Please ensure compliance with the licenses of the included packages.

## Support

For issues or questions:
1. Check the troubleshooting section above.
2. Verify package installations.
3. Review console messages for information about class ordering, selection rules, and preprocessing.
4. Ensure your data meets the format requirements described above.
