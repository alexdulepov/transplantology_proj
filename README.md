[![DOI](https://zenodo.org/badge/1022210779.svg)](https://doi.org/10.5281/zenodo.20277106)

# Nested Cross-Validation with VSURF and Elastic Net

This repository contains R utilities for repeated outer cross-validation that compare direct VSURF predictions, tuned elastic-net/lasso models, and constant reference models. Workflows are provided for binary, multiclass, and continuous outcomes, together with helpers for summarising outer performance and extracting final coefficient sets.

## Repository Contents

- `nested_cv_vsurf_elastic_net.R` - core implementation. Includes:
  - `nested_elastic_binary_outcome()`, `nested_elastic_multiclass_outcome()`, and `nested_elastic_continuous_outcome()` for repeated nested resampling.
  - `outer_perf_nested_*()` helpers that turn resampling objects into tidy performance summaries.
  - `final_model_with_coefs()` to fit full-data VSURF, elastic-net, and lasso models, return direct penalized glmnet coefficients, and optionally run post-selection refits.
- `demo_usage.R` - minimal synthetic example that sources the core functions and runs the binary workflow end-to-end.
- `covid_ev_don_models.R` - end-to-end analysis script for the COVID EV donor dataset. It imports and cleans the source workbook, performs missingness/correlation/variance/outlier checks, visualises predictor distributions, runs the binary nested-CV workflow, evaluates ROC/calibration, and fits the final full-data model.
- `requirements.R` - installs and loads the core package stack used by the reusable modelling functions and demo workflow.

## Installation

Install the required packages once per R environment:

```r
source("requirements.R")
```

`requirements.R` installs CRAN dependencies (VSURF, caret, glmnet, tidyverse, recipes, ModelMetrics, yardstick, MLeval, CalibrationCurves, dcurves, readxl, statip, pheatmap, janitor, doParallel, devtools, ranger, kknn, foreach, iterators) and the `prg` package from GitHub before loading everything into the session.

## Usage

### 1. Binary outcome workflow

```r
source("nested_cv_vsurf_elastic_net.R")

# Data frame 'df' must contain the outcome column plus predictors
# If using data_transformation = "log", numeric predictors must be non-negative.

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
    hyperparam_search = c("grid", "random"),
    random_search_size = 25,
    ntree = 1000,
    nforests = 20,
    inner_cv_method = c("LOOCV", "repeatedcv"),
    inner_cv_repeats = 20,
    inner_cv_folds = 5,
    selection_rule = c("best", "oneSE"),
    data_transformation = c("YeoJohnson", "log"),
    numeric_imputation = c("median", "bag", "linear"),
    nominal_imputation = c("mode", "unknown"),
    add_missing_indicators = c("no", "yes"),
    remove_zero_variance = c("no", "yes"),
    remove_correlated_predictors = c("no", "yes"),
    correlation_threshold = 0.9,
    add_interactions = c("no", "yes"),
    interaction_vars = NULL,
    add_polynomials = c("no", "yes"),
    poly_vars = NULL,
    poly_degree = 2L
)

# Summary
outer_perf  <- outer_perf_nested_binary(results, positive_class = "Yes")
```

`nested_elastic_binary_outcome()` returns nested CV predictions and metadata, including:

- `outer_predictions` - one row per original observation after averaging repeated outer-fold predictions. Columns include `row_id`, `outcome`, `pred_elastic_net`, `pred_lasso`, `pred_prevalence`, `pred_vsurf_interpretation`, and `pred_vsurf_prediction`.
- `prevalence` - the single prevalence value repeated for every binary reference prediction.
- `selected_variables` - for each outer resample, the variables chosen by VSURF, elastic net, and lasso.
- `tuning` - the actual elastic-net and lasso tuning candidates used in the run.
- `preprocessing` - the imputation, missing-indicator, optional zero-variance removal, optional correlated-predictor removal, interaction, and polynomial settings used in the run.

The helper `outer_perf_nested_binary()` also prints stability tables showing how frequently each predictor is selected across outer folds and returns a tidy tibble of outer performance metrics.

### 2. Multiclass outcome workflow

```r
multi_results <- nested_elastic_multiclass_outcome(
    df,
    outcome_var = "outcome",
    cv_outer_folds = 5,
    cv_outer_repeats = 20,
    seed = 1,
    alpha_grid  = seq(0, 1, by = 0.1),
    lambda_grid = 10^seq(-4, 1, length.out = 50),
    hyperparam_search = c("grid", "random"),
    random_search_size = 25,
    ntree = 1000,
    nforests = 20,
    inner_cv_method = c("LOOCV", "repeatedcv"),
    inner_cv_repeats = 20,
    inner_cv_folds = 5,
    selection_rule = c("best", "oneSE"),
    data_transformation = c("YeoJohnson", "log"),
    numeric_imputation = c("median", "bag", "linear"),
    nominal_imputation = c("mode", "unknown"),
    add_missing_indicators = c("no", "yes"),
    remove_zero_variance = c("no", "yes"),
    remove_correlated_predictors = c("no", "yes"),
    correlation_threshold = 0.9,
    add_interactions = c("no", "yes"),
    interaction_vars = NULL,
    add_polynomials = c("no", "yes"),
    poly_vars = NULL,
    poly_degree = 2L
)

outer_multi <- outer_perf_nested_multiclass(multi_results)
```

`nested_elastic_multiclass_outcome()` returns class-probability predictions in long format: one row per original observation and class. The performance helper reports accuracy, kappa, multiclass log loss, Hand-Till ROC AUC, and macro-averaged classification metrics.

### 3. Continuous outcome workflow

```r
reg_results <- nested_elastic_continuous_outcome(
    df,
    outcome_var = "outcome",
    cv_outer_folds = 5,
    cv_outer_repeats = 20,
    seed = 1,
    alpha_grid  = seq(0, 1, by = 0.1),
    lambda_grid = 10^seq(-4, 1, length.out = 50),
    hyperparam_search = c("grid", "random"),
    random_search_size = 25,
    ntree = 1000,
    nforests = 20,
    inner_cv_method = c("LOOCV", "repeatedcv"),
    inner_cv_repeats = 20,
    inner_cv_folds = 5,
    selection_rule = c("best", "oneSE"),
    optim_metric = c("RMSE", "MAE"),
    data_transformation = c("YeoJohnson", "log"),
    numeric_imputation = c("median", "bag", "linear"),
    nominal_imputation = c("mode", "unknown"),
    add_missing_indicators = c("no", "yes"),
    remove_zero_variance = c("no", "yes"),
    remove_correlated_predictors = c("no", "yes"),
    correlation_threshold = 0.9,
    add_interactions = c("no", "yes"),
    interaction_vars = NULL,
    add_polynomials = c("no", "yes"),
    poly_vars = NULL,
    poly_degree = 2L
)

outer_reg  <- outer_perf_nested_continuous(reg_results)
```

Continuous workflows mirror the binary setup but report RMSE, MAE, R2 (traditional and modern), concordance correlation, MAPE/SMAPE, and calibration diagnostics.

`nested_elastic_continuous_outcome()` returns the same high-level object names as the binary workflow. Its `outer_predictions` table uses `pred_baseline_mean` for the constant mean reference model.

### 4. Fit final models on the full dataset

Use `final_model_with_coefs()` when you are ready to fit the full-data VSURF, elastic-net, and lasso models after the resampling study:

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
hyperparam_search = c("grid", "random"),
random_search_size = 25,
ntree = 1000,
nforests = 20,
selection_rule = c("best", "oneSE"),
cont_optim_metric = c("RMSE", "MAE"),
data_transformation = c("log", "YeoJohnson"),
numeric_imputation = c("median", "bag", "linear"),
nominal_imputation = c("mode", "unknown"),
add_missing_indicators = c("no", "yes"),
remove_zero_variance = c("no", "yes"),
remove_correlated_predictors = c("no", "yes"),
correlation_threshold = 0.9,
add_interactions = c("no", "yes"),
interaction_vars = NULL,
add_polynomials = c("no", "yes"),
poly_vars = NULL,
poly_degree = 2L,
post_selection_refit = c("vsurf_only", "all", "none")
)

# Direct selected variables and penalized glmnet coefficients
final_bin$selected_variables
final_bin$coefficient_estimates$elastic_net_penalized
final_bin$coefficient_estimates$lasso_penalized

# Optional post-selection refit coefficients
final_bin$coefficient_estimates$post_selection
```

Set `family = "gaussian"` to obtain analogous continuous-outcome fits.

`final_model_with_coefs()` returns `post_selection_models`, `selected_variables`, `coefficient_estimates`, `vsurf`, `elastic_net`, `lasso`, `tuning`, and `preprocessing`. With the default `post_selection_refit = "vsurf_only"`, ElasticNet and Lasso are reported as their original full-data penalized glmnet fits rather than being refit after their own selection.

### 5. Demo script

`demo_usage.R` generates a toy binary dataset, runs `nested_elastic_binary_outcome()`, and produces summaries/plots illustrating how to interrogate the returned object.  Use it as a quick sanity check that all dependencies are installed.

### 6. COVID EV donor analysis script

`covid_ev_don_models.R` is the repository's current project-specific workflow. It:

- reads `covid_ev_don_orig.xlsx` from a local Excel path,
- recodes the binary outcome to `No` / `Yes` and sex to `Male` / `Female`,
- removes rows with extensive missingness and performs correlation, near-zero-variance, distribution, and outlier checks,
- runs the binary nested-CV workflow with Yeo-Johnson preprocessing, median numeric imputation, unknown-level nominal imputation, and optional zero-variance removal,
- evaluates Lasso predictions with ROC and calibration plots,
- fits a final binomial model using `final_model_with_coefs()`.

The input path is currently hard-coded inside the script, so update the `read_xlsx()` path before running it on another machine.

## Customisation Highlights

- **Variable selection** - toggle VSURF tree counts (`ntree`, `nforests`) and let elastic net or lasso perform their own selection by adjusting `alpha_grid`/`lambda_grid`.
- **Hyperparameter search** - use `hyperparam_search = "grid"` for the full alpha/lambda grid or `"random"` to sample `random_search_size` candidate combinations. Lasso is fit separately with `alpha = 1`; the elastic-net branch uses values in `alpha_grid` below 1.
- **Imputation** - use `numeric_imputation = "median"`, `"bag"`, or `"linear"` for numeric predictors; use `nominal_imputation = "mode"` or `"unknown"` for nominal predictors; set `add_missing_indicators = "yes"` to add missingness indicators before imputation. Linear imputation uses complete numeric predictors as auxiliary variables and falls back to median imputation when none are available in a resample.
- **Zero-variance predictors** - `remove_zero_variance = "no"` by default keeps all predictors in the recipe. Use `"yes"` only when you deliberately want the resampling recipe to drop predictors that are constant inside the training split.
- **Correlated predictors** - set `remove_correlated_predictors = "yes"` to remove highly correlated predictors among the original numeric input columns only. Use `correlation_threshold` to control the absolute-correlation cutoff; the default is `0.9`.
- **Interactions** - set `add_interactions = "yes"` and provide `interaction_vars = c("x1", "group", "x3")` to add pairwise interactions only among those predictors. Numeric and categorical predictors are supported; categorical predictors are dummy-expanded for `recipes::step_interact()` after imputation and before final scaling.
- **Polynomials** - set `add_polynomials = "yes"` and provide `poly_vars = c("x1", "x2")` to create polynomial basis terms with `recipes::step_poly()`. Use `poly_degree` to control the degree; polynomial variables must be numeric.
- **Post-selection refits** - use `post_selection_refit = "vsurf_only"` to refit only after VSURF selection, `"all"` to also refit after ElasticNet/Lasso selection, or `"none"` to return only the direct selected variables and penalized glmnet coefficients.
- **Hyperparameters selection rule** - choose which rule is applied for the hyperparameters selection during the cross-validation (`best` or `oneSE`).
- **Inner resampling** - choose between leave-one-out (`inner_cv_method = "LOOCV"`) or repeated K-fold cross-validation. The `oneSE` selection rule is available for repeated CV but not LOOCV.
- **Optimisation metric** - binary models optimise LogLoss; continuous models can optimise RMSE or MAE via `optim_metric`.
- **Parallelism** - parallelism incorporated by `allowParallel = TRUE` in trainControl and initializing parallel clusters.

## Data Preparation Notes

- Numeric predictors must be >= 0 to avoid failures in the `log` preprocessing step used by the recipes. If `YeoJohnson` is selected then numeric values can have any sign.
- Character predictors are converted to factors internally; ensure categorical values are coded consistently across rows.
- Inspect `covid_ev_don_models.R` for the current project-specific cleaning workflow. It includes explicit data-cleaning decisions before modelling, while the reusable modelling functions still control recipe-based imputation, optional zero-variance removal, and optional correlated-predictor removal during resampling.

## Data availability

The raw data are not included in this repository because they contain study participant data. Data are available from the corresponding author on reasonable request. The repository provides the analysis workflow and reusable modelling functions.

## Troubleshooting & Performance Tips

1. Ensure the outcome column contains both class labels specified via `positive_class`/`negative_class` before fitting binary models.
2. Large nested CV runs (many repeats or high-dimensional predictors) can be computationally intensive; reduce outer repeats or shrink the hyper-parameter grid when prototyping.
3. If VSURF fails because of limited variability, adjust or remove near-zero-variance predictors (see `covid_ev_don_models.R` for an example workflow).
4. When using the GitHub-hosted `prg` package, confirm that `devtools` is available and that you have internet access during installation.

## License

This code is provided for educational and research purposes.  Please ensure compliance with the licenses of the included packages.

## Support

For issues or questions:
1. Check the troubleshooting section above.
2. Verify package installations.
3. Review console messages for information about class ordering, selection rules, and preprocessing.
4. Ensure your data meets the format requirements described above.
