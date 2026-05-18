# ---- Setup ---------------------------------------------------------------
library(caret)
library(tidyverse)
library(ModelMetrics)
library(recipes)
library(glmnet)
library(yardstick)
library(MLeval)    # optional: compare against yardstick)
library(doParallel)
library(VSURF)
library(prg)
library(yardstick)
library(MLeval)    # optional: compare against yardstick
library(doParallel)
library(VSURF)
library(prg)
library(CalibrationCurves)
library(dcurves)
library(readxl)
library(statip)
library(pheatmap)
library(janitor)
set.seed(123)

source("Functions_nested_CV(VSURF+elastic_net).R")

# Read xlsx and recode the outcome

df_ev = read_xlsx("D:/Packages/artur/covid_ev_don_orig.xlsx",
               sheet = 1) 

df_cov_clean = df_ev %>%
  clean_names() %>%
  mutate(outcome = factor(outcome, levels = c(4, 2), labels = c("No", "Yes")),
         sex = factor(sex, levels = c(1, 2), labels = c("Male", "Female"))
  ) %>%
  dplyr::select(outcome,everything(), -x1) %>%
  #drop rows where NA in more than 36 variables
  filter(rowSums(is.na(.)) <= 26) %>%
  #remove row #17
  filter(row_number() != 20)

df_cov_clean %>%
  is.na() %>%
  colSums() #NA in 16 variables

df_cov_clean %>%
  is.na() %>%
  rowSums() #NA in X rows

summary(df_cov_clean)

##################################################Correlation analysis
# ---- Helper: correlation filter ----------------------------------------------
# Log-transforms and standardizes numerics, computes Pearson correlation,
# then removes predictors with |r| >= 0.9. Returns the filtered *original*
# (untransformed) data frame so the modeling functions apply their own transforms.

remove_high_cor <- function(df, cutoff = 0.9, show_heatmap = TRUE, title = "") {
  
  # Scaled copy for correlation only
  df_scaled <- df |>
    mutate(across(where(is.numeric), ~ log(. + 0.001))) |>
    mutate(across(where(is.numeric), ~ scale(.) |> as.vector()))
  
  cor_mat <- df_scaled |>
    dplyr::select(where(is.numeric)) |>
    cor(method = "pearson", use = "pairwise.complete.obs")
  
  if (show_heatmap) {
    pheatmap(
      cor_mat,
      cluster_rows   = TRUE,
      cluster_cols   = TRUE,
      display_numbers = TRUE,
      number_format  = "%.2f",
      color = colorRampPalette(c("blue", "white", "red"))(50),
      main  = paste("Correlation Heatmap —", title)
    )
  }
  
  to_drop <- findCorrelation(cor_mat, cutoff = cutoff, names = TRUE, exact = TRUE, verbose = T)
  
  if (length(to_drop) > 0) {
    message(sprintf("[%s] Removing %d predictor(s) with |r| >= %.1f: %s",
                    title, length(to_drop), cutoff, paste(to_drop, collapse = ", ")))
    df <- df |> dplyr::select(-all_of(to_drop))
  } else {
    message(sprintf("[%s] No predictors exceed |r| >= %.1f — nothing removed.", title, cutoff))
  }
  
  df
}


# ---- Section 3: Correlation filter -------------------------------------------
covid_mod <- remove_high_cor(df_cov_clean, cutoff = 0.8, title = "B cells")

cat("Covid predictors after filter:", ncol(covid_mod) - 1, "\n")


# ---- Section 3b: Zero / near-zero variance check -----------------------------
# freqRatio     : ratio of most to second-most common value (higher = more skewed)
# percentUnique : % of unique values (lower = less informative)
# zeroVar / nzv : flags set by caret (TRUE = problematic)

check_variance <- function(df, title = "") {
  nzv <- df |>
    dplyr::select(where(is.numeric)) |>
    caret::nearZeroVar(saveMetrics = TRUE) |>
    tibble::rownames_to_column("variable") |>
    arrange(percentUnique)
  cat(sprintf("\n=== %s: Near-Zero Variance Check ===\n", title))
  print(nzv)
  invisible(nzv)
}

nzv_b <- check_variance(covid_mod, "Covid")

# select predictors with %unique more than 30
covid_mod <- covid_mod %>%
  dplyr::select(-all_of(nzv_b %>% filter(percentUnique < 30) %>% pull(variable)))

nzv_b <- check_variance(covid_mod, "Covid")

# ---- Section 3c: Distributions by outcome group (boxplots) -------------------

plot_distributions <- function(df, title = "") {
  df |>
    pivot_longer(where(is.numeric), names_to = "variable", values_to = "value") |>
    ggplot(aes(x = outcome, y = value)) +
    geom_boxplot(outlier.shape = 1) +
    facet_wrap(~ variable, scales = "free_y") +
    labs(title = paste("Predictor distributions by outcome —", title),
         x = NULL, y = NULL) +
    theme_bw(base_size = 11)
}

plot_distributions(covid_mod, "Covid")

#plot boxplots but numeric vars are yeo-johnson transformed and scaled
plot_distributions <- function(df, title = "") {
  
  rec <- recipes::recipe(~ ., data = df) %>%
    recipes::update_role(outcome, new_role = "outcome") %>%
    recipes::step_YeoJohnson(
      recipes::all_numeric_predictors()
    ) %>%
    recipes::step_scale(all_numeric_predictors()) %>%
    recipes::step_center(all_numeric_predictors())
  
  rec_prep <- recipes::prep(rec, training = df)
  
  df_transformed <- recipes::bake(rec_prep, new_data = df)
  
  df_transformed %>%
    dplyr::select(where(is.numeric), -outcome) %>%
    tidyr::pivot_longer(
      everything(),
      names_to = "variable",
      values_to = "value"
    ) %>%
    ggplot2::ggplot(ggplot2::aes(x = factor(variable), y = value)) +
    ggplot2::geom_boxplot(outlier.shape = 1) +
    ggplot2::labs(
      title = paste("Yeo-Johnson transformed predictor distributions —", title),
      x = NULL,
      y = NULL
    ) +
    ggplot2::theme_bw(base_size = 11)
}

plot_distributions(covid_mod, "Covid")

#winsorize outliers that are 3*IQR away
# Function to winsorize outliers
# winsorize_outliers <- function(df, outcome = NULL, k = 3) {
#   
#   num_cols <- df |>
#     dplyr::select(dplyr::where(is.numeric)) |>
#     names()
#   
#   if (!is.null(outcome)) {
#     num_cols <- setdiff(num_cols, outcome)
#   }
#   
#   df |>
#     dplyr::mutate(
#       dplyr::across(
#         dplyr::all_of(num_cols),
#         ~ {
#           q1 <- stats::quantile(.x, 0.25, na.rm = TRUE)
#           q3 <- stats::quantile(.x, 0.75, na.rm = TRUE)
#           iqr <- stats::IQR(.x, na.rm = TRUE)
#           
#           lower <- q1 - k * iqr
#           upper <- q3 + k * iqr
#           
#           dplyr::case_when(
#             .x < lower ~ lower,
#             .x > upper ~ upper,
#             TRUE ~ .x
#           )
#         }
#       )
#     )
# }
# 
# covid_mod_win <- winsorize_outliers(covid_mod, outcome = "outcome", k = 3)
# 
# plot_distributions(covid_mod_win)

#check what % of values out of the bounds
check_all_winsor_percent <- function(df, k = 3, outcome = NULL) {
  num_cols <- df |>
    dplyr::select(dplyr::where(is.numeric)) |>
    names()
  
  if (!is.null(outcome)) {
    num_cols <- setdiff(num_cols, outcome)
  }
  
  dplyr::bind_rows(lapply(num_cols, function(col) {
    x <- df[[col]]
    
    q1 <- quantile(x, 0.25, na.rm = TRUE)
    q3 <- quantile(x, 0.75, na.rm = TRUE)
    iqr <- IQR(x, na.rm = TRUE)
    
    lower <- q1 - k * iqr
    upper <- q3 + k * iqr
    
    data.frame(
      variable = col,
      lower_bound = as.numeric(lower),
      upper_bound = as.numeric(upper),
      percent_below_lower = mean(x < lower, na.rm = TRUE) * 100,
      percent_above_upper = mean(x > upper, na.rm = TRUE) * 100,
      total_percent_winsorized = mean(x < lower | x > upper, na.rm = TRUE) * 100
    )
  }))
}

check_all_winsor_percent(covid_mod, k =3, outcome = "outcome")

# ---- Section 3d: Univariate outlier detection (IQR-based) --------------------
# Flags observations > k*IQR beyond the quartiles on any single predictor.
# Applied to both datasets. MCD (multivariate) is not applicable to B cells
# because the all_19_* columns are compositional (sum to 100%), making the
# covariance matrix singular.

detect_outliers_univariate <- function(df, title = "", k = 3) {
  num_vars <- df |> dplyr::select(where(is.numeric)) |> names()
  
  flags <- map_dfr(num_vars, function(v) {
    x   <- df[[v]]
    q   <- quantile(x, c(0.25, 0.75), na.rm = TRUE)
    iqr <- q[2] - q[1]
    tibble(
      variable   = v,
      n_outliers = sum(x < q[1] - k * iqr | x > q[2] + k * iqr, na.rm = TRUE),
      rows       = list(which(x < q[1] - k * iqr | x > q[2] + k * iqr))
    )
  }) |>
    filter(n_outliers > 0)
  
  cat(sprintf("\n=== %s: Univariate outliers (> %d × IQR from quartiles) ===\n", title, k))
  if (nrow(flags) == 0) {
    cat("  No outliers detected.\n")
  } else {
    flags |>
      mutate(rows = map_chr(rows, ~ paste(.x, collapse = ", "))) |>
      print()
  }
  invisible(flags)
}

out_b <- detect_outliers_univariate(covid_mod, "covid")


# ---- Section 3e: MCD multivariate outlier detection (T cells only) -----------
# B cells: MCD not applicable due to compositional constraint in all_19_* cols.
# T cells: full-rank predictor matrix — MCD applicable.
# NOTE: with n=84 complete cases and 12 predictors, the chi-sq(0.975) cutoff
# is conservative. Use results to *inspect* flagged observations, not to
# automatically exclude them.

# detect_outliers_mcd <- function(df, title = "") {
#   data_num <- df |>
#     dplyr::select(where(is.numeric)) |>
#     tidyr::drop_na() |>
#     as.data.frame()
#   
#   n_dropped <- nrow(df) - nrow(data_num)
#   if (n_dropped > 0)
#     message(sprintf("[%s] MCD: %d row(s) with NAs excluded.", title, n_dropped))
#   
#   mcd    <- covMcd(data_num)
#   dists  <- sqrt(mcd$mah)
#   cutoff <- sqrt(qchisq(0.975, df = ncol(data_num)))
#   rows   <- which(complete.cases(df |> dplyr::select(where(is.numeric))))
#   
#   result <- tibble(row = rows, distance = dists, outlier = dists > cutoff)
#   
#   p <- ggplot(result, aes(x = row, y = distance, color = outlier)) +
#     geom_point(size = 2.5) +
#     geom_hline(yintercept = cutoff, linetype = "dashed", color = "red") +
#     scale_color_manual(values = c("FALSE" = "grey40", "TRUE" = "red"),
#                        labels = c("Inlier", "Outlier")) +
#     labs(title = paste("MCD Multivariate Outlier Detection —", title),
#          x = "Row index", y = "Robust Mahalanobis distance", color = NULL) +
#     theme_bw(base_size = 11)
#   print(p)
#   
#   n_out <- sum(result$outlier)
#   if (n_out > 0)
#     cat(sprintf("[%s] %d multivariate outlier(s) detected (rows: %s)\n",
#                 title, n_out, paste(result$row[result$outlier], collapse = ", ")))
#   else
#     cat(sprintf("[%s] No multivariate outliers detected.\n", title))
#   
#   invisible(result)
# }
# 
# # B cells: skip MCD (singular covariance — compositional predictors)
# mcd_b_res <- NULL
# # T cells: MCD applicable
# mcd_t_res <- detect_outliers_mcd(covid_mod, "T cells")


# ---- Section 3f: Missing data pattern (T cells) ------------------------------
# All 5 NAs in `th` fall in the "No" (donors) group — 14% of donors are missing
# this value, 0% of COVID patients. Missingness is not random (not MCAR).
# Median imputation (inside the modeling function) will use the full-sample
# median, which may slightly overestimate the donors' true `th` value.
# Consider whether this pattern is biological or administrative before modeling.

# cat("\n=== T cells: NA pattern in `th` by outcome ===\n")
# t_cells_mod |>
#   mutate(th_missing = is.na(th)) |>
#   count(outcome, th_missing) |>
#   pivot_wider(names_from = th_missing, values_from = n,
#               names_prefix = "NA_", values_fill = 0) |>
#   rename(present = NA_FALSE, missing = NA_TRUE) |>
#   mutate(pct_missing = round(100 * missing / (present + missing), 1)) |>
#   print()


# ---- Section 4a: Nested CV — B cells -----------------------------------------
results_b <- nested_elastic_binary_outcome(
  df                  = df_cov_clean,
  outcome_var         = "outcome",
  positive_class      = "Yes",
  negative_class      = "No",
  cv_outer_folds      = 5,
  cv_outer_repeats    = 20,
  seed                = 1,
  inner_cv_method     = "repeatedcv",
  inner_cv_folds      = 5,
  inner_cv_repeats    = 20,
  alpha_grid          = seq(0, 1, by = 0.1),
  lambda_grid         = 10^seq(-4, 1, length.out = 50),
  hyperparam_search   = "grid",
  ntree               = 1000,
  nforests            = 20,
  selection_rule      = "best",
  data_transformation = "YeoJohnson",
  numeric_imputation  = "median",
  nominal_imputation =  "unknown",
  add_interactions    = "no",
  add_polynomials     = "no",
  remove_zero_variance = "yes",
  add_missing_indicators = "no",
  remove_correlated_predictors = "no",
  correlation_threshold = 0.8,
)

# ---- Section 5a: Evaluate — B cells ------------------------------------------
perf_cv_5_20 <- outer_perf_nested_binary(
  trained_object = results_b,
  positive_class = "Yes",
  negative_class = "No"
)

# ROC curve
roc_b <- pROC::roc(
  response  = results_b$outer_predictions$outcome,
  predictor = results_b$outer_predictions$pred_lasso,
  levels    = c("No", "Yes"),
  direction = "<"
)
plot(roc_b, col = "blue", main = "ROC — B cells (Lasso)")
cat("B cells AUC:", round(auc(roc_b), 3), "\n")

# Calibration
y01_b <- as.integer(results_b$outer_predictions$outcome == "Yes")
cal_b <- val.prob.ci.2(
  results_b$outer_predictions$pred_lasso,
  y01_b,
  logistic  = TRUE,
  col.log   = "blue",
  main      = "Calibration — B cells",
  smooth = 'none'
)
cal_b


################################################################################Final model#######################################################################

# ---- Section 6a: Final model — B cells ---------------------------------------
final_b <- final_model_with_coefs(
  df                   = df_cov_clean,
  outcome_var          = "outcome",
  positive_class       = "Yes",
  negative_class       = "No",
  family               = "binomial",
  cv_method            = "repeatedcv",
  cv_folds             = 5,
  cv_repeats           = 20,
  alpha_grid           = seq(0, 1, by = 0.1),
  lambda_grid          = 10^seq(-4, 1, length.out = 50),
  ntree                = 1000,
  nforests             = 20,
  selection_rule      = "best",
  data_transformation = "YeoJohnson",
  numeric_imputation  = "median",
  nominal_imputation =  "unknown",
  add_interactions    = "no",
  add_polynomials     = "no",
  remove_zero_variance = "yes",
  add_missing_indicators = "no",
  remove_correlated_predictors = "no",
  correlation_threshold = 0.8,
  post_selection_refit = "none"
)
