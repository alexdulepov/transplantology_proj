# Minimal demo for the binary outer-CV workflow.

source("Functions_nested_CV(VSURF+elastic_net).R")

set.seed(123)
n_samples <- 80
n_features <- 12

x <- as.data.frame(matrix(rnorm(n_samples * n_features), nrow = n_samples))
names(x) <- paste0("x", seq_len(n_features))

logit_prob <- 0.6 * x$x1 - 0.5 * x$x2 + 0.4 * x$x3
prob <- plogis(logit_prob)

demo_df <- x
demo_df$outcome <- factor(
  ifelse(rbinom(n_samples, 1, prob) == 1, "Yes", "No"),
  levels = c("No", "Yes")
)
demo_df$x2[sample(n_samples, 6)] <- NA_real_
demo_df$x5[sample(n_samples, 5)] <- NA_real_
demo_df$group <- sample(c("low", "high", NA), n_samples, replace = TRUE, prob = c(0.45, 0.45, 0.10))

cat("Class distribution:\n")
print(table(demo_df$outcome))

results <- nested_elastic_binary_outcome(
  demo_df,
  outcome_var = "outcome",
  positive_class = "Yes",
  negative_class = "No",
  cv_outer_folds = 3,
  cv_outer_repeats = 1,
  seed = 123,
  alpha_grid = c(0, 0.5, 1),
  lambda_grid = 10^seq(-3, 0, length.out = 6),
  hyperparam_search = "random",
  random_search_size = 4,
  ntree = 100,
  nforests = 5,
  inner_cv_method = "repeatedcv",
  inner_cv_repeats = 1,
  inner_cv_folds = 3,
  selection_rule = "best",
  data_transformation = "YeoJohnson",
  numeric_imputation = "median",
  nominal_imputation = "unknown",
  add_missing_indicators = "yes",
  remove_zero_variance = "no",
  remove_correlated_predictors = "no",
  correlation_threshold = 0.9,
  add_interactions = "yes",
  interaction_vars = c("x1", "x2", "group"),
  add_polynomials = "yes",
  poly_vars = c("x1", "x2"),
  poly_degree = 2
)

outer_perf <- outer_perf_nested_binary(
  results,
  positive_class = "Yes",
  negative_class = "No"
)

print(outer_perf)

cat("\nBaseline prevalence used for every row:", results$prevalence, "\n")
cat("Unique baseline predictions:", unique(results$outer_predictions$pred_prevalence), "\n")
