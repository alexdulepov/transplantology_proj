make_glmnet_tune_grid <- function(alpha_values,
                                  lambda_values,
                                  hyperparam_search,
                                  random_search_size,
                                  model_name = "glmnet") {
  # Keep all glmnet tuning candidates in the same caret-friendly shape.
  # Random search samples lambda on the log scale because useful penalties
  # usually span orders of magnitude rather than equal-width intervals.
  alpha_values <- sort(unique(alpha_values))
  lambda_values <- sort(unique(lambda_values))

  if (!length(alpha_values)) {
    stop(sprintf("No alpha values supplied for %s.", model_name))
  }
  if (!length(lambda_values) || any(lambda_values <= 0, na.rm = TRUE)) {
    stop("lambda_grid must contain positive values.")
  }

  if (hyperparam_search == "grid") {
    return(expand.grid(alpha = alpha_values, lambda = lambda_values))
  }

  random_search_size <- as.integer(random_search_size)
  if (is.na(random_search_size) || random_search_size < 1) {
    stop("random_search_size must be a positive integer.")
  }

  if (length(alpha_values) == 1) {
    alpha <- rep(alpha_values, random_search_size)
  } else {
    alpha <- stats::runif(random_search_size, min(alpha_values), max(alpha_values))
  }

  lambda <- 10^stats::runif(
    random_search_size,
    log10(min(lambda_values)),
    log10(max(lambda_values))
  )

  data.frame(alpha = alpha, lambda = lambda)
}

add_imputation_steps <- function(rec,
                                 data,
                                 outcome_var,
                                 numeric_imputation,
                                 nominal_imputation,
                                 add_missing_indicators) {
  predictor_names <- setdiff(names(data), outcome_var)
  is_numeric_predictor <- vapply(data[predictor_names], is.numeric, TRUE)
  complete_numeric_predictors <- predictor_names[
    is_numeric_predictor &
      !vapply(data[predictor_names], function(x) any(is.na(x)), TRUE)
  ]
  missing_numeric_predictors <- predictor_names[
    is_numeric_predictor &
      vapply(data[predictor_names], function(x) any(is.na(x)), TRUE)
  ]

  # Missingness indicators must be created before any imputation removes the NAs.
  if (add_missing_indicators == "yes") {
    rec <- rec |>
      recipes::step_indicate_na(recipes::all_predictors())
  }

  rec <- if (nominal_imputation == "unknown") {
    rec |>
      recipes::step_unknown(recipes::all_nominal_predictors(), new_level = "unknown")
  } else {
    rec |>
      recipes::step_impute_mode(recipes::all_nominal_predictors())
  }

  if (numeric_imputation == "median") {
    rec <- rec |>
      recipes::step_impute_median(recipes::all_numeric_predictors())
  } else if (numeric_imputation == "bag") {
    rec <- rec |>
      recipes::step_impute_bag(recipes::all_numeric_predictors())
  } else if (numeric_imputation == "linear") {
    if (!length(missing_numeric_predictors)) {
      rec <- rec |>
        recipes::step_impute_median(recipes::all_numeric_predictors())
    } else if (!length(complete_numeric_predictors)) {
      warning(
        "linear numeric imputation needs complete numeric auxiliary predictors; ",
        "falling back to median imputation."
      )
      rec <- rec |>
        recipes::step_impute_median(recipes::all_numeric_predictors())
    } else {
      # Linear imputation cannot use predictors that still contain missing values.
      # Using complete numeric predictors keeps this step stable across resamples.
      rec <- rec |>
        recipes::step_impute_linear(
          tidyselect::any_of(missing_numeric_predictors),
          impute_with = tidyselect::any_of(complete_numeric_predictors)
        ) |>
        # Median is a safety net for any new-data missingness in numeric columns
        # that were complete in this training split.
        recipes::step_impute_median(recipes::all_numeric_predictors())
    }
  } else {
    stop("numeric_imputation must be one of 'median', 'bag', or 'linear'.")
  }

  rec
}

validate_correlation_threshold <- function(correlation_threshold) {
  if (!is.numeric(correlation_threshold) ||
      length(correlation_threshold) != 1 ||
      is.na(correlation_threshold) ||
      correlation_threshold <= 0 ||
      correlation_threshold > 1) {
    stop("correlation_threshold must be a single numeric value in the interval (0, 1].")
  }

  as.numeric(correlation_threshold)
}

add_original_numeric_correlation_filter <- function(rec,
                                                    data,
                                                    outcome_var,
                                                    remove_correlated_predictors,
                                                    correlation_threshold) {
  if (remove_correlated_predictors == "no") {
    return(rec)
  }

  # Restrict correlation filtering to predictors that were numeric in the
  # original input data. Engineered terms are intentionally not candidates.
  original_numeric_predictors <- numeric_predictor_names(data, outcome_var)
  if (!length(original_numeric_predictors)) {
    return(rec)
  }

  rec |>
    recipes::step_corr(
      tidyselect::any_of(original_numeric_predictors),
      threshold = correlation_threshold
    )
}

finish_recipe_preprocessing <- function(rec,
                                        data,
                                        outcome_var,
                                        remove_zero_variance,
                                        remove_correlated_predictors,
                                        correlation_threshold) {
  # Zero-variance removal changes the design matrix, so keep it explicit.
  # If disabled, recipes may warn when a resample has constant engineered columns.
  if (remove_zero_variance == "yes") {
    rec <- rec |>
      recipes::step_zv(recipes::all_predictors())
  }

  rec <- add_original_numeric_correlation_filter(
    rec,
    data,
    outcome_var,
    remove_correlated_predictors,
    correlation_threshold
  )

  rec |>
    recipes::step_center(recipes::all_numeric_predictors()) |>
    recipes::step_scale(recipes::all_numeric_predictors())
}

backtick_name <- function(x) {
  paste0("`", gsub("`", "\\\\`", x), "`")
}

single_quote_string <- function(x) {
  paste0("'", gsub("'", "\\\\'", x, fixed = TRUE), "'")
}

is_nominal_vector <- function(x) {
  is.factor(x) || is.character(x)
}

nominal_predictor_names <- function(data, outcome_var) {
  predictor_names <- setdiff(names(data), outcome_var)
  predictor_names[vapply(data[predictor_names], is_nominal_vector, TRUE)]
}

numeric_predictor_names <- function(data, outcome_var) {
  predictor_names <- setdiff(names(data), outcome_var)
  predictor_names[vapply(data[predictor_names], is.numeric, TRUE)]
}

interaction_token <- function(variable, categorical_interaction_vars) {
  if (variable %in% categorical_interaction_vars) {
    return(sprintf("starts_with(%s)", single_quote_string(paste0(variable, "_"))))
  }

  backtick_name(variable)
}

interaction_terms_formula <- function(interaction_vars, categorical_interaction_vars) {
  interaction_pairs <- utils::combn(interaction_vars, 2, simplify = FALSE)
  terms <- vapply(
    interaction_pairs,
    function(pair) {
      paste(
        interaction_token(pair[1], categorical_interaction_vars),
        interaction_token(pair[2], categorical_interaction_vars),
        sep = ":"
      )
    },
    character(1)
  )

  stats::as.formula(
    sprintf("~ %s", paste(terms, collapse = " + ")),
    env = baseenv()
  )
}

validate_interaction_vars <- function(data,
                                      outcome_var,
                                      add_interactions,
                                      interaction_vars) {
  if (add_interactions == "no") {
    return(NULL)
  }

  if (is.null(interaction_vars) || !length(interaction_vars)) {
    stop("When add_interactions = 'yes', provide interaction_vars = c(...) with at least two predictors.")
  }
  if (!is.character(interaction_vars)) {
    stop("interaction_vars must be a character vector of predictor names.")
  }

  interaction_vars <- unique(interaction_vars)
  predictor_names <- setdiff(names(data), outcome_var)
  missing_vars <- setdiff(interaction_vars, predictor_names)
  if (length(missing_vars)) {
    stop(sprintf(
      "interaction_vars contains predictors not found in the data: %s",
      paste(missing_vars, collapse = ", ")
    ))
  }

  allowed_predictors <- c(
    numeric_predictor_names(data, outcome_var),
    nominal_predictor_names(data, outcome_var)
  )
  invalid_vars <- setdiff(interaction_vars, allowed_predictors)
  if (length(invalid_vars)) {
    stop(sprintf(
      "interaction_vars must contain numeric or categorical predictors only: %s",
      paste(invalid_vars, collapse = ", ")
    ))
  }
  if (length(interaction_vars) < 2) {
    stop("interaction_vars must contain at least two predictors.")
  }

  interaction_vars
}

categorical_interaction_vars <- function(data, outcome_var, interaction_vars) {
  if (is.null(interaction_vars)) {
    return(character(0))
  }

  intersect(interaction_vars, nominal_predictor_names(data, outcome_var))
}

format_interaction_setting <- function(add_interactions, interaction_vars) {
  if (add_interactions == "no") {
    return("no")
  }

  sprintf("yes (%s)", paste(interaction_vars, collapse = ", "))
}

format_correlation_setting <- function(remove_correlated_predictors, correlation_threshold) {
  if (remove_correlated_predictors == "no") {
    return("no")
  }

  sprintf("yes (threshold=%.3f)", correlation_threshold)
}

add_interaction_dummy_steps <- function(rec, add_interactions, categorical_interaction_vars) {
  if (add_interactions == "yes" && length(categorical_interaction_vars)) {
    # VSURF recipes keep raw factors. Dummy copies are added only as ingredients
    # for categorical interactions, while the original factor columns stay
    # available as predictors.
    rec <- do.call(
      recipes::step_dummy,
      c(
        list(recipe = rec),
        as.list(categorical_interaction_vars),
        list(one_hot = FALSE, keep_original_cols = TRUE)
      )
    )
  }

  rec
}

add_interaction_steps <- function(rec,
                                  add_interactions,
                                  interaction_vars,
                                  categorical_interaction_vars) {
  if (add_interactions == "yes") {
    # do.call() passes the evaluated formula into step_interact(). Passing a
    # formula object through terms = <object> works sequentially but can break
    # when caret trains recipes on PSOCK workers.
    rec <- do.call(
      recipes::step_interact,
      list(
        recipe = rec,
        terms = interaction_terms_formula(interaction_vars, categorical_interaction_vars)
      )
    )
  }
  rec
}

validate_poly_vars <- function(data,
                               outcome_var,
                               add_polynomials,
                               poly_vars,
                               poly_degree) {
  if (add_polynomials == "no") {
    return(NULL)
  }

  if (is.null(poly_vars) || !length(poly_vars)) {
    stop("When add_polynomials = 'yes', provide poly_vars = c(...) with at least one numeric predictor.")
  }
  if (!is.character(poly_vars)) {
    stop("poly_vars must be a character vector of numeric predictor names.")
  }

  poly_vars <- unique(poly_vars)
  predictor_names <- setdiff(names(data), outcome_var)
  missing_vars <- setdiff(poly_vars, predictor_names)
  if (length(missing_vars)) {
    stop(sprintf(
      "poly_vars contains predictors not found in the data: %s",
      paste(missing_vars, collapse = ", ")
    ))
  }

  invalid_vars <- setdiff(poly_vars, numeric_predictor_names(data, outcome_var))
  if (length(invalid_vars)) {
    stop(sprintf(
      "poly_vars must contain numeric predictors only: %s",
      paste(invalid_vars, collapse = ", ")
    ))
  }

  poly_degree <- as.integer(poly_degree)
  if (is.na(poly_degree) || poly_degree < 2) {
    stop("poly_degree must be an integer of 2 or greater.")
  }

  poly_vars
}

format_poly_setting <- function(add_polynomials, poly_vars, poly_degree) {
  if (add_polynomials == "no") {
    return("no")
  }

  sprintf("yes degree=%d (%s)", poly_degree, paste(poly_vars, collapse = ", "))
}

add_polynomial_steps <- function(rec, add_polynomials, poly_vars, poly_degree) {
  if (add_polynomials == "yes") {
    # step_poly() replaces each selected numeric predictor with its polynomial
    # basis columns. Interactions are already created before this step.
    rec <- do.call(
      recipes::step_poly,
      c(list(recipe = rec), as.list(poly_vars), list(degree = poly_degree))
    )
  }

  rec
}

interaction_base_terms <- function(selected_vars) {
  if (is.null(selected_vars) || !length(selected_vars)) {
    return(character(0))
  }

  unique(unlist(strsplit(selected_vars, "_x_", fixed = TRUE), use.names = FALSE))
}

selected_vars_before_interactions <- function(selected_vars, add_interactions, interaction_vars) {
  if (is.null(selected_vars)) {
    return(NULL)
  }
  if (add_interactions != "yes") {
    return(selected_vars)
  }

  # Selected interaction columns do not exist until step_interact() runs.
  # Keep their component variables first, create the interactions, then prune
  # back to the exact selected columns.
  unique(c(selected_vars, interaction_vars, interaction_base_terms(selected_vars)))
}

#################################################################BINARY OUTCOME#################################################################
#' Nested cross-validation with variable selection (VSURF and ElasticNet)
nested_elastic_binary_outcome <- function(
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
) {
  # ---- Reproducible RNG ----
  RNGkind("L'Ecuyer-CMRG")
  set.seed(seed)
  cv_method <- match.arg(inner_cv_method)
  selection_rule <- match.arg(selection_rule)
  transformation_rule <- match.arg(data_transformation)
  hyperparam_search <- match.arg(hyperparam_search)
  numeric_imputation <- match.arg(numeric_imputation)
  nominal_imputation <- match.arg(nominal_imputation)
  add_missing_indicators <- match.arg(add_missing_indicators)
  remove_zero_variance <- match.arg(remove_zero_variance)
  remove_correlated_predictors <- match.arg(remove_correlated_predictors)
  correlation_threshold <- validate_correlation_threshold(correlation_threshold)
  add_interactions <- match.arg(add_interactions)
  add_polynomials <- match.arg(add_polynomials)
  interaction_vars <- validate_interaction_vars(df, outcome_var, add_interactions, interaction_vars)
  interaction_categorical_vars <- categorical_interaction_vars(df, outcome_var, interaction_vars)
  poly_vars <- validate_poly_vars(df, outcome_var, add_polynomials, poly_vars, poly_degree)
  poly_degree <- as.integer(poly_degree)

  message(sprintf("IMPORTANT: All performance metrics and models are built around the predicted probability of the positive class: '%s'", positive_class))

  # ---- Input checks ----
  if (cv_method == "LOOCV") {
    message("Inner CV method: LOOCV (Leave-One-Out Cross-Validation)")
  } else if (cv_method == "repeatedcv") {
    message(sprintf("Inner CV method: %d-fold CV, repeated %d times", inner_cv_folds, inner_cv_repeats))
  } else {
    stop("inner_cv_method must be either 'LOOCV' or 'repeatedcv'.")
  }

  #Transformation rule
  if (transformation_rule == "YeoJohnson") {
    message("Data transformation: Yeo-Johnson")
  } else if (transformation_rule == "log") {
    message("Data transformation: Natural log-transform (0.0001 offset)")
  } else {
    stop("data_transformation must be either 'YeoJohnson' or 'log'.")
  }

  # Validate selection_rule
  if (selection_rule == "best") {
    message("Final model will use the hyperparameters with the best inner CV performance.")
  } else if (selection_rule == "oneSE") {
    message("Final model will use the most regularized hyperparameters (parsimonious model) within 1 SE of the best inner CV performance.")
  } else  {
    stop("selection_rule must be either 'best' or 'oneSE'.")
  }

  if (selection_rule == "oneSE" && cv_method == "LOOCV") {
    stop("selection_rule = 'oneSE' is not compatible with LOOCV. Choose cv_method = 'repeatedcv' instead.")
  }

  message(sprintf(
    "Hyperparameter search: %s%s",
    hyperparam_search,
    if (hyperparam_search == "random") sprintf(" (%d random candidates per glmnet model)", random_search_size) else ""
  ))
  message(sprintf(
    "Imputation: numeric=%s, nominal=%s, missing indicators=%s, remove zero variance=%s, remove correlated predictors=%s, interactions=%s, polynomials=%s",
    numeric_imputation,
    nominal_imputation,
    add_missing_indicators,
    remove_zero_variance,
    format_correlation_setting(remove_correlated_predictors, correlation_threshold),
    format_interaction_setting(add_interactions, interaction_vars),
    format_poly_setting(add_polynomials, poly_vars, poly_degree)
  ))

  # ---- Outcome factor with explicit order: negative first ----
  if (!all(c(negative_class, positive_class) %in% unique(df[[outcome_var]]))) {
    stop("Outcome does not contain both specified classes.")
  }

  # Only examine predictors; avoid evaluating < 0 on non-numerics
  pred_names <- setdiff(names(df), outcome_var)
  num_pred   <- vapply(df[pred_names], is.numeric, TRUE)
  if (transformation_rule == "log" &&
      any(vapply(df[pred_names][num_pred], function(x) any(x < 0, na.rm = TRUE), TRUE))) {
    stop("Some numeric predictors have negative values; log-transform will fail.")
  }

  if (any(vapply(df[pred_names], is.character, TRUE))) {
    message("Note: character predictors found; converting via step_string2factor().")
  }

  df[[outcome_var]] <- factor(df[[outcome_var]], levels = c(negative_class, positive_class))

  # ---- Outer resampling indices (training indices per resample) ----
  cv_outer_train_folds_rows <-
    caret::createMultiFolds(df[[outcome_var]], k = cv_outer_folds, times = cv_outer_repeats)

  # ---- Preprocessing recipe factory ----
  # VSURF keeps categorical predictors as factors.
  make_recipe_vsurf <- function(data, outcome_var,transformation_rule) {
    all_pred <- setdiff(names(data), outcome_var)
    rec <- recipes::recipe(stats::reformulate(all_pred, response = outcome_var), data = data) |>
      recipes::step_string2factor(recipes::all_nominal_predictors())

    # conditional transformation
    if (transformation_rule == "log") {
      rec <- rec |>
        recipes::step_log(recipes::all_numeric_predictors(),
                          offset = 0.0001)
    } else if (transformation_rule == "YeoJohnson") {
      rec <- rec |>
        recipes::step_YeoJohnson(recipes::all_numeric_predictors())
    }

    # Interactions are created after imputation and before scaling so products are scaled too.
    rec <- add_imputation_steps(
      rec,
      data,
      outcome_var,
      numeric_imputation,
      nominal_imputation,
      add_missing_indicators
    )
    rec <- add_interaction_dummy_steps(rec, add_interactions, interaction_categorical_vars)
    rec <- add_interaction_steps(rec, add_interactions, interaction_vars, interaction_categorical_vars)
    rec <- add_polynomial_steps(rec, add_polynomials, poly_vars, poly_degree)
    rec <- finish_recipe_preprocessing(
      rec,
      data,
      outcome_var,
      remove_zero_variance,
      remove_correlated_predictors,
      correlation_threshold
    )

    rec
  }

  # Glmnet needs a numeric model matrix, so categorical predictors are
  # dummy-expanded before elastic-net or lasso fitting.
  make_recipe_glmnet_design <- function(data, outcome_var,transformation_rule) {
    all_pred <- setdiff(names(data), outcome_var)
    rec <- recipes::recipe(stats::reformulate(all_pred, response = outcome_var), data = data) |>
      recipes::step_string2factor(recipes::all_nominal_predictors())

    # conditional transformation
    if (transformation_rule == "log") {
      rec <- rec |>
        recipes::step_log(recipes::all_numeric_predictors(),
                          offset = 0.0001)
    } else if (transformation_rule == "YeoJohnson") {
      rec <- rec |>
        recipes::step_YeoJohnson(recipes::all_numeric_predictors())
    }

    # Dummy variables are created before the requested interactions
    # so the glmnet recipe stays fully numeric.
    rec <- add_imputation_steps(
      rec,
      data,
      outcome_var,
      numeric_imputation,
      nominal_imputation,
      add_missing_indicators
    ) |>
      recipes::step_dummy(recipes::all_nominal_predictors(), one_hot = FALSE)
    rec <- add_interaction_steps(rec, add_interactions, interaction_vars, interaction_categorical_vars)
    rec <- add_polynomial_steps(rec, add_polynomials, poly_vars, poly_degree)
    rec <- finish_recipe_preprocessing(
      rec,
      data,
      outcome_var,
      remove_zero_variance,
      remove_correlated_predictors,
      correlation_threshold
    )

    rec
  }

  # Lasso is reported as its own model, so the elastic-net branch excludes alpha = 1.
  elastic_alpha_grid <- alpha_grid[alpha_grid < 1]
  if (!length(elastic_alpha_grid)) {
    stop("alpha_grid must include at least one value below 1; lasso uses alpha = 1 as a separate model.")
  }
  tg_elastic <- make_glmnet_tune_grid(
    elastic_alpha_grid,
    lambda_grid,
    hyperparam_search,
    random_search_size,
    "elastic net"
  )
  tg_lasso <- make_glmnet_tune_grid(
    1,
    lambda_grid,
    hyperparam_search,
    random_search_size,
    "lasso"
  )

  # ---- Inner CV control (optimizes logLoss) ----
  if (selection_rule == "best") {
    ctrl_inner <- caret::trainControl(
      method = cv_method,
      number = if (cv_method == "repeatedcv") inner_cv_folds else NA,
      repeats = if (cv_method == "repeatedcv") inner_cv_repeats else NA,
      savePredictions = "final",
      classProbs = TRUE,
      summaryFunction = caret::mnLogLoss,
      verboseIter = FALSE,
      allowParallel = TRUE
    )
  } else if (selection_rule == "oneSE") {
    ctrl_inner <- caret::trainControl(
      method = cv_method,
      number = if (cv_method == "repeatedcv") inner_cv_folds else NA,
      repeats = if (cv_method == "repeatedcv") inner_cv_repeats else NA,
      savePredictions = "final",
      classProbs = TRUE,
      summaryFunction = caret::mnLogLoss,
      selectionFunction = "oneSE",
      verboseIter = FALSE,
      allowParallel = TRUE
    )
  }

  n_outer <- length(cv_outer_train_folds_rows)

  sel_vars_df <- tibble::tibble(
    elastic_sel_vars = vector("list", n_outer),
    lasso_sel_vars   = vector("list", n_outer),
    vsurf_sel_vars_interp   = vector("list", n_outer),
    vsurf_sel_vars_pred   = vector("list", n_outer),
    elastic_coefs      = vector("list", n_outer),
    lasso_coefs        = vector("list", n_outer)
  )

  outer_fold_predictions <- vector("list", n_outer)
  eps <- 1e-12
  # The reference classifier is intentionally constant. This keeps its ROC AUC at 0.5.
  baseline_prevalence <- mean(df[[outcome_var]] == positive_class)
  baseline_prevalence <- pmin(pmax(baseline_prevalence, eps), 1 - eps)

  # ---- Parallel backend (once) ----
  cl <- parallel::makePSOCKcluster(max(1L, parallel::detectCores() - 1L))
  doParallel::registerDoParallel(cl)
  on.exit({
    try(parallel::stopCluster(cl), silent = TRUE)
    foreach::registerDoSEQ()
  }, add = TRUE)

  # ---- Outer CV ---------------------------------------------------------
  for (i in seq_along(cv_outer_train_folds_rows)) {
    outer_train_idx <- cv_outer_train_folds_rows[[i]]
    outer_test_idx  <- setdiff(seq_len(nrow(df)), outer_train_idx)
    outer_d_train   <- df[outer_train_idx, , drop = FALSE]
    outer_d_test    <- df[outer_test_idx,  , drop = FALSE]

    tab <- table(outer_d_train[[outcome_var]])
    message(sprintf("OUTER %d - class counts in the training fold: %s",
                    i, paste(sprintf("%s=%d", names(tab), tab), collapse=", ")))

    # Keep negative first consistently
    outer_d_train[[outcome_var]] <- factor(outer_d_train[[outcome_var]],
                                           levels = c(negative_class, positive_class))
    outer_d_test[[outcome_var]]  <- factor(outer_d_test[[outcome_var]],
                                           levels = c(negative_class, positive_class))

    # Bake training data for VSURF
    rec_vsurf_full <- make_recipe_vsurf(outer_d_train, outcome_var,transformation_rule)
    prep_rec_vsurf_full <- recipes::prep(rec_vsurf_full, training = outer_d_train, retain = TRUE)
    d_train_baked <- recipes::bake(prep_rec_vsurf_full, new_data = NULL)
    d_test_baked  <- recipes::bake(prep_rec_vsurf_full, new_data = outer_d_test)

    x_train <- d_train_baked |> dplyr::select(-dplyr::all_of(outcome_var)) |> as.data.frame()
    y_train <- d_train_baked[[outcome_var]]
    x_test  <- d_test_baked  |> dplyr::select(-dplyr::all_of(outcome_var)) |> as.data.frame()


    # ---- VSURF ----
    vs <- VSURF::VSURF(
      y = y_train,
      x = x_train,
      ntree.thres   = ntree, nfor.thres  = nforests,
      ntree.interp  = ntree, nfor.interp = nforests,
      ntree.pred    = ntree, nfor.pred   = nforests,
      RFimplem = "ranger",
      parallel = FALSE
    )

    sel_idx_inter <- vs$varselect.interp
    sel_idx_pred <- vs$varselect.pred
    vsurf_sel_vars_in <- if (length(sel_idx_inter)) colnames(x_train)[sel_idx_inter] else character(0)
    vsurf_sel_vars_pr <- if (length(sel_idx_pred)) colnames(x_train)[sel_idx_pred] else character(0)
    sel_vars_df$vsurf_sel_vars_interp[[i]] <- vsurf_sel_vars_in
    sel_vars_df$vsurf_sel_vars_pred[[i]] <- vsurf_sel_vars_pr
    message(sprintf("OUTER %02d - VSURF selected at the interpretation step (%d): %s",
                    i, length(vsurf_sel_vars_in), paste(vsurf_sel_vars_in, collapse=", ")))
    message(sprintf("OUTER %02d - VSURF selected at the prediction step (%d): %s",
                    i, length(vsurf_sel_vars_pr), paste(vsurf_sel_vars_pr, collapse=", ")))

    # Keep the original training data inside the VSURF object so predict()
    # remains self-contained when called later in the loop.
    vs$x_train_saved <- x_train
    vs$y_train_saved <- y_train
    vs$call$x <- quote(object$x_train_saved)
    vs$call$y <- quote(object$y_train_saved)

    # probabilities at interpretation step (always available if vsurf ran)
    p_vsurf_inter <- predict(vs, newdata = x_test, step = "interp", type = "prob")[, positive_class]    #predictions at the interpretation step

    # prediction step may be unavailable
    if (is.null(vs$varselect.pred) || length(vs$varselect.pred) == 0) {
      # fallback: use interpretation step as "pred"
      p_vsurf_pred <- p_vsurf_inter
      # or: p_vsurf_pred <- rep(NA_real_, nrow(x_test))
    } else {
      p_vsurf_pred <- predict(vs, newdata = x_test, step = "pred", type = "prob")[, positive_class]     #predictions at the prediction step
    }


    p_vsurf_inter <- pmin(pmax(p_vsurf_inter, eps), 1 - eps)
    p_vsurf_pred  <- pmin(pmax(p_vsurf_pred,  eps), 1 - eps)

    # ---- Elastic-net fit and coefficient selection on the glmnet design ----
    rec_glmnet <- make_recipe_glmnet_design(outer_d_train, outcome_var,transformation_rule)

    elastic_sel <- caret::train(
      rec_glmnet,
      data = outer_d_train,
      method = "glmnet",
      metric = "logLoss",
      maximize = FALSE,
      tuneGrid = tg_elastic,
      trControl = ctrl_inner,
      family = "binomial",
      standardize = FALSE
    )
    coefs <- as.matrix(coef(elastic_sel$finalModel,
                            s = elastic_sel$bestTune$lambda[1]))[, 1]
    el_sel_vars <- setdiff(names(coefs)[coefs != 0], "(Intercept)")
    sel_vars_df$elastic_sel_vars[[i]] <- el_sel_vars
    sel_vars_df$elastic_coefs[[i]] <- coefs
    message(sprintf("OUTER %02d - Elastic selected (%d): %s",
                    i, length(el_sel_vars), paste(el_sel_vars, collapse=", ")))

    # Predictions from the tuned elastic-net model for this outer training fold.
    p_elas <- predict(elastic_sel, newdata = outer_d_test, type = "prob")[[positive_class]]
    p_elas <- pmin(pmax(p_elas, eps), 1 - eps)

    # ---- Lasso fit and coefficient selection on the glmnet design ----
    lasso_sel <- caret::train(
      rec_glmnet,
      data = outer_d_train,
      method = "glmnet",
      metric = "logLoss",
      maximize = FALSE,
      tuneGrid = tg_lasso,
      trControl = ctrl_inner,
      family = "binomial",
      standardize = FALSE
    )
    lasso_coefs <- as.matrix(coef(lasso_sel$finalModel,
                                  s = lasso_sel$bestTune$lambda[1]))[, 1]
    lasso_sel_vars <- setdiff(names(lasso_coefs)[lasso_coefs != 0], "(Intercept)")
    sel_vars_df$lasso_sel_vars[[i]] <- lasso_sel_vars
    sel_vars_df$lasso_coefs[[i]] <- lasso_coefs
    message(sprintf("OUTER %02d - Lasso selected (%d): %s",
                    i, length(lasso_sel_vars), paste(lasso_sel_vars, collapse=", ")))

    p_lasso <- predict(lasso_sel, newdata = outer_d_test, type = "prob")[[positive_class]]
    p_lasso <- pmin(pmax(p_lasso, eps), 1 - eps)

    baseline_pred <- rep(baseline_prevalence, nrow(outer_d_test))

    # Store per-row predictions (align to indices)
    outer_fold_predictions[[i]] <- list(
      row_id = outer_test_idx,
      outcome = outer_d_test[[outcome_var]],
      pred_elastic_net = p_elas,
      pred_lasso = p_lasso,
      pred_prevalence = baseline_pred,
      pred_vsurf_interpretation = p_vsurf_inter,
      pred_vsurf_prediction = p_vsurf_pred
    )

  } # End of outer loop
  # ---- Aggregate predictions over outer folds ----
  outer_predictions_long <- tibble::tibble()

  for (i in seq_along(outer_fold_predictions)) {
    outer_predictions_long <- dplyr::bind_rows(
      outer_predictions_long,
      tibble::tibble(
        row_id = outer_fold_predictions[[i]]$row_id,
        pred_elastic_net = outer_fold_predictions[[i]]$pred_elastic_net,
        pred_lasso = outer_fold_predictions[[i]]$pred_lasso,
        outcome = outer_fold_predictions[[i]]$outcome,
        pred_prevalence = outer_fold_predictions[[i]]$pred_prevalence,
        pred_vsurf_interpretation = outer_fold_predictions[[i]]$pred_vsurf_interpretation,
        pred_vsurf_prediction = outer_fold_predictions[[i]]$pred_vsurf_prediction
      )
    )
  }

  # Average repeated outer-CV predictions to one row per original observation.
  outer_predictions <- outer_predictions_long |>
    dplyr::group_by(row_id) |>
    dplyr::summarise(
      outcome = dplyr::first(as.character(outcome)),
      pred_elastic_net = mean(pred_elastic_net, na.rm = TRUE),
      pred_lasso = mean(pred_lasso, na.rm = TRUE),
      pred_prevalence = mean(pred_prevalence, na.rm = TRUE),
      pred_vsurf_interpretation = mean(pred_vsurf_interpretation, na.rm = TRUE),
      pred_vsurf_prediction = mean(pred_vsurf_prediction, na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::arrange(row_id)

  # Return one prediction table instead of parallel vectors with repeated ids/outcomes.
  return(list(
    outer_predictions = outer_predictions,
    prevalence = baseline_prevalence,
    selected_variables = sel_vars_df,
    tuning = list(
      search = hyperparam_search,
      random_search_size = if (hyperparam_search == "random") random_search_size else NA_integer_,
      elastic_net_candidates = tg_elastic,
      lasso_candidates = tg_lasso
    ),
    preprocessing = list(
      numeric_imputation = numeric_imputation,
      nominal_imputation = nominal_imputation,
      add_missing_indicators = add_missing_indicators,
      remove_zero_variance = remove_zero_variance,
      remove_correlated_predictors = remove_correlated_predictors,
      correlation_threshold = correlation_threshold,
      add_interactions = add_interactions,
      interaction_vars = interaction_vars,
      add_polynomials = add_polynomials,
      poly_vars = poly_vars,
      poly_degree = poly_degree
    )
  ))
}

# ---- Outer performance on averaged predictions ----
outer_perf_nested_binary <- function(trained_object,
                                           positive_class = "Yes",
                                           negative_class = "No") {

  message(sprintf("IMPORTANT: All performance metrics are built around the predicted probability of the positive class: '%s'", positive_class))

  preds <- trained_object$outer_predictions
  selected_variables <- trained_object$selected_variables

  # ---- Variable selection stability ----
  stab_table <- function(sel_list) {
    selected <- unlist(sel_list)
    if (!length(selected)) {
      return(tibble::tibble(variable = character(), times_selected = integer(), freq = numeric()))
    }
    tbl <- sort(table(selected), decreasing = TRUE)
    tibble::tibble(variable = names(tbl), times_selected = as.integer(tbl)) |>
      dplyr::mutate(freq = times_selected / nrow(selected_variables))
  }

  cat("\n--- Elastic Net: selection stability ---\n")
  print(stab_table(selected_variables$elastic_sel_vars))
  cat("\n--- Lasso: selection stability ---\n")
  print(stab_table(selected_variables$lasso_sel_vars))
  cat("\n--- VSURF: selection stability at the interpretation step ---\n")
  print(stab_table(selected_variables$vsurf_sel_vars_interp))
  cat("\n--- VSURF: selection stability at the prediction step ---\n")
  print(stab_table(selected_variables$vsurf_sel_vars_pred))

  get_metrics <- function(y_true, y_pred) {
    # set factor with negative first, positive second (matches event_level='second')
    truth <- factor(y_true, levels = c(negative_class, positive_class))

    # clamp to avoid Inf in qlogis and to be safe for metrics
    eps <- 1e-12
    estimate <- pmin(pmax(y_pred, eps), 1 - eps)

    truth_num <- as.integer(truth == positive_class)
    hard_pred <- factor(ifelse(estimate >= 0.5, positive_class, negative_class),
                        levels = levels(truth))

    roc_auc    <- yardstick::roc_auc_vec(truth, estimate, event_level = "second")
    pr_auc     <- yardstick::pr_auc_vec (truth, estimate, event_level = "second")
    prg_curve  <- prg::create_prg_curve(truth_num, estimate)
    auprg      <- prg::calc_auprg(prg_curve)
    logloss    <- yardstick::mn_log_loss_vec(truth, estimate, event_level = "second")
    brier      <- ModelMetrics::brier(actual = truth_num, predicted = estimate)
    mcc        <- yardstick::mcc_vec(truth, hard_pred, event_level = "second")

    logit_p <- qlogis(estimate)
    fit_citl <- stats::glm(truth_num ~ 1 + offset(logit_p), family = binomial())
    citl <- unname(coef(fit_citl)[1])

    if (stats::sd(logit_p, na.rm = TRUE) == 0) {
      cal_int <- NA_real_
      cal_slope <- NA_real_
    } else {
      fit_cal <- stats::glm(truth_num ~ logit_p, family = binomial())
      cal_int <- unname(coef(fit_cal)[1])
      cal_slope <- unname(coef(fit_cal)[2])
    }

    tibble::tibble(
      roc_auc = roc_auc,
      pr_auc  = pr_auc,
      auprg   = auprg,
      logloss = logloss,
      brier   = brier,
      mcc     = mcc,
      cal_intercept = cal_int,
      cal_slope     = cal_slope,
      citl_intercept= citl
    )
  }

  outer_single_elas_metrics <- get_metrics(preds$outcome, preds$pred_elastic_net)
  lasso_metrics <- get_metrics(preds$outcome, preds$pred_lasso)
  baseline_mod  <- get_metrics(preds$outcome, preds$pred_prevalence)
  vsurf_mod_inter  <- get_metrics(preds$outcome, preds$pred_vsurf_interpretation)
  vsurf_mod_pred  <- get_metrics(preds$outcome, preds$pred_vsurf_prediction)

  single_elas_summary <- dplyr::as_tibble(outer_single_elas_metrics) |> dplyr::mutate(method = "Elastic net (outer)")
  lasso_summary <- dplyr::as_tibble(lasso_metrics) |> dplyr::mutate(method = "Lasso (outer)")
  vsurf_inter_summary <- dplyr::as_tibble(vsurf_mod_inter) |> dplyr::mutate(method = "VSURF interpretation step (outer)")
  vsurf_pred_summary <- dplyr::as_tibble(vsurf_mod_pred) |> dplyr::mutate(method = "VSURF prediction step (outer)")
  base_summary <- dplyr::as_tibble(baseline_mod) |> dplyr::mutate(method = "Reference prevalence model (constant)")

  dplyr::bind_rows(single_elas_summary, lasso_summary, vsurf_inter_summary, vsurf_pred_summary, base_summary) |>
    dplyr::select(method, dplyr::everything())
}

#################################################################MULTICLASS OUTCOME#################################################################

nested_elastic_multiclass_outcome <- function(
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
) {
  RNGkind("L'Ecuyer-CMRG")
  set.seed(seed)
  cv_method <- match.arg(inner_cv_method)
  selection_rule <- match.arg(selection_rule)
  transformation_rule <- match.arg(data_transformation)
  hyperparam_search <- match.arg(hyperparam_search)
  numeric_imputation <- match.arg(numeric_imputation)
  nominal_imputation <- match.arg(nominal_imputation)
  add_missing_indicators <- match.arg(add_missing_indicators)
  remove_zero_variance <- match.arg(remove_zero_variance)
  remove_correlated_predictors <- match.arg(remove_correlated_predictors)
  correlation_threshold <- validate_correlation_threshold(correlation_threshold)
  add_interactions <- match.arg(add_interactions)
  add_polynomials <- match.arg(add_polynomials)
  interaction_vars <- validate_interaction_vars(df, outcome_var, add_interactions, interaction_vars)
  interaction_categorical_vars <- categorical_interaction_vars(df, outcome_var, interaction_vars)
  poly_vars <- validate_poly_vars(df, outcome_var, add_polynomials, poly_vars, poly_degree)
  poly_degree <- as.integer(poly_degree)

  df[[outcome_var]] <- factor(df[[outcome_var]])
  class_levels <- levels(df[[outcome_var]])
  if (length(class_levels) < 3) {
    stop("Multiclass outcome requires at least three classes. Use nested_elastic_binary_outcome() for two classes.")
  }

  pred_names <- setdiff(names(df), outcome_var)
  num_pred   <- vapply(df[pred_names], is.numeric, TRUE)
  if (transformation_rule == "log" &&
      any(vapply(df[pred_names][num_pred], function(x) any(x < 0, na.rm = TRUE), TRUE))) {
    stop("Some numeric predictors have negative values; log-transform will fail.")
  }

  if (selection_rule == "oneSE" && cv_method == "LOOCV") {
    stop("selection_rule = 'oneSE' is not compatible with LOOCV. Choose cv_method = 'repeatedcv' instead.")
  }

  message(sprintf("IMPORTANT: Multiclass metrics use the class probability columns for: %s", paste(class_levels, collapse = ", ")))
  message(sprintf(
    "Inner CV method: %s",
    if (cv_method == "LOOCV") "LOOCV (Leave-One-Out Cross-Validation)" else sprintf("%d-fold CV, repeated %d times", inner_cv_folds, inner_cv_repeats)
  ))
  message(sprintf(
    "Hyperparameter search: %s%s",
    hyperparam_search,
    if (hyperparam_search == "random") sprintf(" (%d random candidates per glmnet model)", random_search_size) else ""
  ))
  message(sprintf(
    "Imputation: numeric=%s, nominal=%s, missing indicators=%s, remove zero variance=%s, remove correlated predictors=%s, interactions=%s, polynomials=%s",
    numeric_imputation,
    nominal_imputation,
    add_missing_indicators,
    remove_zero_variance,
    format_correlation_setting(remove_correlated_predictors, correlation_threshold),
    format_interaction_setting(add_interactions, interaction_vars),
    format_poly_setting(add_polynomials, poly_vars, poly_degree)
  ))

  make_recipe_vsurf <- function(data, outcome_var, transformation_rule) {
    all_pred <- setdiff(names(data), outcome_var)
    rec <- recipes::recipe(stats::reformulate(all_pred, response = outcome_var), data = data) |>
      recipes::step_string2factor(recipes::all_nominal_predictors())

    if (transformation_rule == "log") {
      rec <- rec |>
        recipes::step_log(recipes::all_numeric_predictors(), offset = 0.0001)
    } else if (transformation_rule == "YeoJohnson") {
      rec <- rec |>
        recipes::step_YeoJohnson(recipes::all_numeric_predictors())
    }

    rec <- add_imputation_steps(
      rec,
      data,
      outcome_var,
      numeric_imputation,
      nominal_imputation,
      add_missing_indicators
    )
    rec <- add_interaction_dummy_steps(rec, add_interactions, interaction_categorical_vars)
    rec <- add_interaction_steps(rec, add_interactions, interaction_vars, interaction_categorical_vars)
    rec <- add_polynomial_steps(rec, add_polynomials, poly_vars, poly_degree)
    finish_recipe_preprocessing(
      rec,
      data,
      outcome_var,
      remove_zero_variance,
      remove_correlated_predictors,
      correlation_threshold
    )
  }

  make_recipe_glmnet_design <- function(data, outcome_var, transformation_rule) {
    all_pred <- setdiff(names(data), outcome_var)
    rec <- recipes::recipe(stats::reformulate(all_pred, response = outcome_var), data = data) |>
      recipes::step_string2factor(recipes::all_nominal_predictors())

    if (transformation_rule == "log") {
      rec <- rec |>
        recipes::step_log(recipes::all_numeric_predictors(), offset = 0.0001)
    } else if (transformation_rule == "YeoJohnson") {
      rec <- rec |>
        recipes::step_YeoJohnson(recipes::all_numeric_predictors())
    }

    rec <- add_imputation_steps(
      rec,
      data,
      outcome_var,
      numeric_imputation,
      nominal_imputation,
      add_missing_indicators
    ) |>
      recipes::step_dummy(recipes::all_nominal_predictors(), one_hot = FALSE)
    rec <- add_interaction_steps(rec, add_interactions, interaction_vars, interaction_categorical_vars)
    rec <- add_polynomial_steps(rec, add_polynomials, poly_vars, poly_degree)
    finish_recipe_preprocessing(
      rec,
      data,
      outcome_var,
      remove_zero_variance,
      remove_correlated_predictors,
      correlation_threshold
    )
  }

  normalize_prob_matrix <- function(prob) {
    prob <- as.data.frame(prob)
    missing_classes <- setdiff(class_levels, names(prob))
    for (class_name in missing_classes) {
      prob[[class_name]] <- 0
    }
    prob <- as.matrix(prob[, class_levels, drop = FALSE])
    prob <- pmin(pmax(prob, 1e-12), 1 - 1e-12)
    prob / rowSums(prob)
  }

  extract_multiclass_glmnet_selection <- function(fit) {
    coefs_by_class <- coef(fit$finalModel, s = fit$bestTune$lambda[1])
    if (!is.list(coefs_by_class)) {
      coefs_by_class <- list(model = coefs_by_class)
    }
    selected <- unique(unlist(lapply(coefs_by_class, function(coefs) {
      coef_vec <- as.matrix(coefs)[, 1]
      setdiff(names(coef_vec)[coef_vec != 0], "(Intercept)")
    })))
    list(selected = selected, coefs = coefs_by_class)
  }

  elastic_alpha_grid <- alpha_grid[alpha_grid < 1]
  if (!length(elastic_alpha_grid)) {
    stop("alpha_grid must include at least one value below 1; lasso uses alpha = 1 as a separate model.")
  }
  tg_elastic <- make_glmnet_tune_grid(
    elastic_alpha_grid,
    lambda_grid,
    hyperparam_search,
    random_search_size,
    "elastic net"
  )
  tg_lasso <- make_glmnet_tune_grid(
    1,
    lambda_grid,
    hyperparam_search,
    random_search_size,
    "lasso"
  )

  if (selection_rule == "best") {
    ctrl_inner <- caret::trainControl(
      method = cv_method,
      number = if (cv_method == "repeatedcv") inner_cv_folds else NA,
      repeats = if (cv_method == "repeatedcv") inner_cv_repeats else NA,
      savePredictions = "final",
      classProbs = TRUE,
      summaryFunction = caret::mnLogLoss,
      verboseIter = FALSE,
      allowParallel = TRUE
    )
  } else {
    ctrl_inner <- caret::trainControl(
      method = cv_method,
      number = if (cv_method == "repeatedcv") inner_cv_folds else NA,
      repeats = if (cv_method == "repeatedcv") inner_cv_repeats else NA,
      savePredictions = "final",
      classProbs = TRUE,
      summaryFunction = caret::mnLogLoss,
      selectionFunction = "oneSE",
      verboseIter = FALSE,
      allowParallel = TRUE
    )
  }

  cv_outer_train_folds_rows <- caret::createMultiFolds(df[[outcome_var]], k = cv_outer_folds, times = cv_outer_repeats)
  n_outer <- length(cv_outer_train_folds_rows)

  sel_vars_df <- tibble::tibble(
    elastic_sel_vars = vector("list", n_outer),
    lasso_sel_vars   = vector("list", n_outer),
    vsurf_sel_vars_interp = vector("list", n_outer),
    vsurf_sel_vars_pred = vector("list", n_outer),
    elastic_coefs = vector("list", n_outer),
    lasso_coefs = vector("list", n_outer)
  )
  outer_fold_predictions <- vector("list", n_outer)

  baseline_prevalence <- as.numeric(prop.table(table(df[[outcome_var]]))[class_levels])
  names(baseline_prevalence) <- class_levels
  baseline_prevalence <- pmin(pmax(baseline_prevalence, 1e-12), 1 - 1e-12)
  baseline_prevalence <- baseline_prevalence / sum(baseline_prevalence)

  cl <- parallel::makePSOCKcluster(max(1L, parallel::detectCores() - 1L))
  doParallel::registerDoParallel(cl)
  on.exit({
    try(parallel::stopCluster(cl), silent = TRUE)
    foreach::registerDoSEQ()
  }, add = TRUE)

  for (i in seq_along(cv_outer_train_folds_rows)) {
    outer_train_idx <- cv_outer_train_folds_rows[[i]]
    outer_test_idx  <- setdiff(seq_len(nrow(df)), outer_train_idx)
    outer_d_train   <- df[outer_train_idx, , drop = FALSE]
    outer_d_test    <- df[outer_test_idx,  , drop = FALSE]
    outer_d_train[[outcome_var]] <- factor(outer_d_train[[outcome_var]], levels = class_levels)
    outer_d_test[[outcome_var]] <- factor(outer_d_test[[outcome_var]], levels = class_levels)

    tab <- table(outer_d_train[[outcome_var]])
    message(sprintf("OUTER %d - class counts in the training fold: %s",
                    i, paste(sprintf("%s=%d", names(tab), tab), collapse = ", ")))

    rec_vsurf <- make_recipe_vsurf(outer_d_train, outcome_var, transformation_rule)
    prep_rec_vsurf <- recipes::prep(rec_vsurf, training = outer_d_train, retain = TRUE)
    d_train_baked <- recipes::bake(prep_rec_vsurf, new_data = NULL)
    d_test_baked  <- recipes::bake(prep_rec_vsurf, new_data = outer_d_test)

    x_train <- d_train_baked |> dplyr::select(-dplyr::all_of(outcome_var)) |> as.data.frame()
    y_train <- factor(d_train_baked[[outcome_var]], levels = class_levels)
    x_test  <- d_test_baked |> dplyr::select(-dplyr::all_of(outcome_var)) |> as.data.frame()

    vs <- VSURF::VSURF(
      y = y_train,
      x = x_train,
      ntree.thres = ntree, nfor.thres = nforests,
      ntree.interp = ntree, nfor.interp = nforests,
      ntree.pred = ntree, nfor.pred = nforests,
      RFimplem = "ranger",
      parallel = FALSE
    )

    sel_idx_inter <- vs$varselect.interp
    sel_idx_pred <- vs$varselect.pred
    vsurf_sel_vars_in <- if (length(sel_idx_inter)) colnames(x_train)[sel_idx_inter] else character(0)
    vsurf_sel_vars_pr <- if (length(sel_idx_pred)) colnames(x_train)[sel_idx_pred] else character(0)
    sel_vars_df$vsurf_sel_vars_interp[[i]] <- vsurf_sel_vars_in
    sel_vars_df$vsurf_sel_vars_pred[[i]] <- vsurf_sel_vars_pr

    vs$x_train_saved <- x_train
    vs$y_train_saved <- y_train
    vs$call$x <- quote(object$x_train_saved)
    vs$call$y <- quote(object$y_train_saved)

    p_vsurf_inter <- normalize_prob_matrix(predict(vs, newdata = x_test, step = "interp", type = "prob"))
    p_vsurf_pred <- if (is.null(vs$varselect.pred) || length(vs$varselect.pred) == 0) {
      p_vsurf_inter
    } else {
      normalize_prob_matrix(predict(vs, newdata = x_test, step = "pred", type = "prob"))
    }

    rec_glmnet <- make_recipe_glmnet_design(outer_d_train, outcome_var, transformation_rule)
    elastic_sel <- caret::train(
      rec_glmnet,
      data = outer_d_train,
      method = "glmnet",
      metric = "logLoss",
      maximize = FALSE,
      tuneGrid = tg_elastic,
      trControl = ctrl_inner,
      family = "multinomial",
      standardize = FALSE
    )
    elastic_info <- extract_multiclass_glmnet_selection(elastic_sel)
    sel_vars_df$elastic_sel_vars[[i]] <- elastic_info$selected
    sel_vars_df$elastic_coefs[[i]] <- elastic_info$coefs
    p_elas <- normalize_prob_matrix(predict(elastic_sel, newdata = outer_d_test, type = "prob"))

    lasso_sel <- caret::train(
      rec_glmnet,
      data = outer_d_train,
      method = "glmnet",
      metric = "logLoss",
      maximize = FALSE,
      tuneGrid = tg_lasso,
      trControl = ctrl_inner,
      family = "multinomial",
      standardize = FALSE
    )
    lasso_info <- extract_multiclass_glmnet_selection(lasso_sel)
    sel_vars_df$lasso_sel_vars[[i]] <- lasso_info$selected
    sel_vars_df$lasso_coefs[[i]] <- lasso_info$coefs
    p_lasso <- normalize_prob_matrix(predict(lasso_sel, newdata = outer_d_test, type = "prob"))

    n_test <- nrow(outer_d_test)
    outer_fold_predictions[[i]] <- tibble::tibble(
      row_id = rep(outer_test_idx, each = length(class_levels)),
      outcome = rep(as.character(outer_d_test[[outcome_var]]), each = length(class_levels)),
      class = rep(class_levels, times = n_test),
      pred_elastic_net = as.vector(t(p_elas)),
      pred_lasso = as.vector(t(p_lasso)),
      pred_prevalence = rep(baseline_prevalence[class_levels], times = n_test),
      pred_vsurf_interpretation = as.vector(t(p_vsurf_inter)),
      pred_vsurf_prediction = as.vector(t(p_vsurf_pred))
    )
  }

  outer_predictions <- dplyr::bind_rows(outer_fold_predictions) |>
    dplyr::group_by(row_id, outcome, class) |>
    dplyr::summarise(
      pred_elastic_net = mean(pred_elastic_net, na.rm = TRUE),
      pred_lasso = mean(pred_lasso, na.rm = TRUE),
      pred_prevalence = mean(pred_prevalence, na.rm = TRUE),
      pred_vsurf_interpretation = mean(pred_vsurf_interpretation, na.rm = TRUE),
      pred_vsurf_prediction = mean(pred_vsurf_prediction, na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::arrange(row_id, factor(class, levels = class_levels))

  return(list(
    outer_predictions = outer_predictions,
    classes = class_levels,
    prevalence = baseline_prevalence,
    selected_variables = sel_vars_df,
    tuning = list(
      search = hyperparam_search,
      random_search_size = if (hyperparam_search == "random") random_search_size else NA_integer_,
      elastic_net_candidates = tg_elastic,
      lasso_candidates = tg_lasso
    ),
    preprocessing = list(
      numeric_imputation = numeric_imputation,
      nominal_imputation = nominal_imputation,
      add_missing_indicators = add_missing_indicators,
      remove_zero_variance = remove_zero_variance,
      remove_correlated_predictors = remove_correlated_predictors,
      correlation_threshold = correlation_threshold,
      add_interactions = add_interactions,
      interaction_vars = interaction_vars,
      add_polynomials = add_polynomials,
      poly_vars = poly_vars,
      poly_degree = poly_degree
    )
  ))
}

outer_perf_nested_multiclass <- function(trained_object) {
  preds <- trained_object$outer_predictions
  class_levels <- trained_object$classes
  selected_variables <- trained_object$selected_variables

  stab_table <- function(sel_list) {
    selected <- unlist(sel_list)
    if (!length(selected)) {
      return(tibble::tibble(variable = character(), times_selected = integer(), freq = numeric()))
    }
    tbl <- sort(table(selected), decreasing = TRUE)
    tibble::tibble(variable = names(tbl), times_selected = as.integer(tbl)) |>
      dplyr::mutate(freq = times_selected / nrow(selected_variables))
  }

  cat("\n--- Elastic Net: selection stability ---\n")
  print(stab_table(selected_variables$elastic_sel_vars))
  cat("\n--- Lasso: selection stability ---\n")
  print(stab_table(selected_variables$lasso_sel_vars))
  cat("\n--- VSURF: selection stability at the interpretation step ---\n")
  print(stab_table(selected_variables$vsurf_sel_vars_interp))
  cat("\n--- VSURF: selection stability at the prediction step ---\n")
  print(stab_table(selected_variables$vsurf_sel_vars_pred))

  class_prob_cols <- paste0(".pred_", make.names(class_levels))
  names(class_prob_cols) <- class_levels

  get_metrics <- function(method_col) {
    prob_wide <- preds |>
      dplyr::mutate(prob_col = unname(class_prob_cols[class])) |>
      dplyr::select(row_id, outcome, prob_col, prob = dplyr::all_of(method_col)) |>
      tidyr::pivot_wider(names_from = prob_col, values_from = prob) |>
      dplyr::arrange(row_id)

    for (prob_col in unname(class_prob_cols)) {
      if (!prob_col %in% names(prob_wide)) {
        prob_wide[[prob_col]] <- NA_real_
      }
    }

    truth <- factor(prob_wide$outcome, levels = class_levels)
    prob_mat <- as.matrix(prob_wide[, unname(class_prob_cols), drop = FALSE])
    prob_mat <- pmin(pmax(prob_mat, 1e-12), 1 - 1e-12)
    prob_mat <- prob_mat / rowSums(prob_mat)
    hard_pred <- factor(class_levels[max.col(prob_mat, ties.method = "first")], levels = class_levels)
    truth_idx <- match(as.character(truth), class_levels)
    logloss <- -mean(log(prob_mat[cbind(seq_along(truth_idx), truth_idx)]), na.rm = TRUE)

    prob_wide[, unname(class_prob_cols)] <- prob_mat
    prob_wide$outcome <- truth

    roc_auc_hand_till <- tryCatch(
      yardstick::roc_auc(
        prob_wide,
        outcome,
        dplyr::all_of(unname(class_prob_cols)),
        estimator = "hand_till"
      ) |> dplyr::pull(.estimate),
      error = function(e) NA_real_
    )

    tibble::tibble(
      accuracy = yardstick::accuracy_vec(truth, hard_pred),
      kap = yardstick::kap_vec(truth, hard_pred),
      mn_log_loss = logloss,
      roc_auc_hand_till = roc_auc_hand_till,
      macro_f_meas = yardstick::f_meas_vec(truth, hard_pred, estimator = "macro"),
      macro_sensitivity = yardstick::sens_vec(truth, hard_pred, estimator = "macro"),
      macro_precision = yardstick::precision_vec(truth, hard_pred, estimator = "macro")
    )
  }

  dplyr::bind_rows(
    get_metrics("pred_elastic_net") |> dplyr::mutate(method = "Elastic net (outer)"),
    get_metrics("pred_lasso") |> dplyr::mutate(method = "Lasso (outer)"),
    get_metrics("pred_vsurf_interpretation") |> dplyr::mutate(method = "VSURF interpretation step (outer)"),
    get_metrics("pred_vsurf_prediction") |> dplyr::mutate(method = "VSURF prediction step (outer)"),
    get_metrics("pred_prevalence") |> dplyr::mutate(method = "Reference prevalence model (constant)")
  ) |>
    dplyr::select(method, dplyr::everything())
}

#################################################################CONTINUOUS OUTCOME#################################################################

nested_elastic_continuous_outcome <- function(
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
) {
  # ---- Reproducible RNG ----
  RNGkind("L'Ecuyer-CMRG")
  set.seed(seed)
  cv_method <- match.arg(inner_cv_method)
  selection_rule <- match.arg(selection_rule)
  metric <- match.arg(optim_metric)
  transformation_rule <- match.arg(data_transformation)
  hyperparam_search <- match.arg(hyperparam_search)
  numeric_imputation <- match.arg(numeric_imputation)
  nominal_imputation <- match.arg(nominal_imputation)
  add_missing_indicators <- match.arg(add_missing_indicators)
  remove_zero_variance <- match.arg(remove_zero_variance)
  remove_correlated_predictors <- match.arg(remove_correlated_predictors)
  correlation_threshold <- validate_correlation_threshold(correlation_threshold)
  add_interactions <- match.arg(add_interactions)
  add_polynomials <- match.arg(add_polynomials)
  interaction_vars <- validate_interaction_vars(df, outcome_var, add_interactions, interaction_vars)
  interaction_categorical_vars <- categorical_interaction_vars(df, outcome_var, interaction_vars)
  poly_vars <- validate_poly_vars(df, outcome_var, add_polynomials, poly_vars, poly_degree)
  poly_degree <- as.integer(poly_degree)


  # ---- Pre-checks ----
  pred_names <- setdiff(names(df), outcome_var)
  num_pred   <- vapply(df[pred_names], is.numeric, TRUE)

  if (transformation_rule == "log" &&
      any(vapply(df[pred_names][num_pred], function(x) any(x < 0, na.rm = TRUE), TRUE))) {
    stop("Some numeric predictors have negative values; log-transform will fail.")
  }

  if (any(vapply(df[pred_names], is.character, TRUE))) {
    message("Note: character predictors found; converting via step_string2factor().")
  }

  # ---- Input checks ----
  if (cv_method == "LOOCV") {
    message("Inner CV method: LOOCV (Leave-One-Out Cross-Validation)")
  } else if (cv_method == "repeatedcv") {
    message(sprintf("Inner CV method: %d-fold CV, repeated %d times", inner_cv_folds, inner_cv_repeats))
  } else {
    stop("inner_cv_method must be either 'LOOCV' or 'repeatedcv'.")
  }

  # Validate selection_rule
  if (selection_rule == "best") {
    message("Final model will use the hyperparameters with the best inner CV performance.")
  } else if (selection_rule == "oneSE") {
    message("Final model will use the most regularized hyperparameters (parsimonious model) within 1 SE of the best inner CV performance.")
  } else  {
    stop("selection_rule must be either 'best' or 'oneSE'.")
  }

  if (selection_rule == "oneSE" && cv_method == "LOOCV") {
    stop("selection_rule = 'oneSE' is not compatible with LOOCV. Choose cv_method = 'repeatedcv' instead.")
  }

  message(sprintf(
    "Hyperparameter search: %s%s",
    hyperparam_search,
    if (hyperparam_search == "random") sprintf(" (%d random candidates per glmnet model)", random_search_size) else ""
  ))
  message(sprintf(
    "Imputation: numeric=%s, nominal=%s, missing indicators=%s, remove zero variance=%s, remove correlated predictors=%s, interactions=%s, polynomials=%s",
    numeric_imputation,
    nominal_imputation,
    add_missing_indicators,
    remove_zero_variance,
    format_correlation_setting(remove_correlated_predictors, correlation_threshold),
    format_interaction_setting(add_interactions, interaction_vars),
    format_poly_setting(add_polynomials, poly_vars, poly_degree)
  ))

  #Transformation rule
  if (transformation_rule == "YeoJohnson") {
    message("Data transformation: Yeo-Johnson")
  } else if (transformation_rule == "log") {
    message("Data transformation: Natural log-transform (0.0001 offset)")
  } else {
    stop("data_transformation must be either 'YeoJohnson' or 'log'.")
  }

  # ---- Outer resampling indices ----
  cv_outer_train_folds_rows <-
    caret::createMultiFolds(df[[outcome_var]], k = cv_outer_folds, times = cv_outer_repeats)

  # ---- Recipes ----
  make_recipe_vsurf <- function(data, outcome_var,transformation_rule) {
    all_pred <- setdiff(names(data), outcome_var)
    rec <- recipes::recipe(stats::reformulate(all_pred, response = outcome_var), data = data) |>
      recipes::step_string2factor(recipes::all_nominal_predictors())

    # conditional transformation
    if (transformation_rule == "log") {
      rec <- rec |>
        recipes::step_log(recipes::all_numeric_predictors(),
                          offset = 0.0001)
    } else if (transformation_rule == "YeoJohnson") {
      rec <- rec |>
        recipes::step_YeoJohnson(recipes::all_numeric_predictors())
    }

    # Interactions are created after imputation and before scaling so products are scaled too.
    rec <- add_imputation_steps(
      rec,
      data,
      outcome_var,
      numeric_imputation,
      nominal_imputation,
      add_missing_indicators
    )
    rec <- add_interaction_dummy_steps(rec, add_interactions, interaction_categorical_vars)
    rec <- add_interaction_steps(rec, add_interactions, interaction_vars, interaction_categorical_vars)
    rec <- add_polynomial_steps(rec, add_polynomials, poly_vars, poly_degree)
    rec <- finish_recipe_preprocessing(
      rec,
      data,
      outcome_var,
      remove_zero_variance,
      remove_correlated_predictors,
      correlation_threshold
    )

    rec
  }

  make_recipe_glmnet_design <- function(data, outcome_var,transformation_rule) {
    all_pred <- setdiff(names(data), outcome_var)
    rec <- recipes::recipe(stats::reformulate(all_pred, response = outcome_var), data = data) |>
      recipes::step_string2factor(recipes::all_nominal_predictors())

    # conditional transformation
    if (transformation_rule == "log") {
      rec <- rec |>
        recipes::step_log(recipes::all_numeric_predictors(),
                          offset = 0.0001)
    } else if (transformation_rule == "YeoJohnson") {
      rec <- rec |>
        recipes::step_YeoJohnson(recipes::all_numeric_predictors())
    }

    # Dummy variables are created before the requested interactions
    # so the glmnet recipe stays fully numeric.
    rec <- add_imputation_steps(
      rec,
      data,
      outcome_var,
      numeric_imputation,
      nominal_imputation,
      add_missing_indicators
    ) |>
      recipes::step_dummy(recipes::all_nominal_predictors(), one_hot = FALSE)
    rec <- add_interaction_steps(rec, add_interactions, interaction_vars, interaction_categorical_vars)
    rec <- add_polynomial_steps(rec, add_polynomials, poly_vars, poly_degree)
    rec <- finish_recipe_preprocessing(
      rec,
      data,
      outcome_var,
      remove_zero_variance,
      remove_correlated_predictors,
      correlation_threshold
    )
    rec
  }

  elastic_alpha_grid <- alpha_grid[alpha_grid < 1]
  if (!length(elastic_alpha_grid)) {
    stop("alpha_grid must include at least one value below 1; lasso uses alpha = 1 as a separate model.")
  }
  tg_elastic <- make_glmnet_tune_grid(
    elastic_alpha_grid,
    lambda_grid,
    hyperparam_search,
    random_search_size,
    "elastic net"
  )
  tg_lasso <- make_glmnet_tune_grid(
    1,
    lambda_grid,
    hyperparam_search,
    random_search_size,
    "lasso"
  )

  # ---- Inner CV control (regression) ----
  if (selection_rule == "best") {
    ctrl_inner <- caret::trainControl(
      method = cv_method,
      number = if (cv_method == "repeatedcv") inner_cv_folds else NA,
      repeats = if (cv_method == "repeatedcv") inner_cv_repeats else NA,
      savePredictions = "final",
      summaryFunction = caret::defaultSummary,
      verboseIter = FALSE,
      allowParallel = TRUE
    )
  } else if (selection_rule == "oneSE") {
    ctrl_inner <- caret::trainControl(
      method = cv_method,
      number = if (cv_method == "repeatedcv") inner_cv_folds else NA,
      repeats = if (cv_method == "repeatedcv") inner_cv_repeats else NA,
      savePredictions = "final",
      summaryFunction = caret::defaultSummary,
      selectionFunction = "oneSE",
      verboseIter = FALSE,
      allowParallel = TRUE
    )
  }

  n_outer <- length(cv_outer_train_folds_rows)

  sel_vars_df <- tibble::tibble(
    elastic_sel_vars = vector("list", n_outer),
    lasso_sel_vars   = vector("list", n_outer),
    vsurf_sel_vars_interp   = vector("list", n_outer),
    vsurf_sel_vars_pred   = vector("list", n_outer),
    elastic_coefs      = vector("list", n_outer),
    lasso_coefs        = vector("list", n_outer)
  )

  outer_fold_predictions <- vector("list", n_outer)
  # The mean baseline is intentionally constant across every assessment row.
  baseline_mean <- mean(df[[outcome_var]], na.rm = TRUE)

  # ---- Parallel backend (once) ----
  cl <- parallel::makePSOCKcluster(max(1L, parallel::detectCores() - 1L))
  doParallel::registerDoParallel(cl)
  on.exit({
    try(parallel::stopCluster(cl), silent = TRUE)
    foreach::registerDoSEQ()
  }, add = TRUE)

  # ---- Outer CV ----
  for (i in seq_along(cv_outer_train_folds_rows)) {
    outer_train_idx <- cv_outer_train_folds_rows[[i]]
    outer_test_idx  <- setdiff(seq_len(nrow(df)), outer_train_idx)
    outer_d_train   <- df[outer_train_idx, , drop = FALSE]
    outer_d_test    <- df[outer_test_idx,  , drop = FALSE]

    # VSURF bake
    rec_vsurf_full <- make_recipe_vsurf(outer_d_train, outcome_var,transformation_rule)
    prep_rec_vsurf_full <- recipes::prep(rec_vsurf_full, training = outer_d_train, retain = TRUE)
    d_train_baked <- recipes::bake(prep_rec_vsurf_full, new_data = NULL)
    d_test_baked  <- recipes::bake(prep_rec_vsurf_full, new_data = outer_d_test)

    x_train <- d_train_baked |> dplyr::select(-dplyr::all_of(outcome_var)) |> as.data.frame()
    y_train <- d_train_baked[[outcome_var]]
    x_test  <- d_test_baked  |> dplyr::select(-dplyr::all_of(outcome_var)) |> as.data.frame()


    # VSURF
    vs <- VSURF::VSURF(
      y = y_train,
      x = x_train,
      ntree.thres   = ntree, nfor.thres  = nforests,
      ntree.interp  = ntree, nfor.interp = nforests,
      ntree.pred    = ntree, nfor.pred   = nforests,
      RFimplem = "ranger",
      parallel = FALSE
    )
    sel_idx_inter <- vs$varselect.interp
    sel_idx_pred <- vs$varselect.pred
    vsurf_sel_vars_in <- if (length(sel_idx_inter)) colnames(x_train)[sel_idx_inter] else character(0)
    vsurf_sel_vars_pr <- if (length(sel_idx_pred)) colnames(x_train)[sel_idx_pred] else character(0)
    sel_vars_df$vsurf_sel_vars_interp[[i]] <- vsurf_sel_vars_in
    sel_vars_df$vsurf_sel_vars_pred[[i]] <- vsurf_sel_vars_pr
    message(sprintf("OUTER %02d - VSURF selected at the interpretation step (%d): %s",
                    i, length(vsurf_sel_vars_in), paste(vsurf_sel_vars_in, collapse=", ")))
    message(sprintf("OUTER %02d - VSURF selected at the prediction step (%d): %s",
                    i, length(vsurf_sel_vars_pr), paste(vsurf_sel_vars_pr, collapse=", ")))

    # Keep the original training data inside the VSURF object so predict()
    # remains self-contained when called later in the loop.
    vs$x_train_saved <- x_train
    vs$y_train_saved <- y_train
    vs$call$x <- quote(object$x_train_saved)
    vs$call$y <- quote(object$y_train_saved)

    # Predictions
    p_vsurf_inter <- predict(vs, newdata = x_test, step = "interp")

    if (is.null(vs$varselect.pred) || length(vs$varselect.pred) == 0) {
      # fallback if pred-step selects nothing
      # option A: use interp preds
      p_vsurf_pred <- p_vsurf_inter

      # option B instead (baseline mean on outer train):
      # p_vsurf_pred <- rep(mean(y_train, na.rm = TRUE), nrow(x_test))
    } else {
      p_vsurf_pred <- predict(vs, newdata = x_test, step = "pred")
    }

    # Elastic-net fit and coefficient selection on the glmnet design.
    rec_glmnet <- make_recipe_glmnet_design(outer_d_train, outcome_var,transformation_rule)
    elastic_sel <- caret::train(
      rec_glmnet,
      data = outer_d_train,
      method = "glmnet",
      metric = metric,
      maximize = FALSE,
      tuneGrid = tg_elastic,
      trControl = ctrl_inner,
      family = "gaussian",
      standardize = FALSE
    )
    coefs <- as.matrix(coef(elastic_sel$finalModel,
                            s = elastic_sel$bestTune$lambda[1]))[, 1]
    el_sel_vars <- setdiff(names(coefs)[coefs != 0], "(Intercept)")
    sel_vars_df$elastic_sel_vars[[i]] <- el_sel_vars
    sel_vars_df$elastic_coefs[[i]] <- coefs
    message(sprintf("OUTER %02d - Elastic selected (%d): %s",
                    i, length(el_sel_vars), paste(el_sel_vars, collapse=", ")))

    # One-pass elastic predictions on outer test
    p_elas <- predict(elastic_sel, newdata = outer_d_test)

    # Lasso selection on ALL predictors (dummy'd recipe)
    lasso_sel <- caret::train(
      rec_glmnet,
      data = outer_d_train,
      method = "glmnet",
      metric = metric,
      maximize = FALSE,
      tuneGrid = tg_lasso,
      trControl = ctrl_inner,
      family = "gaussian",
      standardize = FALSE
    )
    lasso_coefs <- as.matrix(coef(lasso_sel$finalModel,
                                  s = lasso_sel$bestTune$lambda[1]))[, 1]
    lasso_sel_vars <- setdiff(names(lasso_coefs)[lasso_coefs != 0], "(Intercept)")
    sel_vars_df$lasso_sel_vars[[i]] <- lasso_sel_vars
    sel_vars_df$lasso_coefs[[i]] <- lasso_coefs
    message(sprintf("OUTER %02d - Lasso selected (%d): %s",
                    i, length(lasso_sel_vars), paste(lasso_sel_vars, collapse=", ")))

    p_lasso <- predict(lasso_sel, newdata = outer_d_test)

    baseline_pred <- rep(baseline_mean, nrow(outer_d_test))

    # Store per-row predictions (align to indices)
    outer_fold_predictions[[i]] <- list(
      row_id = outer_test_idx,
      outcome = outer_d_test[[outcome_var]],
      pred_elastic_net = p_elas,
      pred_lasso = p_lasso,
      pred_baseline_mean = baseline_pred,
      pred_vsurf_interpretation = p_vsurf_inter,
      pred_vsurf_prediction = p_vsurf_pred
    )

  } # end outer

  # ---- Aggregate predictions over outer folds ----
  outer_predictions_long <- tibble::tibble()

  for (i in seq_along(outer_fold_predictions)) {
    outer_predictions_long <- dplyr::bind_rows(
      outer_predictions_long,
      tibble::tibble(
        row_id = outer_fold_predictions[[i]]$row_id,
        pred_elastic_net = outer_fold_predictions[[i]]$pred_elastic_net,
        pred_lasso = outer_fold_predictions[[i]]$pred_lasso,
        outcome = outer_fold_predictions[[i]]$outcome,
        pred_baseline_mean = outer_fold_predictions[[i]]$pred_baseline_mean,
        pred_vsurf_interpretation = outer_fold_predictions[[i]]$pred_vsurf_interpretation,
        pred_vsurf_prediction = outer_fold_predictions[[i]]$pred_vsurf_prediction
      )
    )
  }

  # Average repeated outer-CV predictions to one row per original observation.
  outer_predictions <- outer_predictions_long |>
    dplyr::group_by(row_id) |>
    dplyr::summarise(
      outcome = dplyr::first(outcome),
      pred_elastic_net = mean(pred_elastic_net, na.rm = TRUE),
      pred_lasso = mean(pred_lasso, na.rm = TRUE),
      pred_baseline_mean = mean(pred_baseline_mean, na.rm = TRUE),
      pred_vsurf_interpretation = mean(pred_vsurf_interpretation, na.rm = TRUE),
      pred_vsurf_prediction = mean(pred_vsurf_prediction, na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::arrange(row_id)

  # Return one prediction table instead of parallel vectors with repeated ids/outcomes.
  return(list(
    outer_predictions = outer_predictions,
    baseline_mean = baseline_mean,
    selected_variables = sel_vars_df,
    tuning = list(
      search = hyperparam_search,
      random_search_size = if (hyperparam_search == "random") random_search_size else NA_integer_,
      elastic_net_candidates = tg_elastic,
      lasso_candidates = tg_lasso
    ),
    preprocessing = list(
      numeric_imputation = numeric_imputation,
      nominal_imputation = nominal_imputation,
      add_missing_indicators = add_missing_indicators,
      remove_zero_variance = remove_zero_variance,
      remove_correlated_predictors = remove_correlated_predictors,
      correlation_threshold = correlation_threshold,
      add_interactions = add_interactions,
      interaction_vars = interaction_vars,
      add_polynomials = add_polynomials,
      poly_vars = poly_vars,
      poly_degree = poly_degree
    )
  ))
}


outer_perf_nested_continuous <- function(trained_object) {
  preds <- trained_object$outer_predictions
  selected_variables <- trained_object$selected_variables

  stab_table <- function(sel_list) {
    selected <- unlist(sel_list)
    if (!length(selected)) {
      return(tibble::tibble(variable = character(), times_selected = integer(), freq = numeric()))
    }
    tbl <- sort(table(selected), decreasing = TRUE)
    tibble::tibble(variable = names(tbl), times_selected = as.integer(tbl)) |>
      dplyr::mutate(freq = times_selected / nrow(selected_variables))
  }

  cat("\n--- Elastic Net: selection stability ---\n")
  print(stab_table(selected_variables$elastic_sel_vars))
  cat("\n--- Lasso: selection stability ---\n")
  print(stab_table(selected_variables$lasso_sel_vars))
  cat("\n--- VSURF: selection stability at the interpretation step ---\n")
  print(stab_table(selected_variables$vsurf_sel_vars_interp))
  cat("\n--- VSURF: selection stability at the prediction step ---\n")
  print(stab_table(selected_variables$vsurf_sel_vars_pred))

  get_metrics <- function(y_true, y_pred) {
    dat <- tibble::tibble(truth = y_true, pred = y_pred)
    rmse = yardstick::rmse(dat, truth, pred, na_rm = TRUE) |> dplyr::pull(.estimate)
    mae  = yardstick::mae(dat,  truth, pred, na_rm = TRUE) |> dplyr::pull(.estimate)
    r2   = yardstick::rsq(dat,  truth, pred, na_rm = TRUE) |> dplyr::pull(.estimate)
    r2_trad = yardstick::rsq_trad(dat, truth, pred, na_rm = TRUE) |> dplyr::pull(.estimate)
    ccc  = yardstick::ccc(dat,  truth, pred, na_rm = TRUE) |> dplyr::pull(.estimate)
    mape = yardstick::mape(dat,  truth, pred, na_rm = TRUE) |> dplyr::pull(.estimate)
    smape = yardstick::smape(dat, truth, pred, na_rm = TRUE) |> dplyr::pull(.estimate)

    fit_citl  <- stats::lm(y_true ~ 1 + offset(y_pred))
    citl      <- unname(coef(fit_citl)[1])

    if (stats::sd(y_pred, na.rm = TRUE) == 0) {
      cal_int <- NA_real_
      cal_slope <- NA_real_
    } else {
      fit_cal <- stats::lm(y_true ~ y_pred)
      cal_int <- unname(coef(fit_cal)[1])
      cal_slope <- unname(coef(fit_cal)[2])
    }

    tibble::tibble(
      rmse = rmse, mae = mae, r2 = r2, r2_trad = r2_trad, ccc = ccc, mape = mape, smape = smape,
      cal_intercept = cal_int, cal_slope = cal_slope, citl_intercept = citl
    )
  }

  single_elas   <- get_metrics(preds$outcome, preds$pred_elastic_net)
  lasso_mod     <- get_metrics(preds$outcome, preds$pred_lasso)
  baseline_mod  <- get_metrics(preds$outcome, preds$pred_baseline_mean)
  vsurf_mod_inter  <- get_metrics(preds$outcome, preds$pred_vsurf_interpretation)
  vsurf_mod_pred  <- get_metrics(preds$outcome, preds$pred_vsurf_prediction)

  dplyr::bind_rows(
    dplyr::as_tibble(single_elas) |> dplyr::mutate(method = "Elastic net (outer)"),
    dplyr::as_tibble(lasso_mod) |> dplyr::mutate(method = "Lasso (outer)"),
    dplyr::as_tibble(vsurf_mod_inter) |> dplyr::mutate(method = "VSURF interpretation step (outer)"),
    dplyr::as_tibble(vsurf_mod_pred) |> dplyr::mutate(method = "VSURF prediction step (outer)"),
    dplyr::as_tibble(baseline_mod) |> dplyr::mutate(method = "Reference mean model (constant)")
  ) |> dplyr::select(method, dplyr::everything())
}

#################################################################FINAL MODEL WITH COEFS#################################################################

final_model_with_coefs <- function(df,
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
                                   ) {

  family <- match.arg(family)
  cv_method <- match.arg(cv_method)
  selection_rule <- match.arg(selection_rule)
  metric <- match.arg(cont_optim_metric)
  transformation_rule <- match.arg(data_transformation)
  hyperparam_search <- match.arg(hyperparam_search)
  numeric_imputation <- match.arg(numeric_imputation)
  nominal_imputation <- match.arg(nominal_imputation)
  add_missing_indicators <- match.arg(add_missing_indicators)
  remove_zero_variance <- match.arg(remove_zero_variance)
  remove_correlated_predictors <- match.arg(remove_correlated_predictors)
  correlation_threshold <- validate_correlation_threshold(correlation_threshold)
  add_interactions <- match.arg(add_interactions)
  add_polynomials <- match.arg(add_polynomials)
  post_selection_refit <- match.arg(post_selection_refit)
  interaction_vars <- validate_interaction_vars(df, outcome_var, add_interactions, interaction_vars)
  interaction_categorical_vars <- categorical_interaction_vars(df, outcome_var, interaction_vars)
  poly_vars <- validate_poly_vars(df, outcome_var, add_polynomials, poly_vars, poly_degree)
  poly_degree <- as.integer(poly_degree)
  set.seed(1)

  if (selection_rule == "best") {
    message("Final model will use the hyperparameters with the best inner CV performance.")
  } else if (selection_rule == "oneSE") {
    message("Final model will use the most regularized hyperparameters (parsimonious model) within 1 SE of the best inner CV performance.")
  } else  {
    stop("selection_rule must be either 'best' or 'oneSE'.")
  }

  if (selection_rule == "oneSE" && cv_method == "LOOCV") {
    stop("selection_rule = 'oneSE' is not compatible with LOOCV. Choose cv_method = 'repeatedcv' instead.")
  }

  message(sprintf(
    "Hyperparameter search: %s%s",
    hyperparam_search,
    if (hyperparam_search == "random") sprintf(" (%d random candidates per glmnet model)", random_search_size) else ""
  ))
  message(sprintf("Post-selection refit: %s", post_selection_refit))
  message(sprintf(
    "Imputation: numeric=%s, nominal=%s, missing indicators=%s, remove zero variance=%s, remove correlated predictors=%s, interactions=%s, polynomials=%s",
    numeric_imputation,
    nominal_imputation,
    add_missing_indicators,
    remove_zero_variance,
    format_correlation_setting(remove_correlated_predictors, correlation_threshold),
    format_interaction_setting(add_interactions, interaction_vars),
    format_poly_setting(add_polynomials, poly_vars, poly_degree)
  ))

  #Transformation rule
  if (transformation_rule == "YeoJohnson") {
    message("Data transformation: Yeo-Johnson")
  } else if (transformation_rule == "log") {
    message("Data transformation: Natural log-transform (0.0001 offset)")
  } else {
    stop("data_transformation must be either 'YeoJohnson' or 'log'.")
  }

  elastic_alpha_grid <- alpha_grid[alpha_grid < 1]
  if (!length(elastic_alpha_grid)) {
    stop("alpha_grid must include at least one value below 1; lasso uses alpha = 1 as a separate model.")
  }
  tg_elastic <- make_glmnet_tune_grid(
    elastic_alpha_grid,
    lambda_grid,
    hyperparam_search,
    random_search_size,
    "elastic net"
  )
  tg_lasso <- make_glmnet_tune_grid(
    1,
    lambda_grid,
    hyperparam_search,
    random_search_size,
    "lasso"
  )

  post_selection_targets <- switch(
    post_selection_refit,
    vsurf_only = "VSURF",
    all = c("VSURF", "ElasticNet", "Lasso"),
    none = character(0)
  )

  if (family == "binomial") {
    message(sprintf("IMPORTANT: All variables and coefficients are selected based on the predicted probability of the positive class: '%s'", positive_class))

    df[[outcome_var]] <- relevel(factor(df[[outcome_var]]), ref = negative_class)

    make_recipe_vsurf <- function(selected_vars, data, outcome_var,transformation_rule) {
      all_pred <- setdiff(names(data), outcome_var)
      rec <- recipes::recipe(stats::reformulate(all_pred, response = outcome_var), data = data) |>
        recipes::step_string2factor(recipes::all_nominal_predictors())

      # conditional transformation
      if (transformation_rule == "log") {
        rec <- rec |>
          recipes::step_log(recipes::all_numeric_predictors(),
                            offset = 0.0001)
      } else if (transformation_rule == "YeoJohnson") {
        rec <- rec |>
          recipes::step_YeoJohnson(recipes::all_numeric_predictors())
      }

      # Keep selected predictors plus any ingredients needed to recreate
      # selected interaction columns.
      rec <- add_imputation_steps(
        rec,
        data,
        outcome_var,
        numeric_imputation,
        nominal_imputation,
        add_missing_indicators
      )

      if (!is.null(selected_vars) && add_interactions == "no" && add_polynomials == "no") {
        rec <- rec |>
          recipes::step_rm(
            recipes::all_predictors(),
            -tidyselect::any_of(selected_vars_before_interactions(selected_vars, add_interactions, interaction_vars))
          )
      }
      rec <- add_interaction_dummy_steps(rec, add_interactions, interaction_categorical_vars)
      rec <- add_interaction_steps(rec, add_interactions, interaction_vars, interaction_categorical_vars)
      rec <- add_polynomial_steps(rec, add_polynomials, poly_vars, poly_degree)
      if (!is.null(selected_vars)) {
        rec <- rec |>
          recipes::step_rm(recipes::all_predictors(), -tidyselect::any_of(selected_vars))
      }
      rec <- finish_recipe_preprocessing(
        rec,
        data,
        outcome_var,
        remove_zero_variance,
        remove_correlated_predictors,
        correlation_threshold
      )
      rec
    }

    # VSURF returns raw predictor names; glmnet post-selection refits need
    # those predictors dummy-expanded after selection.
    make_recipe_vsurf_glmnet <- function(selected_vars, data, outcome_var, transformation_rule) {
      stopifnot(!is.null(selected_vars))  # only for selected-var refits
      rec <- make_recipe_vsurf(selected_vars, data, outcome_var, transformation_rule) |>
        recipes::step_dummy(recipes::all_nominal_predictors(), one_hot = FALSE)
      rec <- finish_recipe_preprocessing(
        rec,
        data,
        outcome_var,
        remove_zero_variance,
        remove_correlated_predictors,
        correlation_threshold
      )
      rec
    }


    # Glmnet fits use the full numeric design matrix: factors are dummy-expanded
    # before elastic-net or lasso coefficients are estimated.
    make_recipe_glmnet_design <- function(selected_vars, data, outcome_var,transformation_rule) {
      all_pred <- setdiff(names(data), outcome_var)
      rec <- recipes::recipe(stats::reformulate(all_pred, response = outcome_var), data = data) |>
        recipes::step_string2factor(recipes::all_nominal_predictors())

      # conditional transformation
      if (transformation_rule == "log") {
        rec <- rec |>
          recipes::step_log(recipes::all_numeric_predictors(),
                            offset = 0.0001)
      } else if (transformation_rule == "YeoJohnson") {
        rec <- rec |>
          recipes::step_YeoJohnson(recipes::all_numeric_predictors())
      }

      # Dummy categorical predictors before interaction creation because glmnet
      # works with the baked numeric design matrix.
      rec <- add_imputation_steps(
        rec,
        data,
        outcome_var,
        numeric_imputation,
        nominal_imputation,
        add_missing_indicators
      ) |>
        recipes::step_dummy(recipes::all_nominal_predictors(), one_hot = FALSE)

      if (!is.null(selected_vars) && add_interactions == "no" && add_polynomials == "no") {
        rec <- rec |>
          recipes::step_rm(
            recipes::all_predictors(),
            -tidyselect::any_of(selected_vars_before_interactions(selected_vars, add_interactions, interaction_vars))
          )
      }
      rec <- add_interaction_steps(rec, add_interactions, interaction_vars, interaction_categorical_vars)
      rec <- add_polynomial_steps(rec, add_polynomials, poly_vars, poly_degree)
      if (!is.null(selected_vars)) {
        rec <- rec |>
          recipes::step_rm(recipes::all_predictors(), -tidyselect::any_of(selected_vars))
      }
      rec <- finish_recipe_preprocessing(
        rec,
        data,
        outcome_var,
        remove_zero_variance,
        remove_correlated_predictors,
        correlation_threshold
      )
      rec
    }

    rec_vsurf_full <- make_recipe_vsurf(NULL, df, outcome_var,transformation_rule)
    d_train_baked <- recipes::bake(recipes::prep(rec_vsurf_full, training = df), new_data = NULL)

    x_train <- d_train_baked |>
      dplyr::select(-dplyr::all_of(outcome_var)) |>
      as.data.frame()

    y_train <- d_train_baked[[outcome_var]]

    # ---- VSURF ----
    vs <- VSURF::VSURF(
      y = y_train,
      x = x_train,
      ntree.thres   = ntree, nfor.thres  = nforests,
      ntree.interp  = ntree, nfor.interp = nforests,
      ntree.pred    = ntree, nfor.pred   = nforests,
      RFimplem = "ranger",
      parallel = FALSE
    )
    sel_idx_int <- vs$varselect.interp
    sel_idx_pred <- vs$varselect.pred
    vsurf_sel_vars_int <- if (length(sel_idx_int)) colnames(x_train)[sel_idx_int] else character(0)
    vsurf_sel_vars_pred <- if (length(sel_idx_pred)) colnames(x_train)[sel_idx_pred] else character(0)
    message(sprintf("VSURF at the interpretation step selected (%d): %s",
            length(vsurf_sel_vars_int), paste(vsurf_sel_vars_int, collapse=", ")))
    message(sprintf("VSURF at the prediction step selected (%d): %s",
            length(vsurf_sel_vars_pred), paste(vsurf_sel_vars_pred, collapse=", ")))

    # ---- Glmnet tuning control ----
    if (selection_rule == "best") {
      ctrl_inner <- caret::trainControl(
        method = cv_method,
        number = if (cv_method == "repeatedcv") cv_folds else NA,
        repeats = if (cv_method == "repeatedcv") cv_repeats else NA,
        savePredictions = "final",
        classProbs = TRUE,
        summaryFunction = caret::mnLogLoss,
        verboseIter = FALSE,
        allowParallel = TRUE
      )
    } else if (selection_rule == "oneSE") {
      ctrl_inner <- caret::trainControl(
        method = cv_method,
        number = if (cv_method == "repeatedcv") cv_folds else NA,
        repeats = if (cv_method == "repeatedcv") cv_repeats else NA,
        savePredictions = "final",
        classProbs = TRUE,
        summaryFunction = caret::mnLogLoss,
        selectionFunction = "oneSE",
        verboseIter = FALSE,
        allowParallel = TRUE
      )
    }

    # ---- Elastic-net and lasso selection on the full glmnet design ----
    rec_glmnet_full <- make_recipe_glmnet_design(NULL, df, outcome_var,transformation_rule)

    elas_pre <- caret::train(
      rec_glmnet_full,
      data = df,
      method = "glmnet",
      metric = "logLoss",
      maximize = FALSE,
      tuneGrid = tg_elastic,
      trControl = ctrl_inner,
      family = "binomial",
      standardize = FALSE
    )

    # Extract the original penalized coefficients from the full-data fit.
    coefs <- as.matrix(coef(elas_pre$finalModel,
                            s = elas_pre$bestTune$lambda[1]))[, 1]
    el_sel_vars <- setdiff(names(coefs)[coefs != 0], "(Intercept)")
    coef_elas_only <- coefs[coefs != 0]
    message(sprintf("ElasticNet selected (%d): %s",
                    length(el_sel_vars), paste(el_sel_vars, collapse=", ")))
    message("ElasticNet non-zero penalized coefficients from the full-data fit:")
    print(coef_elas_only)

    lasso_pre <- caret::train(
      rec_glmnet_full,
      data = df,
      method = "glmnet",
      metric = "logLoss",
      maximize = FALSE,
      tuneGrid = tg_lasso,
      trControl = ctrl_inner,
      family = "binomial",
      standardize = FALSE
    )

    lasso_coefs <- as.matrix(coef(lasso_pre$finalModel,
                                  s = lasso_pre$bestTune$lambda[1]))[, 1]
    lasso_sel_vars <- setdiff(names(lasso_coefs)[lasso_coefs != 0], "(Intercept)")
    coef_lasso_only <- lasso_coefs[lasso_coefs != 0]
    message(sprintf("Lasso selected (%d): %s",
                    length(lasso_sel_vars), paste(lasso_sel_vars, collapse=", ")))
    message("Lasso non-zero penalized coefficients from the full-data fit:")
    print(coef_lasso_only)

    selected_variables <- tibble::tibble(
      vsurf_interpretation = list(vsurf_sel_vars_int),
      vsurf_prediction = list(vsurf_sel_vars_pred),
      elastic_net = list(el_sel_vars),
      lasso = list(lasso_sel_vars)
    )

    post_selection_fits <- list()
    post_selection_coefs <- stats::setNames(
      vector("list", length(post_selection_targets)),
      post_selection_targets
    )

    if (!length(post_selection_targets)) {
      message("No post-selection refits requested; returning selected variables and direct penalized glmnet coefficients.")
    }

    for (v in post_selection_targets) {
      message(sprintf("\nFitting post-selection model using %s-selected variables...", v))
      if (v == "VSURF") {
        sel_vars <- vsurf_sel_vars_int
      } else if (v == "Lasso") {
        sel_vars <- lasso_sel_vars
      } else {
        sel_vars <- el_sel_vars
      }

      if (length(sel_vars) == 0) {
        rec_main <- recipes::recipe(stats::as.formula(paste(outcome_var, "~ 1")), data = df)
        final_fit <- caret::train(
          rec_main, data = df,
          method = "glm",
          metric = "logLoss", maximize = FALSE,
          trControl = caret::trainControl(method = "none"),  # no CV needed
          family = "binomial"
        )

        post_selection_fits[[v]] <- final_fit
        post_selection_coefs[[v]] <- coef(final_fit$finalModel)

      } else if (length(sel_vars) == 1) {
        # VSURF names raw predictors; elastic net and lasso name post-dummy predictors.
        rec_main <- if (v == "VSURF") {
          make_recipe_vsurf(sel_vars, df, outcome_var,transformation_rule)
        } else {
          make_recipe_glmnet_design(sel_vars, df, outcome_var,transformation_rule)
        }
        final_fit <- caret::train(
          rec_main, data = df,
          method = "glm",
          metric = "logLoss", maximize = FALSE,
          trControl = caret::trainControl(method = "none"),  # no CV needed
          family = "binomial"
        )

        post_selection_fits[[v]] <- final_fit
        post_selection_coefs[[v]] <- coef(final_fit$finalModel)

      } else {
        if (v == "VSURF") {
          rec_main <- make_recipe_vsurf_glmnet(sel_vars, df, outcome_var,transformation_rule)
        } else {
          rec_main <- make_recipe_glmnet_design(sel_vars, df, outcome_var,transformation_rule)
        }
        tune_grid <- if (v == "Lasso") tg_lasso else tg_elastic

        final_fit <- caret::train(
          rec_main, data = df,
          method = "glmnet",
          metric = "logLoss", maximize = FALSE,
          tuneGrid = tune_grid,
          trControl = ctrl_inner,
          family = "binomial",
          standardize = FALSE
        )

        post_selection_fits[[v]] <- final_fit
        coefs_2 <- as.matrix(coef(final_fit$finalModel,
                                  s = final_fit$bestTune$lambda[1]))[, 1]
        post_selection_coefs[[v]] <- coefs_2[coefs_2 != 0]
      }
      message(sprintf("Post-selection coefficients after %s selection:\n", v))
      print(post_selection_coefs[[v]])
    }

    coefficient_estimates <- tibble::tibble(
      elastic_net_penalized = list(coef_elas_only),
      lasso_penalized = list(coef_lasso_only),
      post_selection = list(post_selection_coefs)
    )

    return(list(
      post_selection_models = post_selection_fits,
      selected_variables = selected_variables,
      coefficient_estimates = coefficient_estimates,
      vsurf = vs,
      elastic_net = elas_pre,
      lasso = lasso_pre,
      tuning = list(
        search = hyperparam_search,
        random_search_size = if (hyperparam_search == "random") random_search_size else NA_integer_,
        elastic_net_candidates = tg_elastic,
        lasso_candidates = tg_lasso
      ),
      preprocessing = list(
        numeric_imputation = numeric_imputation,
        nominal_imputation = nominal_imputation,
        add_missing_indicators = add_missing_indicators,
        remove_zero_variance = remove_zero_variance,
        remove_correlated_predictors = remove_correlated_predictors,
        correlation_threshold = correlation_threshold,
        add_interactions = add_interactions,
        interaction_vars = interaction_vars,
        add_polynomials = add_polynomials,
        poly_vars = poly_vars,
        poly_degree = poly_degree,
        post_selection_refit = post_selection_refit
      )
    ))
  }
  else if (family == "gaussian") {

    make_recipe_vsurf <- function(selected_vars, data, outcome_var,transformation_rule) {
      all_pred <- setdiff(names(data), outcome_var)
      rec <- recipes::recipe(stats::reformulate(all_pred, response = outcome_var), data = data) |>
        recipes::step_string2factor(recipes::all_nominal_predictors())

      # conditional transformation
      if (transformation_rule == "log") {
        rec <- rec |>
          recipes::step_log(recipes::all_numeric_predictors(),
                            offset = 0.0001)
      } else if (transformation_rule == "YeoJohnson") {
        rec <- rec |>
          recipes::step_YeoJohnson(recipes::all_numeric_predictors())
      }

      # Keep selected predictors plus any ingredients needed to recreate
      # selected interaction columns.
      rec <- add_imputation_steps(
        rec,
        data,
        outcome_var,
        numeric_imputation,
        nominal_imputation,
        add_missing_indicators
      )

      if (!is.null(selected_vars) && add_interactions == "no" && add_polynomials == "no") {
        rec <- rec |>
          recipes::step_rm(
            recipes::all_predictors(),
            -tidyselect::any_of(selected_vars_before_interactions(selected_vars, add_interactions, interaction_vars))
          )
      }
      rec <- add_interaction_dummy_steps(rec, add_interactions, interaction_categorical_vars)
      rec <- add_interaction_steps(rec, add_interactions, interaction_vars, interaction_categorical_vars)
      rec <- add_polynomial_steps(rec, add_polynomials, poly_vars, poly_degree)
      if (!is.null(selected_vars)) {
        rec <- rec |>
          recipes::step_rm(recipes::all_predictors(), -tidyselect::any_of(selected_vars))
      }
      rec <- finish_recipe_preprocessing(
        rec,
        data,
        outcome_var,
        remove_zero_variance,
        remove_correlated_predictors,
        correlation_threshold
      )
      rec
    }

    # VSURF returns raw predictor names; glmnet post-selection refits need
    # those predictors dummy-expanded after selection.
    make_recipe_vsurf_glmnet <- function(selected_vars, data, outcome_var, transformation_rule) {
      stopifnot(!is.null(selected_vars))  # only for selected-var refits
      rec <- make_recipe_vsurf(selected_vars, data, outcome_var, transformation_rule) |>
        recipes::step_dummy(recipes::all_nominal_predictors(), one_hot = FALSE)
      rec <- finish_recipe_preprocessing(
        rec,
        data,
        outcome_var,
        remove_zero_variance,
        remove_correlated_predictors,
        correlation_threshold
      )
      rec
    }


    # Glmnet fits use the full numeric design matrix: factors are dummy-expanded
    # before elastic-net or lasso coefficients are estimated.
    make_recipe_glmnet_design <- function(selected_vars, data, outcome_var,transformation_rule) {
      all_pred <- setdiff(names(data), outcome_var)
      rec <- recipes::recipe(stats::reformulate(all_pred, response = outcome_var), data = data) |>
        recipes::step_string2factor(recipes::all_nominal_predictors())

      # conditional transformation
      if (transformation_rule == "log") {
        rec <- rec |>
          recipes::step_log(recipes::all_numeric_predictors(),
                            offset = 0.0001)
      } else if (transformation_rule == "YeoJohnson") {
        rec <- rec |>
          recipes::step_YeoJohnson(recipes::all_numeric_predictors())
      }

      # Dummy categorical predictors before interaction creation because glmnet
      # works with the baked numeric design matrix.
      rec <- add_imputation_steps(
        rec,
        data,
        outcome_var,
        numeric_imputation,
        nominal_imputation,
        add_missing_indicators
      ) |>
        recipes::step_dummy(recipes::all_nominal_predictors(), one_hot = FALSE)

      if (!is.null(selected_vars) && add_interactions == "no" && add_polynomials == "no") {
        rec <- rec |>
          recipes::step_rm(
            recipes::all_predictors(),
            -tidyselect::any_of(selected_vars_before_interactions(selected_vars, add_interactions, interaction_vars))
          )
      }
      rec <- add_interaction_steps(rec, add_interactions, interaction_vars, interaction_categorical_vars)
      rec <- add_polynomial_steps(rec, add_polynomials, poly_vars, poly_degree)
      if (!is.null(selected_vars)) {
        rec <- rec |>
          recipes::step_rm(recipes::all_predictors(), -tidyselect::any_of(selected_vars))
      }
      rec <- finish_recipe_preprocessing(
        rec,
        data,
        outcome_var,
        remove_zero_variance,
        remove_correlated_predictors,
        correlation_threshold
      )
      rec
    }

    rec_vsurf_full <- make_recipe_vsurf(NULL, df, outcome_var,transformation_rule)
    d_train_baked <- recipes::bake(recipes::prep(rec_vsurf_full, training = df), new_data = NULL)

    x_train <- d_train_baked |>
      dplyr::select(-dplyr::all_of(outcome_var)) |>
      as.data.frame()

    y_train <- d_train_baked[[outcome_var]]

    # ---- VSURF ----
    vs <- VSURF::VSURF(
      y = y_train,
      x = x_train,
      ntree.thres   = ntree, nfor.thres  = nforests,
      ntree.interp  = ntree, nfor.interp = nforests,
      ntree.pred    = ntree, nfor.pred   = nforests,
      RFimplem = "ranger",
      parallel = FALSE
    )
    sel_idx_int <- vs$varselect.interp
    sel_idx_pred <- vs$varselect.pred
    vsurf_sel_vars_int <- if (length(sel_idx_int)) colnames(x_train)[sel_idx_int] else character(0)
    vsurf_sel_vars_pred <- if (length(sel_idx_pred)) colnames(x_train)[sel_idx_pred] else character(0)
    message(sprintf("VSURF at the interpretation step selected (%d): %s",
                    length(vsurf_sel_vars_int), paste(vsurf_sel_vars_int, collapse=", ")))
    message(sprintf("VSURF at the prediction step selected (%d): %s",
                    length(vsurf_sel_vars_pred), paste(vsurf_sel_vars_pred, collapse=", ")))

    # ---- Glmnet tuning control ----
    if (selection_rule == "best") {
      ctrl_inner <- caret::trainControl(
        method = cv_method,
        number = if (cv_method == "repeatedcv") cv_folds else NA,
        repeats = if (cv_method == "repeatedcv") cv_repeats else NA,
        savePredictions = "final",
        summaryFunction = caret::defaultSummary,
        verboseIter = FALSE,
        allowParallel = TRUE
      )
    } else if (selection_rule == "oneSE") {
      ctrl_inner <- caret::trainControl(
        method = cv_method,
        number = if (cv_method == "repeatedcv") cv_folds else NA,
        repeats = if (cv_method == "repeatedcv") cv_repeats else NA,
        savePredictions = "final",
        summaryFunction = caret::defaultSummary,
        selectionFunction = "oneSE",
        verboseIter = FALSE,
        allowParallel = TRUE
      )
    }

    # ---- Elastic-net and lasso selection on the full glmnet design ----
    rec_glmnet_full <- make_recipe_glmnet_design(NULL, df, outcome_var,transformation_rule)

    elas_pre <- caret::train(
      rec_glmnet_full,
      data = df,
      method = "glmnet",
      metric = metric,
      maximize = FALSE,
      tuneGrid = tg_elastic,
      trControl = ctrl_inner,
      family = "gaussian",
      standardize = FALSE
    )

    # Extract the original penalized coefficients from the full-data fit.
    coefs <- as.matrix(coef(elas_pre$finalModel,
                            s = elas_pre$bestTune$lambda[1]))[, 1]
    el_sel_vars <- setdiff(names(coefs)[coefs != 0], "(Intercept)")
    coef_elas_only <- coefs[coefs != 0]
    message(sprintf("ElasticNet selected (%d): %s",
                    length(el_sel_vars), paste(el_sel_vars, collapse=", ")))
    message("ElasticNet non-zero penalized coefficients from the full-data fit:")
    print(coef_elas_only)

    lasso_pre <- caret::train(
      rec_glmnet_full,
      data = df,
      method = "glmnet",
      metric = metric,
      maximize = FALSE,
      tuneGrid = tg_lasso,
      trControl = ctrl_inner,
      family = "gaussian",
      standardize = FALSE
    )

    lasso_coefs <- as.matrix(coef(lasso_pre$finalModel,
                                  s = lasso_pre$bestTune$lambda[1]))[, 1]
    lasso_sel_vars <- setdiff(names(lasso_coefs)[lasso_coefs != 0], "(Intercept)")
    coef_lasso_only <- lasso_coefs[lasso_coefs != 0]
    message(sprintf("Lasso selected (%d): %s",
                    length(lasso_sel_vars), paste(lasso_sel_vars, collapse=", ")))
    message("Lasso non-zero penalized coefficients from the full-data fit:")
    print(coef_lasso_only)

    selected_variables <- tibble::tibble(
      vsurf_interpretation = list(vsurf_sel_vars_int),
      vsurf_prediction = list(vsurf_sel_vars_pred),
      elastic_net = list(el_sel_vars),
      lasso = list(lasso_sel_vars)
    )

    post_selection_fits <- list()
    post_selection_coefs <- stats::setNames(
      vector("list", length(post_selection_targets)),
      post_selection_targets
    )

    if (!length(post_selection_targets)) {
      message("No post-selection refits requested; returning selected variables and direct penalized glmnet coefficients.")
    }

    for (v in post_selection_targets) {
      message(sprintf("\nFitting post-selection model using %s-selected variables...", v))
      if (v == "VSURF") {
        sel_vars <- vsurf_sel_vars_int
      } else if (v == "Lasso") {
        sel_vars <- lasso_sel_vars
      } else {
        sel_vars <- el_sel_vars
      }

      if (length(sel_vars) == 0) {
        rec_main <- recipes::recipe(stats::as.formula(paste(outcome_var, "~ 1")), data = df)
        final_fit <- caret::train(
          rec_main, data = df,
          method = "glm",
          metric = metric, maximize = FALSE,
          trControl = caret::trainControl(method = "none"),  # no CV needed
          family = "gaussian"
        )

        post_selection_fits[[v]] <- final_fit
        post_selection_coefs[[v]] <- coef(final_fit$finalModel)

      } else if (length(sel_vars) == 1) {
        # VSURF names raw predictors; elastic net and lasso name post-dummy predictors.
        rec_main <- if (v == "VSURF") {
          make_recipe_vsurf(sel_vars, df, outcome_var,transformation_rule)
        } else {
          make_recipe_glmnet_design(sel_vars, df, outcome_var,transformation_rule)
        }
        final_fit <- caret::train(
          rec_main, data = df,
          method = "glm",
          metric = metric, maximize = FALSE,
          trControl = caret::trainControl(method = "none"),  # no CV needed
          family = "gaussian"
        )

        post_selection_fits[[v]] <- final_fit
        post_selection_coefs[[v]] <- coef(final_fit$finalModel)

      } else {
        if (v == "VSURF") {
          rec_main <- make_recipe_vsurf_glmnet(sel_vars, df, outcome_var,transformation_rule)
        } else {
          rec_main <- make_recipe_glmnet_design(sel_vars, df, outcome_var,transformation_rule)
        }
        tune_grid <- if (v == "Lasso") tg_lasso else tg_elastic

        final_fit <- caret::train(
          rec_main, data = df,
          method = "glmnet",
          metric = metric, maximize = FALSE,
          tuneGrid = tune_grid,
          trControl = ctrl_inner,
          family = "gaussian",
          standardize = FALSE
        )

        post_selection_fits[[v]] <- final_fit
        coefs_2 <- as.matrix(coef(final_fit$finalModel,
                                  s = final_fit$bestTune$lambda[1]))[, 1]
        post_selection_coefs[[v]] <- coefs_2[coefs_2 != 0]
      }
      message(sprintf("Post-selection coefficients after %s selection:\n", v))
      print(post_selection_coefs[[v]])
    }

    coefficient_estimates <- tibble::tibble(
      elastic_net_penalized = list(coef_elas_only),
      lasso_penalized = list(coef_lasso_only),
      post_selection = list(post_selection_coefs)
    )

    return(list(
      post_selection_models = post_selection_fits,
      selected_variables = selected_variables,
      coefficient_estimates = coefficient_estimates,
      vsurf = vs,
      elastic_net = elas_pre,
      lasso = lasso_pre,
      tuning = list(
        search = hyperparam_search,
        random_search_size = if (hyperparam_search == "random") random_search_size else NA_integer_,
        elastic_net_candidates = tg_elastic,
        lasso_candidates = tg_lasso
      ),
      preprocessing = list(
        numeric_imputation = numeric_imputation,
        nominal_imputation = nominal_imputation,
        add_missing_indicators = add_missing_indicators,
        remove_zero_variance = remove_zero_variance,
        remove_correlated_predictors = remove_correlated_predictors,
        correlation_threshold = correlation_threshold,
        add_interactions = add_interactions,
        interaction_vars = interaction_vars,
        add_polynomials = add_polynomials,
        poly_vars = poly_vars,
        poly_degree = poly_degree,
        post_selection_refit = post_selection_refit
      )
    ))
  } else {
    stop("family must be either 'binomial' or 'gaussian'.")
  }
}
