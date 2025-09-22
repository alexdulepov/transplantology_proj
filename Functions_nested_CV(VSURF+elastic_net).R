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
    ntree = 1000,
    nforests = 20,
    inner_cv_method = c("LOOCV", "repeatedcv"),
    inner_cv_repeats = 20,
    inner_cv_folds = 5,
    selection_rule = c("best", "oneSE")
) {
  # ---- Reproducible RNG ----
  RNGkind("L'Ecuyer-CMRG")
  set.seed(seed)
  cv_method <- match.arg(inner_cv_method)
  selection_rule <- match.arg(selection_rule)
  
  message(sprintf("IMPORTANT: All performance metrics and models are built around the predicted probability of the positive class: '%s'", positive_class))
  
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
  
  # ---- Outcome factor with explicit order: negative first ----
  if (!all(c(negative_class, positive_class) %in% unique(df[[outcome_var]]))) {
    stop("Outcome does not contain both specified classes.")
  }
  
  # Only examine predictors; avoid evaluating < 0 on non-numerics
  pred_names <- setdiff(names(df), outcome_var)
  num_pred   <- vapply(df[pred_names], is.numeric, TRUE)
  if (any(vapply(df[pred_names][num_pred], function(x) any(x < 0, na.rm = TRUE), TRUE))) {
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
  # VSURF: keep factors raw; NO dummying
  make_recipe_vsurf <- function(vars_to_keep, data, outcome_var) {
    all_pred <- setdiff(names(data), outcome_var)
    rec <- recipes::recipe(stats::reformulate(all_pred, response = outcome_var), data = data) |>
      recipes::step_string2factor(recipes::all_nominal_predictors()) |>
      recipes::step_impute_mode(recipes::all_nominal_predictors()) |>
      recipes::step_impute_knn(recipes::all_numeric_predictors(), neighbors = 5) |>
      recipes::step_log(recipes::all_numeric_predictors(), offset = 1) |>
      recipes::step_center(recipes::all_numeric_predictors()) |>
      recipes::step_scale(recipes::all_numeric_predictors())
    
    if (!is.null(vars_to_keep)) {
      rec <- rec |>
        recipes::step_rm(recipes::all_predictors(), -tidyselect::any_of(vars_to_keep))
    }
    rec
  }
  
  # Refit after VSURF for GLMNET: select RAW vars first, then dummy (k-1 coding)
  make_recipe_vsurf_glmnet <- function(vars_to_keep, data, outcome_var) {
    stopifnot(!is.null(vars_to_keep))  # only for selected-var refits
    rec <- make_recipe_vsurf(vars_to_keep, data, outcome_var) |>
      recipes::step_dummy(recipes::all_nominal_predictors(), one_hot = FALSE)
    rec
  }
  
  # ElasticNet self-selection: dummy FIRST → then (optionally) select dummy names
  make_recipe_after_elas <- function(vars_to_keep, data, outcome_var) {
    all_pred <- setdiff(names(data), outcome_var)
    rec <- recipes::recipe(stats::reformulate(all_pred, response = outcome_var), data = data) |>
      recipes::step_string2factor(recipes::all_nominal_predictors()) |>
      recipes::step_impute_mode(recipes::all_nominal_predictors()) |>
      recipes::step_impute_knn(recipes::all_numeric_predictors(), neighbors = 5) |>
      recipes::step_log(recipes::all_numeric_predictors(), offset = 1) |>
      recipes::step_center(recipes::all_numeric_predictors()) |>
      recipes::step_scale(recipes::all_numeric_predictors()) |>
      recipes::step_dummy(recipes::all_nominal_predictors(), one_hot = FALSE)
    
    if (!is.null(vars_to_keep)) {
      rec <- rec |>
        recipes::step_rm(recipes::all_predictors(), -tidyselect::any_of(vars_to_keep))
    }
    rec
  }
  
  tg <- expand.grid(alpha = alpha_grid, lambda = lambda_grid)
  
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
    vsurf_sel_vars   = vector("list", n_outer),
    final_coefs      = vector("list", n_outer)
  )
  
  n_alg <- 2L
  inner_perf <- vector("list", n_outer)
  for (i in seq_len(n_outer)) {
    inner_perf[[i]] <- vector("list", n_alg)
    names(inner_perf[[i]]) <- c("VSURF", "ElasticNet")
  }
  
  elas_once <- vector("list", n_outer)
  
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
    rec_all <- make_recipe_vsurf(NULL, outer_d_train, outcome_var)
    d_train_baked <- recipes::bake(recipes::prep(rec_all, training = outer_d_train), new_data = NULL)
    
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
    
    sel_idx <- vs$varselect.interp
    vsurf_sel_vars <- if (length(sel_idx)) colnames(x_train)[sel_idx] else character(0)
    sel_vars_df$vsurf_sel_vars[[i]] <- vsurf_sel_vars
    message(sprintf("OUTER %02d — VSURF selected (%d): %s",
                    i, length(vsurf_sel_vars), paste(vsurf_sel_vars, collapse=", ")))
    
    # ---- Elastic-net selection on all predictors ----
    rec_all <- make_recipe_after_elas(NULL, outer_d_train, outcome_var)
    
    elastic_sel <- caret::train(
      rec_all,
      data = outer_d_train,
      method = "glmnet",
      metric = "logLoss",
      maximize = FALSE,
      tuneGrid = tg,
      trControl = ctrl_inner,
      family = "binomial",
      standardize = FALSE
    )
    coefs <- as.matrix(coef(elastic_sel$finalModel,
                            s = elastic_sel$bestTune$lambda[1]))[, 1]
    el_sel_vars <- setdiff(names(coefs)[coefs != 0], "(Intercept)")
    sel_vars_df$elastic_sel_vars[[i]] <- el_sel_vars
    message(sprintf("OUTER %02d — Elastic selected (%d): %s",
                    i, length(el_sel_vars), paste(el_sel_vars, collapse=", ")))
    
    # predictions from single elastic net model (no variable selection)
    eps <- 1e-12
    p_elas <- predict(elastic_sel, newdata = outer_d_test, type = "prob")[[positive_class]]
    p_elas <- pmin(pmax(p_elas, eps), 1 - eps)
    
    # Intercept-only baseline fitted on OUTER TRAIN
    rec_prev <- recipes::recipe(stats::as.formula(paste(outcome_var, "~ 1")), data = outer_d_train)
    
    fit_prev <- caret::train(
      rec_prev, data = outer_d_train,
      method = "glm",
      metric = "logLoss", maximize = FALSE,
      trControl = caret::trainControl(method = "none"), family = "binomial"
    )
    
    p_prev <- predict(fit_prev, newdata = outer_d_test, type = "prob")[[positive_class]]
    baseline_pred <- pmin(pmax(p_prev, eps), 1 - eps)
    
    # Store per-row predictions (align to indices)
    elas_once[[i]] <- list(
      preds_elas_once    = p_elas,
      pred_idx_elas_once = outer_test_idx,
      y_elas_once        = outer_d_test[[outcome_var]],
      baseline_pred_once = baseline_pred
    )
    
    list_vars_both <- list(vsurf_sel_vars, el_sel_vars)
    
    # Store coefficients of final models for both methods
    sel_vars_df$final_coefs[[i]] <- vector("list", 2L)
    names(sel_vars_df$final_coefs[[i]]) <- c("VSURF", "ElasticNet")
    
    #---------------------- Inner CV with selected vars -------------------
    for (j in seq_along(list_vars_both)) {
      vars <- list_vars_both[[j]]
      
      if (j == 1) {
        # ---------------- VSURF branch ----------------
        if (length(vars) == 0) {
          rec_main <- recipes::recipe(stats::as.formula(paste(outcome_var, "~ 1")), data = outer_d_train)
          fit_inner <- caret::train(
            rec_main, data = outer_d_train,
            method = "glm",
            metric = "logLoss", maximize = FALSE,
            trControl = ctrl_inner, family = "binomial"
          )
          
        } else if (length(vars) == 1) {
          # raw factors, no dummying
          rec_main <- make_recipe_vsurf(vars, outer_d_train, outcome_var)
          fit_inner <- caret::train(
            rec_main, data = outer_d_train,
            method = "glm",
            metric = "logLoss", maximize = FALSE,
            trControl = ctrl_inner, family = "binomial"
          )
          
        } else {
          # penalized refit after VSURF: select raw vars → dummy for glmnet
          rec_main <- make_recipe_vsurf_glmnet(vars, outer_d_train, outcome_var)
          fit_inner <- caret::train(
            rec_main, data = outer_d_train,
            method = "glmnet",
            metric = "logLoss", maximize = FALSE,
            tuneGrid = tg, trControl = ctrl_inner,
            family = "binomial", standardize = FALSE
          )
        }
        
      } else {
        # ---------------- ElasticNet self-selection branch ----------------
        if (length(vars) == 0) {
          rec_main <- recipes::recipe(stats::as.formula(paste(outcome_var, "~ 1")), data = outer_d_train)
          fit_inner <- caret::train(
            rec_main, data = outer_d_train,
            method = "glm",
            metric = "logLoss", maximize = FALSE,
            trControl = ctrl_inner, family = "binomial"
          )
          
        } else if (length(vars) == 1) {
          # dummying is fine; still GLM
          rec_main <- make_recipe_after_elas(vars, outer_d_train, outcome_var)
          fit_inner <- caret::train(
            rec_main, data = outer_d_train,
            method = "glm",
            metric = "logLoss", maximize = FALSE,
            trControl = ctrl_inner, family = "binomial"
          )
          
        } else {
          rec_main <- make_recipe_after_elas(vars, outer_d_train, outcome_var)
          fit_inner <- caret::train(
            rec_main, data = outer_d_train,
            method = "glmnet",
            metric = "logLoss", maximize = FALSE,
            tuneGrid = tg, trControl = ctrl_inner,
            family = "binomial", standardize = FALSE
          )
        }
      }
      
      # Coefficients of final model
      if (fit_inner[["method"]] == "glmnet") {
        lam <- fit_inner$bestTune$lambda[1]
        coefs2 <- as.matrix(coef(fit_inner$finalModel, s = lam))
      } else {  # glm
        coefs2 <- coef(fit_inner$finalModel)
      }
      
      # Store coefficients of final model
      sel_vars_df$final_coefs[[i]][[j]] <- coefs2
      
      #---- Predictions on inner training (resamples) ----
      inner_preds <- fit_inner$pred
      p_in   <- inner_preds[[positive_class]]
      truth  <- inner_preds$obs
      
      inner_dat <- tibble::tibble(
        truth     = truth,
        truth_num = as.integer(truth == positive_class),
        pred_pos = p_in,
        hard_pred = factor(ifelse(p_in >= 0.5, positive_class, negative_class),
                           levels = levels(truth))
      )
      
      # ---- Metrics (yardstick with event = second; PRG for AUPRG) ----
      roc_auc       <- yardstick::roc_auc(inner_dat, truth, pred_pos, event_level = "second") |> dplyr::pull(.estimate)
      pr_auc        <- yardstick::pr_auc(inner_dat,  truth, pred_pos, event_level = "second") |> dplyr::pull(.estimate)
      prg_curve     <- prg::create_prg_curve(inner_dat$truth_num, inner_dat$pred_pos)
      auprg         <- prg::calc_auprg(prg_curve)
      logloss       <- yardstick::mn_log_loss(inner_dat, truth, pred_pos, event_level = "second") |> dplyr::pull(.estimate)
      if (fit_inner[["method"]] == "glmnet") {
        logloss_caret <- fit_inner$results |>
          dplyr::filter(alpha == fit_inner$bestTune$alpha,
                        lambda == fit_inner$bestTune$lambda) |>
          dplyr::pull(logLoss)
      } else {  # glm has no tuning grid; results has a single row
        logloss_caret <- fit_inner$results$logLoss[1]
      }
      brier = ModelMetrics::brier(actual = inner_dat$truth_num, predicted = inner_dat$pred_pos)
      mcc           <- yardstick::mcc(inner_dat, truth, hard_pred,  event_level = "second") |> dplyr::pull(.estimate)
      
      # Predict on outer test (using clamped probabilities)
      eps <- 1e-12
      p_out <- predict(fit_inner, newdata = outer_d_test, type = "prob")[[positive_class]]
      p_out <- pmin(pmax(p_out, eps), 1 - eps)
      
      # Calibration on outer test
      y01 <- as.integer(outer_d_test[[outcome_var]] == positive_class)
      logit_p <- qlogis(p_out)
      fit_cal  <- stats::glm(y01 ~ logit_p, family = binomial())
      cal_int  <- unname(coef(fit_cal)[1])
      cal_slope<- unname(coef(fit_cal)[2])
      fit_citl <- stats::glm(y01 ~ 1 + offset(logit_p), family = binomial())
      citl     <- unname(coef(fit_citl)[1])
      
      inner_perf[[i]][[j]] <- list(
        preds    = p_out,
        pred_idx = outer_test_idx,
        y        = outer_d_test[[outcome_var]],
        inner_preds_biased = p_in,
        inner_obs_biased = truth,
        inner_idx = outer_train_idx[inner_preds$rowIndex],
        metrics  = list(
          roc_auc = roc_auc,
          pr_auc  = pr_auc,
          auprg   = auprg,
          logloss = logloss,
          logloss_caret = logloss_caret,
          brier   = brier,
          mcc     = mcc,
          cal_intercept = cal_int,
          cal_slope     = cal_slope,
          citl_intercept= citl
        )
      )
    }
  } # End of outer loop
  # ---- Aggregate predictions over outer folds ----
  outer_df <- tibble::tibble()
  inner_df <- tibble::tibble()
  elas_once_df = tibble::tibble()
  
  for (i in seq_along(inner_perf)) {
    
    # ----- OUTER TEST FROM SINGLE ELASTIC PREDICTIONS -----
    elas_once_df <- dplyr::bind_rows(
      elas_once_df,
      tibble::tibble(
        idx     = elas_once[[i]]$pred_idx_elas_once,   # original row ids in test
        pred    = elas_once[[i]]$preds_elas_once,
        outcome = elas_once[[i]]$y_elas_once,
        baseline = elas_once[[i]]$baseline_pred_once
      )
    )
    
    for (j in seq_along(inner_perf[[i]])) {
      method <- if (j == 1) "VSURF" else "ElasticNet"
      
      # ----- OUTER TEST PREDICTIONS -----
      outer_df <- dplyr::bind_rows(
        outer_df,
        tibble::tibble(
          idx     = inner_perf[[i]][[j]]$pred_idx,   # original row ids in test
          pred    = inner_perf[[i]][[j]]$preds,
          outcome = inner_perf[[i]][[j]]$y,
          method  = method
        )
      )
      
      # ----- INNER (BIASED) PREDICTIONS -----
      inner_df <- dplyr::bind_rows(
        inner_df,
        tibble::tibble(
          inner_idx         = inner_perf[[i]][[j]]$inner_idx,          # original ids in train
          inner_biased_pred = inner_perf[[i]][[j]]$inner_preds_biased,
          inner_biased_obs  = inner_perf[[i]][[j]]$inner_obs_biased,
          method            = method
        )
      )
    }
  }
  
  # ----- Aggregate outer preds from single elastic-----
  avg_elas <- elas_once_df |>
    dplyr::group_by(idx) |>
    dplyr::summarise(avg_pred = mean(pred, na.rm = TRUE), .groups = "drop") |>
    dplyr::arrange(idx)
  
  outer_y_elas <- elas_once_df |>
    dplyr::group_by(idx) |>
    dplyr::summarise(y = dplyr::first(as.character(outcome)), .groups = "drop") |>
    dplyr::arrange(idx) |>
    dplyr::pull(y)
  
  # ----- Aggregate outer preds -----
  avg_wide <- outer_df |>
    dplyr::group_by(idx, method) |>
    dplyr::summarise(avg_pred = mean(pred, na.rm = TRUE), .groups = "drop") |>
    tidyr::pivot_wider(names_from = method, values_from = avg_pred) |>
    dplyr::arrange(idx)
  
  outer_y <- outer_df |>
    dplyr::group_by(idx) |>
    dplyr::summarise(y = dplyr::first(as.character(outcome)), .groups = "drop") |>
    dplyr::arrange(idx) |>
    dplyr::pull(y)
  
  # Aggregate prevalence (naive) model preds
  
  baseline_avg <- elas_once_df |>
    dplyr::group_by(idx) |>
    dplyr::summarise(avg_baseline = mean(baseline, na.rm = TRUE), .groups = "drop") |>
    dplyr::arrange(idx)
  
  # ----- Aggregate inner (biased) preds -----
  avg_inner_biased <- inner_df |>
    dplyr::filter(!is.na(inner_idx)) |>
    dplyr::group_by(inner_idx, method) |>
    dplyr::summarise(avg_pred = mean(inner_biased_pred, na.rm = TRUE), .groups = "drop") |>
    tidyr::pivot_wider(names_from = method, values_from = avg_pred) |>
    dplyr::arrange(inner_idx)
  
  inner_biased_obs <- inner_df |>
    dplyr::filter(!is.na(inner_idx)) |>
    dplyr::group_by(inner_idx) |>
    dplyr::summarise(y = dplyr::first(as.character(inner_biased_obs)), .groups = "drop") |>
    dplyr::arrange(inner_idx) |>
    dplyr::pull(y)
  
  return(list(
    idx_outer                          = avg_wide$idx,
    avg_final_outer_preds_vs           = avg_wide$VSURF,
    avg_final_outer_preds_elas         = avg_wide$ElasticNet,
    outer_y                            = outer_y,
    idx_inner                          = avg_inner_biased$inner_idx,
    avg_inner_biased_preds_vs          = avg_inner_biased$VSURF,
    avg_inner_biased_preds_elas        = avg_inner_biased$ElasticNet,
    inner_biased_obs                   = inner_biased_obs,
    idx_outer_elas                     = avg_elas$idx,
    avg_final_outer_preds_single_elas  = avg_elas$avg_pred,
    outer_y_single_elas                = outer_y_elas,
    avg_final_outer_preds_prev         = baseline_avg$avg_baseline,
    sel_vars_df                        = sel_vars_df,
    inner_perf                         = inner_perf
  ))
}

# ---- Inner resample summaries ----
inner_perf_nested_binary <- function(trained_object) {
  get_metrics <- function(m_list) {
    dplyr::bind_rows(lapply(m_list, function(x) as.data.frame(x$metrics)))
  }
  
  vs_metrics   <- get_metrics(lapply(trained_object$inner_perf, function(x) x[[1]]))
  elas_metrics <- get_metrics(lapply(trained_object$inner_perf, function(x) x[[2]]))
  
  summarise_block <- function(df) {
    dplyr::summarise(df,
                     roc_auc_mean = mean(roc_auc, na.rm = TRUE),       roc_auc_sd = sd(roc_auc, na.rm = TRUE),
                     pr_auc_mean  = mean(pr_auc,  na.rm = TRUE),       pr_auc_sd  = sd(pr_auc,  na.rm = TRUE),
                     auprg_mean   = mean(auprg,   na.rm = TRUE),       auprg_sd   = sd(auprg,   na.rm = TRUE),
                     logloss_mean = mean(logloss, na.rm = TRUE),       logloss_sd = sd(logloss, na.rm = TRUE),
                     logloss_caret_mean = mean(logloss_caret, na.rm = TRUE), logloss_caret_sd = sd(logloss_caret, na.rm = TRUE),
                     brier_mean   = mean(brier,   na.rm = TRUE),       brier_sd   = sd(brier,   na.rm = TRUE),
                     mcc_mean     = mean(mcc,     na.rm = TRUE),       mcc_sd     = sd(mcc,     na.rm = TRUE),
                     cal_int_mean = mean(cal_intercept, na.rm = TRUE), cal_int_sd = sd(cal_intercept, na.rm = TRUE),
                     cal_slope_mean=mean(cal_slope, na.rm = TRUE),     cal_slope_sd=sd(cal_slope, na.rm = TRUE),
                     citl_int_mean=mean(citl_intercept, na.rm = TRUE), citl_int_sd=sd(citl_intercept, na.rm = TRUE)
    )
  }
  
  vs_summary   <- summarise_block(vs_metrics)   |> dplyr::mutate(method = "VSURF")
  elas_summary <- summarise_block(elas_metrics) |> dplyr::mutate(method = "ElasticNet")
  
  dplyr::bind_rows(vs_summary, elas_summary) |>
    dplyr::select(method, dplyr::everything())
}


# ---- Outer performance on averaged predictions ----
outer_perf_nested_binary <- function(trained_object,
                                           positive_class = "Yes",
                                           negative_class = "No") {
  
  message(sprintf("IMPORTANT: All performance metrics are built around the predicted probability of the positive class: '%s'", positive_class))
  
  # ---- Variable selection stability ----
  stab_table <- function(sel_list) {
    tbl <- sort(table(unlist(sel_list)), decreasing = TRUE)
    as.data.frame(tbl, stringsAsFactors = FALSE) |>
      dplyr::rename(variable = Var1, times_selected = Freq) |>
      dplyr::mutate(freq = times_selected / nrow(trained_object$sel_vars_df))
  }
  
  cat("\n--- Elastic Net: selection stability ---\n")
  print(stab_table(trained_object$sel_vars_df$elastic_sel_vars))
  cat("\n--- VSURF: selection stability ---\n")
  print(stab_table(trained_object$sel_vars_df$vsurf_sel_vars))
  
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
    
    logit_p    <- qlogis(estimate)
    fit_cal    <- stats::glm(truth_num ~ logit_p, family = binomial())
    cal_int    <- unname(coef(fit_cal)[1])
    cal_slope  <- unname(coef(fit_cal)[2])
    fit_citl   <- stats::glm(truth_num ~ 1 + offset(logit_p), family = binomial())
    citl       <- unname(coef(fit_citl)[1])
    
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
  
  vs_metrics   <- get_metrics(trained_object$outer_y, trained_object$avg_final_outer_preds_vs)
  elas_metrics <- get_metrics(trained_object$outer_y, trained_object$avg_final_outer_preds_elas)
  inner_vs_metrics <- get_metrics(trained_object$inner_biased_obs, trained_object$avg_inner_biased_preds_vs)
  inner_elas_metrics <- get_metrics(trained_object$inner_biased_obs, trained_object$avg_inner_biased_preds_elas)
  outer_single_elas_metrics <- get_metrics(trained_object$outer_y_single_elas, trained_object$avg_final_outer_preds_single_elas)
  baseline_mod  <- get_metrics(trained_object$outer_y_single_elas, trained_object$avg_final_outer_preds_prev)
  
  vs_summary   <- dplyr::as_tibble(vs_metrics)   |> dplyr::mutate(method = "VSURF unbiased (outer)")
  elas_summary <- dplyr::as_tibble(elas_metrics) |> dplyr::mutate(method = "ElasticNet unbiased (outer)")
  base_summary <- dplyr::as_tibble(baseline_mod) |> dplyr::mutate(method = "Reference (prevalence pred) model (outer)")
  biased_vs_summary   <- dplyr::as_tibble(inner_vs_metrics)   |> dplyr::mutate(method = "VSURF biased (inner)")
  biased_elas_summary <- dplyr::as_tibble(inner_elas_metrics) |> dplyr::mutate(method = "ElasticNet biased (inner)")
  single_elas_summary <- dplyr::as_tibble(outer_single_elas_metrics) |> dplyr::mutate(method = "ElasticNet without var preselection (outer)")
  
  dplyr::bind_rows(vs_summary, elas_summary, single_elas_summary, base_summary, biased_vs_summary, biased_elas_summary) |>
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
    ntree = 1000,
    nforests = 20,
    inner_cv_method = c("LOOCV", "repeatedcv"),
    inner_cv_repeats = 20,
    inner_cv_folds = 5,
    selection_rule = c("best", "oneSE"),
    optim_metric = c("RMSE", "MAE")
) {
  # ---- Reproducible RNG ----
  RNGkind("L'Ecuyer-CMRG")
  set.seed(seed)
  cv_method <- match.arg(inner_cv_method)
  selection_rule <- match.arg(selection_rule)
  metric <- match.arg(optim_metric)
  
  
  # ---- Pre-checks ----
  pred_names <- setdiff(names(df), outcome_var)
  num_pred   <- vapply(df[pred_names], is.numeric, TRUE)
  
  if (any(vapply(df[pred_names][num_pred], function(x) any(x < 0, na.rm = TRUE), TRUE))) {
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
  
  # ---- Outer resampling indices ----
  cv_outer_train_folds_rows <-
    caret::createMultiFolds(df[[outcome_var]], k = cv_outer_folds, times = cv_outer_repeats)
  
  # ---- Recipes ----
  make_recipe_vsurf <- function(vars_to_keep, data, outcome_var) {
    all_pred <- setdiff(names(data), outcome_var)
    rec <- recipes::recipe(stats::reformulate(all_pred, response = outcome_var), data = data) |>
      recipes::step_string2factor(recipes::all_nominal_predictors()) |>
      recipes::step_impute_mode(recipes::all_nominal_predictors()) |>
      recipes::step_impute_knn(recipes::all_numeric_predictors(), neighbors = 5) |>
      recipes::step_log(recipes::all_numeric_predictors(), offset = 1) |>
      recipes::step_center(recipes::all_numeric_predictors()) |>
      recipes::step_scale(recipes::all_numeric_predictors())
    if (!is.null(vars_to_keep)) {
      rec <- rec |>
        recipes::step_rm(recipes::all_predictors(), -tidyselect::any_of(vars_to_keep))
    }
    rec
  }
  
  make_recipe_vsurf_glmnet <- function(vars_to_keep, data, outcome_var) {
    stopifnot(!is.null(vars_to_keep))
    make_recipe_vsurf(vars_to_keep, data, outcome_var) |>
      recipes::step_dummy(recipes::all_nominal_predictors(), one_hot = FALSE)
  }
  
  make_recipe_after_elas <- function(vars_to_keep, data, outcome_var) {
    all_pred <- setdiff(names(data), outcome_var)
    rec <- recipes::recipe(stats::reformulate(all_pred, response = outcome_var), data = data) |>
      recipes::step_string2factor(recipes::all_nominal_predictors()) |>
      recipes::step_impute_mode(recipes::all_nominal_predictors()) |>
      recipes::step_impute_knn(recipes::all_numeric_predictors(), neighbors = 5) |>
      recipes::step_log(recipes::all_numeric_predictors(), offset = 1) |>
      recipes::step_center(recipes::all_numeric_predictors()) |>
      recipes::step_scale(recipes::all_numeric_predictors()) |>
      recipes::step_dummy(recipes::all_nominal_predictors(), one_hot = FALSE)
    if (!is.null(vars_to_keep)) {
      rec <- rec |>
        recipes::step_rm(recipes::all_predictors(), -tidyselect::any_of(vars_to_keep))
    }
    rec
  }
  
  tg <- expand.grid(alpha = alpha_grid, lambda = lambda_grid)
  
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
    vsurf_sel_vars   = vector("list", n_outer),
    final_coefs      = vector("list", n_outer)
  )
  
  n_alg <- 2L
  inner_perf <- vector("list", n_outer)
  for (i in seq_len(n_outer)) {
    inner_perf[[i]] <- vector("list", n_alg)
    names(inner_perf[[i]]) <- c("VSURF", "ElasticNet")
  }
  
  elas_once <- vector("list", n_outer)
  
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
    rec_all_vsurf <- make_recipe_vsurf(NULL, outer_d_train, outcome_var)
    d_train_baked <- recipes::bake(recipes::prep(rec_all_vsurf, training = outer_d_train), new_data = NULL)
    x_train <- d_train_baked |>
      dplyr::select(-dplyr::all_of(outcome_var)) |>
      as.data.frame()
    y_train <- d_train_baked[[outcome_var]]
    
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
    sel_idx <- vs$varselect.interp
    vsurf_sel_vars <- if (length(sel_idx)) colnames(x_train)[sel_idx] else character(0)
    sel_vars_df$vsurf_sel_vars[[i]] <- vsurf_sel_vars
    message(sprintf("OUTER %02d — VSURF selected (%d): %s",
                    i, length(vsurf_sel_vars), paste(vsurf_sel_vars, collapse=", ")))
    
    # Elastic outer selection on ALL predictors (dummy'd recipe)
    rec_all_elas <- make_recipe_after_elas(NULL, outer_d_train, outcome_var)
    elastic_sel <- caret::train(
      rec_all_elas,
      data = outer_d_train,
      method = "glmnet",
      metric = metric,
      maximize = FALSE,
      tuneGrid = tg,
      trControl = ctrl_inner,
      family = "gaussian",
      standardize = FALSE
    )
    coefs <- as.matrix(coef(elastic_sel$finalModel,
                            s = elastic_sel$bestTune$lambda[1]))[, 1]
    el_sel_vars <- setdiff(names(coefs)[coefs != 0], "(Intercept)")
    sel_vars_df$elastic_sel_vars[[i]] <- el_sel_vars
    message(sprintf("OUTER %02d — Elastic selected (%d): %s",
                    i, length(el_sel_vars), paste(el_sel_vars, collapse=", ")))
    
    # One-pass elastic predictions on outer test
    p_elas <- predict(elastic_sel, newdata = outer_d_test)
    
    # Intercept-only baseline fitted on OUTER TRAIN
    rec_avg <- recipes::recipe(stats::as.formula(paste(outcome_var, "~ 1")), data = outer_d_train)
    
    fit_avg <- caret::train(
      rec_avg, data = outer_d_train,
      method = "glm",
      metric = metric, maximize = FALSE,
      trControl = caret::trainControl(method = "none"),  # no CV needed
      family = gaussian()
    )
    
    baseline_pred <- predict(fit_avg, newdata = outer_d_test)
    
    # Store per-row predictions (align to indices)
    elas_once[[i]] <- list(
      preds_elas_once    = p_elas,
      pred_idx_elas_once = outer_test_idx,
      y_elas_once        = outer_d_test[[outcome_var]],
      baseline_pred_once = baseline_pred
    )
    
    list_vars_both <- list(vsurf_sel_vars, el_sel_vars)
    sel_vars_df$final_coefs[[i]] <- vector("list", 2L)
    names(sel_vars_df$final_coefs[[i]]) <- c("VSURF", "ElasticNet")
    
    # ----- Inner loop (selected vars) -----
    for (j in seq_along(list_vars_both)) {
      vars <- list_vars_both[[j]]
      
      if (j == 1) {
        # VSURF branch
        if (length(vars) == 0) {
          rec_main <- recipes::recipe(stats::as.formula(paste(outcome_var, "~ 1")), data = outer_d_train)
          fit_inner <- caret::train(
            rec_main, data = outer_d_train,
            method = "glm",
            metric = metric, maximize = FALSE,
            trControl = ctrl_inner, family = "gaussian"
          )
        } else if (length(vars) == 1) {
          rec_main <- make_recipe_vsurf(vars, outer_d_train, outcome_var)
          fit_inner <- caret::train(
            rec_main, data = outer_d_train,
            method = "glm",
            metric = metric, maximize = FALSE,
            trControl = ctrl_inner, family = "gaussian"
          )
        } else {
          rec_main <- make_recipe_vsurf_glmnet(vars, outer_d_train, outcome_var)
          fit_inner <- caret::train(
            rec_main, data = outer_d_train,
            method = "glmnet",
            metric = metric, maximize = FALSE,
            tuneGrid = tg, trControl = ctrl_inner,
            family = "gaussian", standardize = FALSE
          )
        }
      } else {
        # Elastic self-selection branch
        if (length(vars) == 0) {
          rec_main <- recipes::recipe(stats::as.formula(paste(outcome_var, "~ 1")), data = outer_d_train)
          fit_inner <- caret::train(
            rec_main, data = outer_d_train,
            method = "glm",
            metric = metric, maximize = FALSE,
            trControl = ctrl_inner, family = "gaussian"
          )
        } else if (length(vars) == 1) {
          rec_main <- make_recipe_after_elas(vars, outer_d_train, outcome_var)
          fit_inner <- caret::train(
            rec_main, data = outer_d_train,
            method = "glm",
            metric = metric, maximize = FALSE,
            trControl = ctrl_inner, family = "gaussian"
          )
        } else {
          rec_main <- make_recipe_after_elas(vars, outer_d_train, outcome_var)
          fit_inner <- caret::train(
            rec_main, data = outer_d_train,
            method = "glmnet",
            metric = metric, maximize = FALSE,
            tuneGrid = tg, trControl = ctrl_inner,
            family = "gaussian", standardize = FALSE
          )
        }
      }
      
      # Coefficients of final model
      if (fit_inner[["method"]] == "glmnet") {
        lam <- fit_inner$bestTune$lambda[1]
        coefs2 <- as.matrix(coef(fit_inner$finalModel, s = lam))
      } else {
        coefs2 <- coef(fit_inner$finalModel)
      }
      sel_vars_df$final_coefs[[i]][[j]] <- coefs2
      
      # ---- Inner predictions & metrics ----
      inner_preds <- fit_inner$pred
      p_in   <- inner_preds$pred
      truth  <- inner_preds$obs
      
      inner_dat <- tibble::tibble(truth = truth, pred = p_in)
      
      rmse_val   <- yardstick::rmse(inner_dat, truth, pred, na_rm = TRUE) |> dplyr::pull(.estimate)
      mae_val    <- yardstick::mae(inner_dat,  truth, pred, na_rm = TRUE) |> dplyr::pull(.estimate)
      r2_val     <- yardstick::rsq(inner_dat,  truth, pred, na_rm = TRUE) |> dplyr::pull(.estimate)
      r2_tradval <- yardstick::rsq_trad(inner_dat, truth, pred, na_rm = TRUE) |> dplyr::pull(.estimate)
      ccc_val    <- yardstick::ccc(inner_dat, truth, pred, na_rm = TRUE) |> dplyr::pull(.estimate)
      
      # ---- Predict on outer test ----
      p_out <- predict(fit_inner, newdata = outer_d_test)
      
      # ---- Calibration (slope & intercept) ----
      y     <- outer_d_test[[outcome_var]]
      yhat  <- p_out
      fit_cal   <- stats::lm(y ~ yhat)
      cal_int   <- unname(coef(fit_cal)[1])
      cal_slope <- unname(coef(fit_cal)[2])
      fit_citl  <- stats::lm(y ~ 1 + offset(yhat))
      citl      <- unname(coef(fit_citl)[1])
      
      inner_perf[[i]][[j]] <- list(
        preds    = p_out,
        pred_idx = outer_test_idx,
        y        = y,
        inner_preds_biased = p_in,
        inner_obs_biased   = truth,
        inner_idx = outer_train_idx[inner_preds$rowIndex],
        metrics  = list(
          rmse = rmse_val,
          mae  = mae_val,
          r2   = r2_val,
          r2_trad = r2_tradval,
          ccc = ccc_val,
          cal_intercept = cal_int,
          cal_slope     = cal_slope,
          citl_intercept= citl
        )
      )
    }
  } # end outer
  
  # ---- Aggregate predictions over outer folds ----
  outer_df <- tibble::tibble()
  inner_df <- tibble::tibble()
  elas_once_df <- tibble::tibble()
  
  for (i in seq_along(inner_perf)) {
    # Single elastic outer preds
    elas_once_df <- dplyr::bind_rows(
      elas_once_df,
      tibble::tibble(
        idx     = elas_once[[i]]$pred_idx_elas_once,
        pred    = elas_once[[i]]$preds_elas_once,
        outcome = elas_once[[i]]$y_elas_once,
        baseline = elas_once[[i]]$baseline_pred_once
      )
    )
    
    for (j in seq_along(inner_perf[[i]])) {
      method <- if (j == 1) "VSURF" else "ElasticNet"
      
      outer_df <- dplyr::bind_rows(
        outer_df,
        tibble::tibble(
          idx     = inner_perf[[i]][[j]]$pred_idx,
          pred    = inner_perf[[i]][[j]]$preds,
          outcome = inner_perf[[i]][[j]]$y,
          method  = method
        )
      )
      inner_df <- dplyr::bind_rows(
        inner_df,
        tibble::tibble(
          inner_idx         = inner_perf[[i]][[j]]$inner_idx,
          inner_biased_pred = inner_perf[[i]][[j]]$inner_preds_biased,
          inner_biased_obs  = inner_perf[[i]][[j]]$inner_obs_biased,
          method            = method
        )
      )
    }
  }
  
  # Aggregate single elastic
  avg_elas <- elas_once_df |>
    dplyr::group_by(idx) |>
    dplyr::summarise(avg_pred = mean(pred, na.rm = TRUE), .groups = "drop") |>
    dplyr::arrange(idx)
  
  outer_y_elas <- elas_once_df |>
    dplyr::group_by(idx) |>
    dplyr::summarise(y = dplyr::first(outcome), .groups = "drop") |>
    dplyr::arrange(idx) |>
    dplyr::pull(y)
  
  # Aggregate outer preds
  avg_wide <- outer_df |>
    dplyr::group_by(idx, method) |>
    dplyr::summarise(avg_pred = mean(pred, na.rm = TRUE), .groups = "drop") |>
    tidyr::pivot_wider(names_from = method, values_from = avg_pred) |>
    dplyr::arrange(idx)
  
  outer_y <- outer_df |>
    dplyr::group_by(idx) |>
    dplyr::summarise(y = dplyr::first(outcome), .groups = "drop") |>
    dplyr::arrange(idx) |>
    dplyr::pull(y)
  
  # Aggregate average model preds
  
  baseline_avg <- elas_once_df |>
    dplyr::group_by(idx) |>
    dplyr::summarise(avg_baseline = mean(baseline, na.rm = TRUE), .groups = "drop") |>
    dplyr::arrange(idx)
  
  # Aggregate inner (biased) preds
  avg_inner_biased <- inner_df |>
    dplyr::filter(!is.na(inner_idx)) |>
    dplyr::group_by(inner_idx, method) |>
    dplyr::summarise(avg_pred = mean(inner_biased_pred, na.rm = TRUE), .groups = "drop") |>
    tidyr::pivot_wider(names_from = method, values_from = avg_pred) |>
    dplyr::arrange(inner_idx)
  
  inner_biased_obs <- inner_df |>
    dplyr::filter(!is.na(inner_idx)) |>
    dplyr::group_by(inner_idx) |>
    dplyr::summarise(y = dplyr::first(inner_biased_obs), .groups = "drop") |>
    dplyr::arrange(inner_idx) |>
    dplyr::pull(y)
  
  return(list(
    idx_outer                          = avg_wide$idx,
    avg_final_outer_preds_vs           = avg_wide$VSURF,
    avg_final_outer_preds_elas         = avg_wide$ElasticNet,
    outer_y                            = outer_y,
    idx_inner                          = avg_inner_biased$inner_idx,
    avg_inner_biased_preds_vs          = avg_inner_biased$VSURF,
    avg_inner_biased_preds_elas        = avg_inner_biased$ElasticNet,
    inner_biased_obs                   = inner_biased_obs,
    idx_outer_elas                     = avg_elas$idx,
    avg_final_outer_preds_single_elas  = avg_elas$avg_pred,
    outer_y_single_elas                = outer_y_elas,
    avg_final_outer_preds_baseline     = baseline_avg$avg_baseline,
    sel_vars_df                        = sel_vars_df,
    inner_perf                         = inner_perf
  ))
}


inner_perf_nested_continuous <- function(trained_object) {
  get_metrics <- function(m_list) {
    dplyr::bind_rows(lapply(m_list, function(x) as.data.frame(x$metrics)))
  }
  
  vs_metrics   <- get_metrics(lapply(trained_object$inner_perf, function(x) x[[1]]))
  elas_metrics <- get_metrics(lapply(trained_object$inner_perf, function(x) x[[2]]))
  
  summarise_block <- function(df) {
    dplyr::summarise(
      df,
      rmse_mean = mean(rmse, na.rm = TRUE),             rmse_sd = sd(rmse, na.rm = TRUE),
      mae_mean  = mean(mae,  na.rm = TRUE),             mae_sd  = sd(mae,  na.rm = TRUE),
      r2_mean   = mean(r2,   na.rm = TRUE),             r2_sd   = sd(r2,   na.rm = TRUE),
      r2_trad_mean = mean(r2_trad, na.rm = TRUE),       r2_trad_sd = sd(r2_trad, na.rm = TRUE),
      ccc_mean  = mean(ccc,  na.rm = TRUE),             ccc_sd  = sd(ccc,  na.rm = TRUE),
      cal_int_mean = mean(cal_intercept, na.rm = TRUE), cal_int_sd = sd(cal_intercept, na.rm = TRUE),
      cal_slope_mean = mean(cal_slope, na.rm = TRUE),   cal_slope_sd = sd(cal_slope, na.rm = TRUE),
      citl_int_mean = mean(citl_intercept, na.rm = TRUE), citl_int_sd = sd(citl_intercept, na.rm = TRUE)
    )
  }
  
  vs_summary   <- summarise_block(vs_metrics)   |> dplyr::mutate(method = "VSURF")
  elas_summary <- summarise_block(elas_metrics) |> dplyr::mutate(method = "ElasticNet")
  
  dplyr::bind_rows(vs_summary, elas_summary) |> dplyr::select(method, dplyr::everything())
}

outer_perf_nested_continuous <- function(trained_object) {
  stab_table <- function(sel_list, n_outer) {
    tbl <- sort(table(unlist(sel_list)), decreasing = TRUE)
    as.data.frame(tbl, stringsAsFactors = FALSE) |>
      dplyr::rename(variable = Var1, times_selected = Freq) |>
      dplyr::mutate(freq = times_selected / n_outer)
  }
  
  cat("\n--- Elastic Net: selection stability ---\n")
  print(stab_table(trained_object$sel_vars_df$elastic_sel_vars, nrow(trained_object$sel_vars_df)))
  cat("\n--- VSURF: selection stability ---\n")
  print(stab_table(trained_object$sel_vars_df$vsurf_sel_vars, nrow(trained_object$sel_vars_df)))
  
  get_metrics <- function(y_true, y_pred) {
    dat <- tibble::tibble(truth = y_true, pred = y_pred)
    rmse = yardstick::rmse(dat, truth, pred, na_rm = TRUE) |> dplyr::pull(.estimate)
    mae  = yardstick::mae(dat,  truth, pred, na_rm = TRUE) |> dplyr::pull(.estimate)
    r2   = yardstick::rsq(dat,  truth, pred, na_rm = TRUE) |> dplyr::pull(.estimate)
    r2_trad = yardstick::rsq_trad(dat, truth, pred, na_rm = TRUE) |> dplyr::pull(.estimate)
    ccc  = yardstick::ccc(dat,  truth, pred, na_rm = TRUE) |> dplyr::pull(.estimate)
    mape = yardstick::mape(dat,  truth, pred, na_rm = TRUE) |> dplyr::pull(.estimate)
    smape = yardstick::smape(dat, truth, pred, na_rm = TRUE) |> dplyr::pull(.estimate)
    
    fit_cal   <- stats::lm(y_true ~ y_pred)
    cal_int   <- unname(coef(fit_cal)[1])
    cal_slope <- unname(coef(fit_cal)[2])
    fit_citl  <- stats::lm(y_true ~ 1 + offset(y_pred))
    citl      <- unname(coef(fit_citl)[1])
    
    tibble::tibble(
      rmse = rmse, mae = mae, r2 = r2, r2_trad = r2_trad, ccc = ccc, mape = mape, smape = smape,
      cal_intercept = cal_int, cal_slope = cal_slope, citl_intercept = citl
    )
  }
  
  vs_metrics    <- get_metrics(trained_object$outer_y, trained_object$avg_final_outer_preds_vs)
  elas_metrics  <- get_metrics(trained_object$outer_y, trained_object$avg_final_outer_preds_elas)
  inner_vs      <- get_metrics(trained_object$inner_biased_obs, trained_object$avg_inner_biased_preds_vs)
  inner_elas    <- get_metrics(trained_object$inner_biased_obs, trained_object$avg_inner_biased_preds_elas)
  single_elas   <- get_metrics(trained_object$outer_y_single_elas, trained_object$avg_final_outer_preds_single_elas)
  baseline_mod  <- get_metrics(trained_object$outer_y_single_elas, trained_object$avg_final_outer_preds_baseline)
  
  dplyr::bind_rows(
    dplyr::as_tibble(vs_metrics)    |> dplyr::mutate(method = "VSURF unbiased (outer)"),
    dplyr::as_tibble(elas_metrics)  |> dplyr::mutate(method = "ElasticNet unbiased (outer)"),
    dplyr::as_tibble(single_elas)   |> dplyr::mutate(method = "ElasticNet no preselection (outer)"),
    dplyr::as_tibble(baseline_mod)  |> dplyr::mutate(method = "Reference (mean pred) model (outer)"),
    dplyr::as_tibble(inner_vs)      |> dplyr::mutate(method = "VSURF biased (inner)"),
    dplyr::as_tibble(inner_elas)    |> dplyr::mutate(method = "ElasticNet biased (inner)")
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
                                    ntree = 1000,
                                    nforests = 20,
                                    selection_rule = c("best", "oneSE"),
                                    cont_optim_metric = c("RMSE", "MAE")) {
  
  family <- match.arg(family)
  cv_method <- match.arg(cv_method)
  selection_rule <- match.arg(selection_rule)
  metric <- match.arg(cont_optim_metric)
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
  
  if (family == "binomial") {
    message(sprintf("IMPORTANT: All variables and coefficients are selected based on the predicted probability of the positive class: '%s'", positive_class))
    
    df[[outcome_var]] <- relevel(factor(df[[outcome_var]]), ref = negative_class)
    
    make_recipe_vsurf <- function(vars_to_keep, data, outcome_var) {
      all_pred <- setdiff(names(data), outcome_var)
      rec <- recipes::recipe(stats::reformulate(all_pred, response = outcome_var), data = data) |>
        recipes::step_string2factor(recipes::all_nominal_predictors()) |>
        recipes::step_impute_mode(recipes::all_nominal_predictors()) |>
        recipes::step_impute_knn(recipes::all_numeric_predictors(), neighbors = 5) |>
        recipes::step_log(recipes::all_numeric_predictors(), offset = 1) |>
        recipes::step_center(recipes::all_numeric_predictors()) |>
        recipes::step_scale(recipes::all_numeric_predictors())
      
      if (!is.null(vars_to_keep)) {
        rec <- rec |>
          recipes::step_rm(recipes::all_predictors(), -tidyselect::any_of(vars_to_keep))
      }
      rec
    }
    
    # Refit after VSURF for GLMNET: select RAW vars first, then dummy (k-1 coding)
    make_recipe_vsurf_glmnet <- function(vars_to_keep, data, outcome_var) {
      stopifnot(!is.null(vars_to_keep))  # only for selected-var refits
      rec <- make_recipe_vsurf(vars_to_keep, data, outcome_var) |>
        recipes::step_dummy(recipes::all_nominal_predictors(), one_hot = FALSE)
      rec
    }
    
    # ElasticNet self-selection: dummy FIRST → then (optionally) select dummy names
    make_recipe_after_elas <- function(vars_to_keep, data, outcome_var) {
      all_pred <- setdiff(names(data), outcome_var)
      rec <- recipes::recipe(stats::reformulate(all_pred, response = outcome_var), data = data) |>
        recipes::step_string2factor(recipes::all_nominal_predictors()) |>
        recipes::step_impute_mode(recipes::all_nominal_predictors()) |>
        recipes::step_impute_knn(recipes::all_numeric_predictors(), neighbors = 5) |>
        recipes::step_log(recipes::all_numeric_predictors(), offset = 1) |>
        recipes::step_center(recipes::all_numeric_predictors()) |>
        recipes::step_scale(recipes::all_numeric_predictors()) |>
        recipes::step_dummy(recipes::all_nominal_predictors(), one_hot = FALSE)
      
      if (!is.null(vars_to_keep)) {
        rec <- rec |>
          recipes::step_rm(recipes::all_predictors(), -tidyselect::any_of(vars_to_keep))
      }
      rec
    }
    
    rec_all <- make_recipe_vsurf(NULL, df, outcome_var)
    d_train_baked <- recipes::bake(recipes::prep(rec_all, training = df), new_data = NULL)
    
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
    
    # ---- Elastic Net ----
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
    
    tg <- expand.grid(
      alpha  = alpha_grid,
      lambda = lambda_grid
    )
    
    #---- Elastic-net selection on all predictors ----
    rec_all <- make_recipe_after_elas(NULL, df, outcome_var)
    
    elas_pre <- caret::train(
      rec_all,
      data = df,
      method = "glmnet",
      metric = "logLoss",
      maximize = FALSE,
      tuneGrid = tg,
      trControl = ctrl_inner,
      family = "binomial",
      standardize = FALSE
    )
    
    #Coef
    coefs <- as.matrix(coef(elas_pre$finalModel,
                            s = elas_pre$bestTune$lambda[1]))[, 1]
    el_sel_vars <- setdiff(names(coefs)[coefs != 0], "(Intercept)")
    coef_elas_only <- coefs[coefs != 0]
    message(sprintf("ElasticNet selected (%d): %s",
                    length(el_sel_vars), paste(el_sel_vars, collapse=", ")))
    message("ElasticNet coefficients (all predictors):")
    print(coef_elas_only)
    
    #---- Final model ----
    sel_vars_df <- tibble::tibble(
      elastic_single_coef      = list(coef_elas_only),
      vsurf_sel_vars_int_step  = list(vsurf_sel_vars_int),
      vsurf_sel_vars_pred_step = list(vsurf_sel_vars_pred),
      final_coefs              = list(list(VSURF = NULL, ElasticNet = NULL)) 
    )
    
    final_fits <- list()
    
    for (v in c("VSURF", "ElasticNet")) {
      message(sprintf("\nFitting final elastic-net model using %s-selected variables...", v))
      if (v == "VSURF") {
        sel_vars <- vsurf_sel_vars_int
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
        
        final_fits[[v]] <- final_fit
        coef_final <- coef(final_fit$finalModel)
        sel_vars_df$final_coefs[[1]][[v]] <- coef_final
        
      } else if (length(sel_vars) == 1) {
        rec_main <- make_recipe_vsurf(sel_vars, df, outcome_var)
        final_fit <- caret::train(
          rec_main, data = df,
          method = "glm",
          metric = "logLoss", maximize = FALSE,
          trControl = caret::trainControl(method = "none"),  # no CV needed
          family = "binomial"
        )
        
        final_fits[[v]] <- final_fit
        coef_final <- coef(final_fit$finalModel)
        sel_vars_df$final_coefs[[1]][[v]] <- coef_final
        
      } else {
        if (v == "VSURF") {
          rec_main <- make_recipe_vsurf_glmnet(sel_vars, df, outcome_var)
        } else {
          rec_main <- make_recipe_after_elas(sel_vars, df, outcome_var)
        }
        
        final_fit <- caret::train(
          rec_main, data = df,
          method = "glmnet",
          metric = "logLoss", maximize = FALSE,
          tuneGrid = tg,
          trControl = ctrl_inner,  
          family = "binomial",
          standardize = FALSE
        )
        
        final_fits[[v]] <- final_fit
        coefs_2 <- as.matrix(coef(final_fit$finalModel,
                                  s = final_fit$bestTune$lambda[1]))[, 1]
        sel_vars_df$final_coefs[[1]][[v]] <- coefs_2[coefs_2 != 0]
      }
      message(sprintf("Final model coefficients (after %s preselection):\n", v))
      print(sel_vars_df$final_coefs[[1]][[v]])
    }
    return(list(
      final_model = final_fits,
      sel_vars_df = sel_vars_df,
      vsurf_object = vs,
      elasticnet_object = elas_pre
    ))
  }
  else if (family == "gaussian") {

    make_recipe_vsurf <- function(vars_to_keep, data, outcome_var) {
      all_pred <- setdiff(names(data), outcome_var)
      rec <- recipes::recipe(stats::reformulate(all_pred, response = outcome_var), data = data) |>
        recipes::step_string2factor(recipes::all_nominal_predictors()) |>
        recipes::step_impute_mode(recipes::all_nominal_predictors()) |>
        recipes::step_impute_knn(recipes::all_numeric_predictors(), neighbors = 5) |>
        recipes::step_log(recipes::all_numeric_predictors(), offset = 1) |>
        recipes::step_center(recipes::all_numeric_predictors()) |>
        recipes::step_scale(recipes::all_numeric_predictors())
      
      if (!is.null(vars_to_keep)) {
        rec <- rec |>
          recipes::step_rm(recipes::all_predictors(), -tidyselect::any_of(vars_to_keep))
      }
      rec
    }
    
    # Refit after VSURF for GLMNET: select RAW vars first, then dummy (k-1 coding)
    make_recipe_vsurf_glmnet <- function(vars_to_keep, data, outcome_var) {
      stopifnot(!is.null(vars_to_keep))  # only for selected-var refits
      rec <- make_recipe_vsurf(vars_to_keep, data, outcome_var) |>
        recipes::step_dummy(recipes::all_nominal_predictors(), one_hot = FALSE)
      rec
    }
    
    # ElasticNet self-selection: dummy FIRST → then (optionally) select dummy names
    make_recipe_after_elas <- function(vars_to_keep, data, outcome_var) {
      all_pred <- setdiff(names(data), outcome_var)
      rec <- recipes::recipe(stats::reformulate(all_pred, response = outcome_var), data = data) |>
        recipes::step_string2factor(recipes::all_nominal_predictors()) |>
        recipes::step_impute_mode(recipes::all_nominal_predictors()) |>
        recipes::step_impute_knn(recipes::all_numeric_predictors(), neighbors = 5) |>
        recipes::step_log(recipes::all_numeric_predictors(), offset = 1) |>
        recipes::step_center(recipes::all_numeric_predictors()) |>
        recipes::step_scale(recipes::all_numeric_predictors()) |>
        recipes::step_dummy(recipes::all_nominal_predictors(), one_hot = FALSE)
      
      if (!is.null(vars_to_keep)) {
        rec <- rec |>
          recipes::step_rm(recipes::all_predictors(), -tidyselect::any_of(vars_to_keep))
      }
      rec
    }
    
    rec_all <- make_recipe_vsurf(NULL, df, outcome_var)
    d_train_baked <- recipes::bake(recipes::prep(rec_all, training = df), new_data = NULL)
    
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
    
    # ---- Elastic Net ----
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
    
    tg <- expand.grid(
      alpha  = alpha_grid,
      lambda = lambda_grid
    )
    
    #---- Elastic-net selection on all predictors ----
    rec_all <- make_recipe_after_elas(NULL, df, outcome_var)
    
    elas_pre <- caret::train(
      rec_all,
      data = df,
      method = "glmnet",
      metric = metric,
      maximize = FALSE,
      tuneGrid = tg,
      trControl = ctrl_inner,
      family = "gaussian",
      standardize = FALSE
    )
    
    #Coef
    coefs <- as.matrix(coef(elas_pre$finalModel,
                            s = elas_pre$bestTune$lambda[1]))[, 1]
    el_sel_vars <- setdiff(names(coefs)[coefs != 0], "(Intercept)")
    coef_elas_only <- coefs[coefs != 0]
    message(sprintf("ElasticNet selected (%d): %s",
                    length(el_sel_vars), paste(el_sel_vars, collapse=", ")))
    message("ElasticNet coefficients (all predictors):")
    print(coef_elas_only)
    
    #---- Final model ----
    sel_vars_df <- tibble::tibble(
      elastic_single_coef      = list(coef_elas_only),
      vsurf_sel_vars_int_step  = list(vsurf_sel_vars_int),
      vsurf_sel_vars_pred_step = list(vsurf_sel_vars_pred),
      final_coefs              = list(list(VSURF = NULL, ElasticNet = NULL))  
    )
    
    final_fits <- list()

    for (v in c("VSURF", "ElasticNet")) {
      message(sprintf("\nFitting final elastic-net model using %s-selected variables...", v))
      if (v == "VSURF") {
        sel_vars <- vsurf_sel_vars_int
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
        
        final_fits[[v]] <- final_fit
        coef_final <- coef(final_fit$finalModel)
        sel_vars_df$final_coefs[[1]][[v]] <- coef_final
        
      } else if (length(sel_vars) == 1) {
        rec_main <- make_recipe_vsurf(sel_vars, df, outcome_var)
        final_fit <- caret::train(
          rec_main, data = df,
          method = "glm",
          metric = metric, maximize = FALSE,
          trControl = caret::trainControl(method = "none"),  # no CV needed
          family = "gaussian"
        )
        
        final_fits[[v]] <- final_fit
        coef_final <- coef(final_fit$finalModel)
        sel_vars_df$final_coefs[[1]][[v]] <- coef_final
        
      } else {
        if (v == "VSURF") {
          rec_main <- make_recipe_vsurf_glmnet(sel_vars, df, outcome_var)
        } else {
          rec_main <- make_recipe_after_elas(sel_vars, df, outcome_var)
        }
        
        final_fit <- caret::train(
          rec_main, data = df,
          method = "glmnet",
          metric = metric, maximize = FALSE,
          tuneGrid = tg,
          trControl = ctrl_inner,  
          family = "gaussian",
          standardize = FALSE
        )
        
        final_fits[[v]] <- final_fit
        coefs_2 <- as.matrix(coef(final_fit$finalModel,
                                  s = final_fit$bestTune$lambda[1]))[, 1]
        sel_vars_df$final_coefs[[1]][[v]] <- coefs_2[coefs_2 != 0]
      }
      message(sprintf("Final model coefficients (after %s preselection):\n", v))
      print(sel_vars_df$final_coefs[[1]][[v]])
    }
    return(list(
      final_model = final_fits,
      sel_vars_df = sel_vars_df,
      vsurf_object = vs,
      elasticnet_object = elas_pre
    ))
  } else {
    stop("family must be either 'binomial' or 'gaussian'.")
  }
}
                          
