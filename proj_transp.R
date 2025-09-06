library(caret)
library(tidyverse)
library(ModelMetrics)
library(recipes)
library(glmnet)
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
library(VIM)
library(janitor)
library(report)
set.seed(123)


# Read xlsx and recode the outcome
df <- read_xlsx("D:/Packages/proj_transpl/Cytokines_heart_data_24.06 1 and 4 groups.xlsx",
                sheet = 1)

df_1 = df %>%
  clean_names() %>%
  mutate(
    outcome = pull(df[,2]),
    outcome = ifelse(outcome == 1, 1, 0),
    patient_id= row_number() # Create a unique patient ID
  ) %>%
  mutate(outcome = factor(outcome, levels = c(0, 1), labels = c("No", "Yes")),
         sex = factor(sex_1_m_2_f, levels = c(1, 2), labels = c("Male", "Female")),
         time_tr = time_from_transplantation_month
         ) %>%
  select(patient_id, outcome, sex, age,time_tr, c(-1,-2,-5, -time_from_transplantation_month))

df_1 %>%
  is.na() %>%
  colSums() #NA in 3 variables: CRP,Troponin,NT-proBNP

###########################KNN imputation
# Sorting variables in the data frame by the number of NAs
vars_by_NAs <- df_1 %>%
  is.na() %>%
  colSums() %>%
  sort(decreasing = FALSE) %>% 
  names()

# Imputation with kNN, k=5, using weighted mean for numeric variables
df_imp_0<- df_1 %>% 
  select(all_of(vars_by_NAs)) %>% 
  kNN(k=5,
      numFun = weighted.mean,
      weightDist = TRUE) %>%
  select(1:64)

df_imp_0 %>%
  is.na() %>%
  colSums() #0 NA values


##################################################' Find near zero variance predictors
zero_var = df_imp_0 %>%
  caret::nearZeroVar(saveMetrics = TRUE) %>%
  filter(percentUnique < 50) #Variables egf_13, il_17e_il_25_51, il_17f_53 , il_22_55 have very low percent unique values => remove
  
df_imp_1 = df_imp_0 %>%
  select(-c(egf_13, il_17e_il_25_51, il_17f_53, il_22_55,time_tr,patient_id))#%>%
  #slice(-c(6))  
#dont delete observation 2 because this is the only one female among cases and you will get the separation issue,
#however this observation is an outlier(low value of if while being rejected), so you can keep it for sensitivity analysis

boxplot(df_imp_1,outline=FALSE,col="cornflowerblue")
##################################################Correlation analysis

cor_mat=cor(df_imp_1[,4:length(df_imp_1)], method = "pearson")

#' Find correlated variable pairs in a correlation matrix

find_correlated <- function(cor_mat, cutoff = 0.70) {
  
  if (!is.matrix(cor_mat) || !is.numeric(cor_mat))
    stop("cor_mat must be a numeric matrix (e.g. the output of cor()).")
  
  if (!identical(rownames(cor_mat), colnames(cor_mat)))
    stop("Row and column names must match (a square correlation matrix).")
  
  ## Keep only the upper triangle so we don’t report duplicates or self‑correlations
  cor_upper <- cor_mat
  cor_upper[lower.tri(cor_upper, diag = TRUE)] <- NA
  
  ## Indices of cells that meet the cutoff
  idx <- which(abs(cor_upper) >= cutoff, arr.ind = TRUE)
  
  ## Build a tidy data frame
  data.frame(
    var1 = rownames(cor_mat)[idx[,"row"]],
    var2 = colnames(cor_mat)[idx[,"col"]],
    r    = cor_upper[idx],
    row.names = NULL
  )
}

high_cor_var=find_correlated(cor_mat)

#heatmap
pheatmap(
  cor_mat,
  cluster_rows = TRUE,    # whether to cluster rows
  cluster_cols = TRUE,    # whether to cluster columns
  display_numbers = TRUE, # show correlation values in each cell
  number_format = "%.2f", # format correlation values to 2 decimals
  color = colorRampPalette(c("blue", "white", "red"))(50),
  main = "Correlation Heatmap"
)

#############################################################################################################VSURF(NON STANDARDIZED DATA)

x_train = df_imp_1 %>%
  select(2:length(df_imp_1)) %>%
  mutate(across(where(is.numeric), ~ log1p(.)))  #convert all numeric vars to numeric

y_train <- df_imp_1 %>% 
  pull(outcome) 

set.seed(999)
vsurf_model = VSURF(y=y_train, x=x_train, ntree.thres = 10000,nfor.thres = 100, 
                    ntree.interp = 10000, nfor.interp=100, 
                    ntree.pred = 10000, nfor.pred = 100,
                    RFimplem = "randomForest", parallel = TRUE)
summary(vsurf_model)
plot(vsurf_model)

plot(vsurf_model, step="thres", imp.sd=F, main="Variable importance plot", var.names = T) # variable importance plot
#green line - a smoothing function using CART
#red line - threshold which is a minimum predicted value from green line (everything below threshold is rejected for further steps)
colnames(x_train[vsurf_model[["varselect.thres"]]])
vsurf_model[["imp.mean.dec"]]*100 #var importance

plot(vsurf_model, step="interp", main="Variable importance plot", var.names = T)
#red line - threshold, the smallest model with OOB less than lowest OOB+ 1 sd
colnames(x_train[vsurf_model[["varselect.interp"]]]) #"if_ng_26"     "mip_1b_67"    "tn_fa_76"     "il_12_p40_43" "time_tr"      "tn_fb_77" - the order of importance.
vsurf_model[["err.interp"]] # OOB error RATE for the model with selected variables after interp step. 
#0.5 error rate means that the model is not better than random guess (=variables in the dataset are random). 
#Thats because OOB error rate for classification(and regression if values are from 0 to 1 from simulation) is [number of wrong classifications/ n] can be from 0 to 1

plot(vsurf_model, step="pred", imp.mean=FALSE, main="Variable importance plot")
colnames(x_train[vsurf_model[["varselect.pred"]]])

ggplot(df_imp_1, aes(x=if_ng_26)) + 
  geom_histogram()

ggplot(df_imp_1, aes(x=il_12_p40_43)) + 
  geom_histogram()

ggplot(df_imp_1, aes(x=tn_fb_77)) + 
  geom_histogram() #different scale

ggplot(df_imp_1, aes(x=mip_1b_67)) + 
  geom_histogram()

ggplot(df_imp_1, aes(x=tn_fa_76)) + 
  geom_histogram()

df_imp_1 %>%
  select(if_ng_26, il_12_p40_43, tn_fb_77, mip_1b_67, tn_fa_76) %>%
  preProcess(
    method = c("expoTrans", "center", "scale"),
    na.remove = TRUE
  ) %>%
  cor()  #transformation is desirable


rf = cforest(outcome ~ if_ng_26 , data = df_imp_1)
estimates(rf)
############################################################################################## Elastic Net Regression with caret

#All numeric columns should be >= 0 for log transform inside the function

nested_elastic_binary_outcome <- function(
    df,
    outcome_var = "outcome",
    positive_class = "Yes",
    negative_class = "No",
    cv_outer_folds = 5,
    cv_outer_repeats = 20,
    seed = 1,
    alpha_grid  = seq(0, 1, by = 0.1),
    lambda_grid = 10^seq(-4, 1, length.out = 50)
) {
  # ---- Reproducible RNG ----
  RNGkind("L'Ecuyer-CMRG")
  set.seed(seed)
  
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
  make_recipe <- function(vars_to_keep, data, outcome_var) {
    # Use all raw predictors in the formula
    all_pred <- setdiff(names(data), outcome_var)
    
    rec <- recipes::recipe(stats::reformulate(all_pred, response = outcome_var), data = data) |>
      recipes::step_string2factor(recipes::all_nominal_predictors()) |>
      recipes::step_impute_mode(recipes::all_nominal_predictors()) |>
      recipes::step_impute_knn(recipes::all_numeric_predictors(), neighbors = 5) |>
      recipes::step_log(recipes::all_numeric_predictors(), offset = 1) |>
      recipes::step_dummy(recipes::all_nominal_predictors(), one_hot = FALSE) |>
      recipes::step_center(recipes::all_numeric_predictors()) |>
      recipes::step_scale(recipes::all_numeric_predictors()) 
    
    
    # After dummying, optionally keep only the selected post-dummy columns
    if (!is.null(vars_to_keep)) {
      rec <- rec |>
        recipes::step_select(recipes::all_outcomes(), tidyselect::any_of(vars_to_keep))
    }
    
    return(rec)
  }
  
  tg <- expand.grid(alpha = alpha_grid, lambda = lambda_grid)
  
  # ---- Inner CV control (optimizes logLoss) ----
  ctrl_inner <- caret::trainControl(
    method = "LOOCV",
    classProbs = TRUE,
    savePredictions = "final",
    summaryFunction = caret::mnLogLoss,
    verboseIter = FALSE,
    allowParallel = FALSE
  )
  
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
    message(sprintf("Outer %d class counts: %s",
                    i, paste(sprintf("%s=%d", names(tab), tab), collapse=", ")))
    
    # Keep negative first consistently
    outer_d_train[[outcome_var]] <- factor(outer_d_train[[outcome_var]],
                                           levels = c(negative_class, positive_class))
    outer_d_test[[outcome_var]]  <- factor(outer_d_test[[outcome_var]],
                                           levels = c(negative_class, positive_class))
    
    # Bake training data for VSURF
    rec_all <- make_recipe(NULL, outer_d_train, outcome_var)
    d_train_baked <- recipes::bake(recipes::prep(rec_all, training = outer_d_train), new_data = NULL)
    
    x_train <- d_train_baked |>
      dplyr::select(-dplyr::all_of(outcome_var)) |>
      as.data.frame()
    
    y_train <- d_train_baked[[outcome_var]]
    
    # ---- VSURF ----
    vs <- VSURF::VSURF(
      y = y_train,
      x = x_train,
      ntree.thres   = 5000, nfor.thres  = 50,
      ntree.interp  = 5000, nfor.interp = 50,
      ntree.pred    = 5000, nfor.pred   = 50,
      RFimplem = "ranger",
      parallel = TRUE
    )
    
    sel_idx <- vs$varselect.interp
    vsurf_sel_vars <- if (length(sel_idx)) colnames(x_train)[sel_idx] else character(0)
    sel_vars_df$vsurf_sel_vars[[i]] <- vsurf_sel_vars
    message(sprintf("OUTER %02d — VSURF selected (%d): %s",
                    i, length(vsurf_sel_vars), paste(vsurf_sel_vars, collapse=", ")))
    
    # ---- Elastic-net selection on all predictors ----
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
    
    elas_once[[i]] <- list(
      preds_elas_once    = p_elas,
      pred_idx_elas_once = outer_test_idx,
      y_elas_once        = outer_d_test[[outcome_var]]
    )
    
    list_vars_both <- list(vsurf_sel_vars, el_sel_vars)
    
    # Store coefficients of final models for both methods
    sel_vars_df$final_coefs[[i]] <- vector("list", 2L)
    names(sel_vars_df$final_coefs[[i]]) <- c("VSURF", "ElasticNet")
    
    #---------------------- Inner CV with selected vars -------------------
    for (j in seq_along(list_vars_both)) {
      vars <- list_vars_both[[j]]
      rec_main <- make_recipe(vars, outer_d_train, outcome_var)
      
      if (length(vars) == 0) {
        # Build an intercept-only recipe
        rec_main <- recipes::recipe(stats::as.formula(paste(outcome_var, "~ 1")), data = outer_d_train)
        fit_inner <- caret::train(
          rec_main, data = outer_d_train,
          method = "glm",
          metric = "logLoss",
          maximize = FALSE,
          trControl = ctrl_inner,
          family = "binomial"
        )
      } else if (length(vars) == 1) {
        # 'glm' path
        fit_inner <- caret::train(
          rec_main,
          data = outer_d_train,
          method = "glm",
          metric = "logLoss",
          maximize = FALSE,
          trControl = ctrl_inner,
          family = "binomial"
        )
      } else {
        #glmnet path
        fit_inner <- caret::train(
          rec_main,
          data = outer_d_train,
          method = "glmnet",
          metric = "logLoss",
          maximize = FALSE,
          tuneGrid = tg,
          trControl = ctrl_inner,
          family = "binomial",
          standardize = FALSE
        )
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
        outcome = elas_once[[i]]$y_elas_once
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
    sel_vars_df                        = sel_vars_df,
    inner_perf                         = inner_perf
  ))
}


results = nested_elastic_binary_outcome(
  df_imp_1, 
  outcome_var = "outcome",
  positive_class = "Yes",
  negative_class = "No",
  cv_outer_folds = 5,
  cv_outer_repeats = 20,
  seed = 1,
  alpha_grid = seq(0, 1, by = 0.1),
  lambda_grid = 10^seq(-4, 1, length.out = 50))

# ---- Inner resample summaries ----
inner_model_perf_nested_binary <- function(trained_object) {
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

inner_res <- inner_model_perf_nested_binary(trained_object = results)


# ---- Outer performance on averaged predictions ----
outer_model_perf_nested_binary <- function(trained_object,
                                           positive_class = "Yes",
                                           negative_class = "No") {
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
  
  vs_summary   <- dplyr::as_tibble(vs_metrics)   |> dplyr::mutate(method = "VSURF unbiased (outer)")
  elas_summary <- dplyr::as_tibble(elas_metrics) |> dplyr::mutate(method = "ElasticNet unbiased (outer)")
  biased_vs_summary   <- dplyr::as_tibble(inner_vs_metrics)   |> dplyr::mutate(method = "VSURF biased (inner)")
  biased_elas_summary <- dplyr::as_tibble(inner_elas_metrics) |> dplyr::mutate(method = "ElasticNet biased (inner)")
  single_elas_summary <- dplyr::as_tibble(outer_single_elas_metrics) |> dplyr::mutate(method = "ElasticNet without var preselection (outer)")
  
  dplyr::bind_rows(vs_summary, elas_summary, single_elas_summary, biased_vs_summary, biased_elas_summary) |>
    dplyr::select(method, dplyr::everything())
}

outer_res <- outer_model_perf_nested_binary(trained_object = results,
                                            positive_class = "Yes",
                                            negative_class = "No")

#AUC ROC-the probability is that a randomly chosen cancer patient is ranked higher (given higher prob) than a randomly chosen healthy patient 
#the AUC judges a correctly assigned biopsy (TP) and an unnecessary biopsy (FP) as equally important.
# ---- Calibration plots (rms) ----
y01 <- as.integer(results$outer_y == "Yes")  # or pass in positive_class if not "Yes"

cal_plot_vs <- val.prob.ci.2(
  results$avg_final_outer_preds_vs,
  y01,
  logistic = TRUE,
  col.log = "blue",
  CL.smooth = FALSE
)

cal_plot_elas <- val.prob.ci.2(
  results$avg_final_outer_preds_elas,
  y01,
  logistic = TRUE,
  col.log = "blue",
  CL.smooth = FALSE
)

cal_plot_elas <- val.prob.ci.2(
  results$avg_final_outer_preds_single_elas,
  y01,
  logistic = TRUE,
  col.log = "blue",
  CL.smooth = FALSE
)

cal_plot_elas <- val.prob.ci.2(
  results$avg_inner_biased_preds_elas,
  y01,
  logistic = TRUE,
  col.log = "blue",
  CL.smooth = FALSE
)

cal_plot_elas <- val.prob.ci.2(
  results$avg_inner_biased_preds_vs,
  y01,
  logistic = TRUE,
  col.log = "blue",
  CL.smooth = FALSE
)

# Optional extra check
val.prob(results$avg_final_outer_preds_vs,   y01)
val.prob(results$avg_final_outer_preds_elas, y01)


# ---- Decision curve analysis (dcurves) ----
res_df <- tibble::tibble(
  pred_elas = results$avg_final_outer_preds_elas,
  pred_vs   = results$avg_final_outer_preds_vs,
  pred_single_elas = results$avg_final_outer_preds_single_elas,
  y01       = y01,
  if_ng_26 = df_imp_1$if_ng_26[results$idx_outer],
  sex = df_imp_1$sex[results$idx_outer]
)

dcurves::dca(
  y01 ~ pred_elas + pred_vs +pred_single_elas,
  data = res_df,
  thresholds = seq(0, 0.6, by = 0.01)
) %>%
  plot(smooth = T)

#The value of 0.16 at a threshold probability of 20% can be interpreted as: “Comparing to conducting no treatment(intervention), 
#intervention on the basis of the elastic model is the equivalent of a strategy that found 16 rejections per hundred patients without conducting any unnecessary interventions”

dca(y01 ~ pred_elas + pred_vs +pred_single_elas,
    data = res_df,
    thresholds = seq(0.05, 1, 0.01)
) %>%
  net_intervention_avoided() %>%
  plot(smooth = T)

#At a probability threshold of 15-30%, the net reduction in interventions is about 0.25. 
#In other words, at this probability threshold, biopsying patients on the basis of the model is the equivalent of
#a strategy that led to an absolute 25% reduction in the number of biopsies without missing any heart rejections (25 unnecessary biopsies avoided per 100 patients). 

dca(y01 ~ pred_elas + pred_vs +pred_single_elas,
    data = res_df,
    thresholds = seq(0.1, 0.15, 0.2)
) %>%
  as_tibble() %>%
  select(label, threshold, net_benefit) %>%
  gt::gt() %>%
  gt::fmt_percent(columns = threshold, decimals = 0) %>%
  gt::cols_label(
    label = "Strategy",
    threshold = "Decision Threshold",
    net_benefit = "Net Benefit"
  ) %>%
  gt::cols_align("left", columns = label)





################################################################################Final model#######################################################################
set.seed(1)
outer_d_train = df_imp_1
outcome_var = "outcome"
positive_class = "Yes"
negative_class = "No"
outer_d_train[[outcome_var]] <- relevel(factor(outer_d_train[[outcome_var]]), ref = "No")


make_recipe <- function(vars_to_keep, data, outcome_var) {
  # Use all raw predictors in the formula
  all_pred <- setdiff(names(data), outcome_var)
  
  rec <- recipes::recipe(stats::reformulate(all_pred, response = outcome_var), data = data) |>
    recipes::step_string2factor(recipes::all_nominal_predictors()) |>
    recipes::step_impute_mode(recipes::all_nominal_predictors()) |>
    recipes::step_impute_knn(recipes::all_numeric_predictors(), neighbors = 5) |>
    recipes::step_log(recipes::all_numeric_predictors(), offset = 1) |>
    recipes::step_dummy(recipes::all_nominal_predictors(), one_hot = FALSE) |>
    recipes::step_center(recipes::all_numeric_predictors()) |>
    recipes::step_scale(recipes::all_numeric_predictors()) 
  
  # After dummying, optionally keep only the selected post-dummy columns
  if (!is.null(vars_to_keep) && length(vars_to_keep) > 0) {
    rec <- rec |>
      step_select(recipes::all_outcomes(), tidyselect::any_of(vars_to_keep))
  }
  
  return(rec)
}

rec_all <- make_recipe(NULL, outer_d_train, outcome_var)
d_train_baked <- recipes::bake(recipes::prep(rec_all, training = outer_d_train), new_data = NULL)

x_train <- d_train_baked |>
  dplyr::select(-dplyr::all_of(outcome_var)) |>
  as.data.frame()

y_train <- d_train_baked[[outcome_var]]

# ---- VSURF ----
vs <- VSURF::VSURF(
  y = y_train,
  x = x_train,
  ntree.thres   = 10000, nfor.thres  = 100,
  ntree.interp  = 10000, nfor.interp = 100,
  ntree.pred    = 10000, nfor.pred   = 100,
  RFimplem = "randomForest",
  parallel = TRUE
)
sel_idx <- vs$varselect.interp
vsurf_sel_vars <- if (length(sel_idx)) colnames(x_train)[sel_idx] else character(0)
#message(sprintf("OUTER %02d — VSURF selected (%d): %s",
#               i, length(vsurf_sel_vars), paste(vsurf_sel_vars, collapse=", ")))

#pred_terms = c("if_ng_26", "mip_1b_67", "tn_fa_76", "il_12_p40_43", "tn_fb_77")

ctrl_inner <- caret::trainControl(
  method = "LOOCV",
  classProbs = TRUE,
  savePredictions = "final",
  summaryFunction = caret::mnLogLoss,
  verboseIter = FALSE,
  allowParallel = FALSE
)

tg <- expand.grid(
  alpha  = seq(0, 1, by = 0.05),
  lambda = 10^seq(-4, 1, length.out = 100)
)

#---- Elastic-net selection on all predictors ----

elas_pre <- caret::train(
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

#Coef
coefs <- as.matrix(coef(elas_pre$finalModel,
                        s = elas_pre$bestTune$lambda[1]))[, 1]
el_sel_vars <- setdiff(names(coefs)[coefs != 0], "(Intercept)")
coef_elas_only <- coef(elas_pre$finalModel, s = elas_pre$bestTune$lambda)

#---- Final model with VSURF selected vars ----
rec_main <- make_recipe(el_sel_vars, outer_d_train, outcome_var)
#rec_main <- make_recipe(el_sel_vars, outer_d_train, outcome_var)

final_fit <- caret::train(
  rec_main,
  data = outer_d_train,
  method = "glmnet",
  metric = "logLoss",
  maximize = FALSE,
  tuneGrid = tg,
  trControl = ctrl_inner,
  family = "binomial",
  standardize = FALSE
)


#Coef
coef_final <- coef(final_fit$finalModel, s = final_fit$bestTune$lambda)
coefs <- as.matrix(coef(final_fit$finalModel,
                        s = final_fit$bestTune$lambda[1]))[, 1]
el_sel_vars <- setdiff(names(coefs)[coefs != 0], "(Intercept)")

# Predict a new case
new_case = df_imp_1 %>%
  slice(1) %>%
  mutate(if_ng_26 = 9, mip_1b_67 = 32, tn_fb_77 = 0.6, tn_fa_76 = 11, il_12_p40_43 = 0)

pred_probs <- predict(final_fit, newdata = new_case, type = "prob")[, positive_class] #Yes=0.96
pred_probs

#manual calculation of prob using coefs
df_imp_1_log <- df_imp_1 %>%
  select(where(is.numeric)) %>%
  mutate(across(everything(), log1p))

new_case_log = new_case %>%
  select(where(is.numeric)) %>%
  mutate(across(where(is.numeric), log1p)) %>%
  select(if_ng_26, mip_1b_67, tn_fa_76, il_12_p40_43, tn_fb_77)

new_case_log_center <- new_case_log - c(mean(df_imp_1_log$if_ng_26), 
                                        mean(df_imp_1_log$mip_1b_67),
                                        mean(df_imp_1_log$tn_fa_76),
                                        mean(df_imp_1_log$il_12_p40_43),
                                        mean(df_imp_1_log$tn_fb_77))

new_case_log_scale <- new_case_log_center / c(sd(df_imp_1_log$if_ng_26), 
                                              sd(df_imp_1_log$mip_1b_67),
                                              sd(df_imp_1_log$tn_fa_76),
                                              sd(df_imp_1_log$il_12_p40_43),
                                              sd(df_imp_1_log$tn_fb_77))

intercept = coef_final[1]
beta_vector = coef_final[-1]
log_odds = intercept + sum(beta_vector * as.numeric(new_case_log_scale))

predicted_prob <- 1 / (1 + exp(-log_odds)) #0.96 - matched the caret prediction
plogis(log_odds) #0.96 - matched the caret prediction

df <- data.frame(
  Yes = results$avg_final_outer_preds_single_elas,           
  No  = 1-results$avg_final_outer_preds_single_elas,      
  obs = results$outer_y_single_elas
)
evalm(df, positive = "Yes", optimise = "MCC") #similar - check


df_1 <- data.frame(
  Yes = results$avg_final_outer_preds_elas,           
  No  = 1-results$avg_final_outer_preds_elas,      
  obs = results$outer_y
)

evalm(df_1, positive = "Yes", optimise = "MCC") #similar - check

df_2 <- data.frame(
  Yes = results$avg_final_outer_preds_vs,           
  No  = 1-results$avg_final_outer_preds_vs,      
  obs = results$outer_y
)

evalm(df_2, positive = "Yes", optimise = "MCC") #similar - check
