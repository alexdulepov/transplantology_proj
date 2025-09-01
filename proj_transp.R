library(tidyverse)
library(readxl)
library(statip)
library(pheatmap)
library(VIM)
library(glmnet)
library(caret)
library(party)
library(flexplot)
library(VSURF)
library(janitor)
library(MLeval)
library(pROC)
library(tibble)
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
  mutate(outcome = factor(outcome, levels = c(0, 1)),
         sex = factor(sex_1_m_2_f, levels = c(1, 2), labels = c("Male", "Female")),
         time_tr = time_from_transplantation_month
         ) %>%
  select(patient_id, outcome, sex, age,time_tr, c(-1,-2,-5, -time_from_transplantation_month)) %>%
  mutate(outcome= factor(outcome, levels = c(1, 0), labels = c("Yes", "No")))

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

# ---- Setup ---------------------------------------------------------------
library(caret)
library(tibble)
library(dplyr)
library(ModelMetrics)
library(recipes)
library(glmnet)
library(yardstick)
library(ggplot2)   # for autoplot()
library(MLeval)    # optional: compare against yardstick
library(lattice)   # for caret::calibration plot
library(doParallel)
library(VSURF)
library(prg)
library(CalibrationCurves)
library(dcurves)


nested_elastic_binary_outcome <- function(
    df,
    outcome_var = "outcome",
    positive_class = "Yes",
    negative_class = "No",
    cv_outer_folds = 5,
    cv_outer_repeats = 20,
    seed = 1,
    alpha_grid  = seq(0, 1, by = 0.05),
    lambda_grid = 10^seq(-4, 1, length.out = 100)
) {
  # ---- Reproducible RNG ----
  RNGkind("L'Ecuyer-CMRG")
  set.seed(seed)
  
  # ---- Outcome factor with explicit order: negative first ----
  if (!all(c(negative_class, positive_class) %in% unique(df[[outcome_var]]))) {
    stop("Outcome does not contain both specified classes.")
  }
  df[[outcome_var]] <- factor(df[[outcome_var]],
                              levels = c(negative_class, positive_class))
  
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
      recipes::step_center(recipes::all_numeric_predictors()) |>
      recipes::step_scale(recipes::all_numeric_predictors()) |>
      recipes::step_dummy(recipes::all_nominal_predictors(), one_hot = FALSE)
    
    # After dummying, optionally keep only the selected post-dummy columns
    if (!is.null(vars_to_keep) && length(vars_to_keep) > 0) {
      rec <- rec |>
        recipes::step_select(recipes::all_outcomes(), tidyselect::any_of(vars_to_keep))
    }
    
    rec
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
    elastic_sel_vars   = vector("list", n_outer),
    vsurf_sel_vars     = vector("list", n_outer),
    final_vars_elastic = vector("list", n_outer)  # stores coefs for both methods
  )
  
  n_alg <- 2L
  inner_perf <- vector("list", n_outer)
  for (i in seq_len(n_outer)) {
    inner_perf[[i]] <- vector("list", n_alg)
  }
  
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
      ntree.thres   = 10000, nfor.thres  = 100,
      ntree.interp  = 10000, nfor.interp = 100,
      ntree.pred    = 10000, nfor.pred   = 100,
      RFimplem = "randomForest",
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
    
    list_vars_both <- list(vsurf_sel_vars, el_sel_vars)
    
    #---------------------- Inner CV with selected vars -------------------
    for (j in seq_along(list_vars_both)) {
      vars <- list_vars_both[[j]]
      rec_main <- make_recipe(vars, outer_d_train, outcome_var)
      
      if (length(vars) > 1) {
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
      } else {
        fit_inner <- caret::train(
          rec_main,
          data = outer_d_train,
          method = "glm",
          metric = "logLoss",
          maximize = FALSE,
          trControl = ctrl_inner,
          family = "binomial"
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
      sel_vars_df$final_vars_elastic[[i]][[j]] <- coefs2
      
      #---- Predictions on inner training (resamples) ----
      inner_preds <- fit_inner$pred
      p_in   <- inner_preds[[positive_class]]
      truth  <- inner_preds$obs
      
      inner_dat <- tibble::tibble(
        truth     = truth,
        truth_num = as.integer(truth == positive_class),
        .pred_Yes = p_in,
        hard_pred = factor(ifelse(p_in >= 0.5, positive_class, negative_class),
                           levels = levels(truth))
      )
      
      # ---- Metrics (yardstick with event = second; PRG for AUPRG) ----
      roc_auc       <- yardstick::roc_auc(inner_dat, truth, .pred_Yes, event_level = "second") |> dplyr::pull(.estimate)
      pr_auc        <- yardstick::pr_auc(inner_dat,  truth, .pred_Yes, event_level = "second") |> dplyr::pull(.estimate)
      prg_curve     <- prg::create_prg_curve(inner_dat$truth_num, inner_dat$.pred_Yes)
      auprg         <- prg::calc_auprg(prg_curve)
      logloss       <- yardstick::mn_log_loss(inner_dat, truth, .pred_Yes, event_level = "second") |> dplyr::pull(.estimate)
      if (fit_inner[["method"]] == "glmnet") {
        logloss_caret <- fit_inner$results %>%
          dplyr::filter(alpha == fit_inner$bestTune$alpha,
                        lambda == fit_inner$bestTune$lambda) %>%
          dplyr::pull(logLoss)
      } else {  # glm has no tuning grid; results has a single row
        logloss_caret <- fit_inner$results$logLoss[1]
      }
      brier = ModelMetrics::brier(actual = inner_dat$truth_num, predicted = inner_dat$.pred_Yes)
      mcc           <- yardstick::mcc(inner_dat,   truth, hard_pred,  event_level = "second") |> dplyr::pull(.estimate)
      
      # Predict on outer test
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
  }
  
  # ---- Aggregate outer predictions (average across repeats) -------------
  pred_df <- tibble::tibble()
  for (i in seq_along(inner_perf)) {
    for (j in seq_along(inner_perf[[i]])) {
      method <- if (j == 1) "VSURF" else "ElasticNet"
      pred_df <- dplyr::bind_rows(
        pred_df,
        tibble::tibble(
          idx     = inner_perf[[i]][[j]]$pred_idx,
          pred    = inner_perf[[i]][[j]]$preds,
          outcome = inner_perf[[i]][[j]]$y,
          method  = method
        )
      )
    }
  }
  
  avg_wide <- pred_df |>
    dplyr::group_by(idx, method) |>
    dplyr::summarise(avg_pred = mean(pred), .groups = "drop") |>
    tidyr::pivot_wider(names_from = method, values_from = avg_pred) |>
    dplyr::arrange(idx)
  
  outer_y <- pred_df |>
    dplyr::group_by(idx) |>
    dplyr::summarise(y = dplyr::first(as.character(outcome)), .groups = "drop") |>
    dplyr::arrange(idx) |>
    dplyr::pull(y)
  
  list(
    avg_final_outer_preds_vs   = avg_wide$VSURF,
    avg_final_outer_preds_elas = avg_wide$ElasticNet,
    outer_y = outer_y,
    sel_vars_df = sel_vars_df,
    inner_perf = inner_perf
  )
}


results = nested_elastic_binary_outcome(
  df_imp_1, 
  outcome_var = "outcome",
  positive_class = "Yes",
  negative_class = "No",
  cv_outer_folds = 5,
  cv_outer_repeats = 20,
  seed = 1,
  alpha_grid = seq(0, 1, by = 0.05),
  lambda_grid = 10^seq(-4, 1, length.out = 100))



stab_table <- function(sel_list) {
  vec <- unlist(sel_list, use.names = FALSE)
  if (!length(vec)) return(tibble(variable = character(0), times_selected = integer(0), freq = numeric(0)))
  dfc <- as.data.frame(table(vec), stringsAsFactors = FALSE)
  names(dfc) <- c("variable", "times_selected")
  dfc <- dfc[order(dfc$times_selected, decreasing = TRUE), , drop = FALSE]
  dfc$freq <- dfc$times_selected / N
  tibble::as_tibble(dfc)
}

stability <- list(
  elasticPreElast = stab_table(res_elasticPreElast$sel_list),
  vsurf           = stab_table(res_vsurf$sel_list),
  elasticOnly     = tibble(variable = character(0), times_selected = integer(0), freq = numeric(0)),
  elasticPreElast_glm = stab_table(res_elasticPreElast_glm$sel_list),
  vsurf_glm           = stab_table(res_vsurf_glm$sel_list)
)

#VSURF metrics

rocauc = yardstick::roc_auc_vec(
  truth = factor(results$outer_y, levels = c("No", "Yes")),
  estimate = results$avg_final_outer_preds_vs,
  event_level = "second"
)

prauc = yardstick::pr_auc_vec(
  truth = factor(results$outer_y, levels = c("No", "Yes")),
  estimate = results$avg_final_outer_preds_vs,
  event_level = "second"
)

logloss = yardstick::mn_log_loss_vec(
  truth = factor(results$outer_y, levels = c("No", "Yes")),
  estimate = results$avg_final_outer_preds_vs,
  event_level = "second"
)

brier = ModelMetrics::brier(actual = as.integer(results$outer_y == "Yes"), predicted = results$avg_final_outer_preds_vs)

mcc = yardstick::mcc_vec(
  truth = factor(results$outer_y, levels = c("No", "Yes")),
  estimate = factor(ifelse(results$avg_final_outer_preds_vs >= 0.5, "Yes", "No"),
                    levels = c("No", "Yes")),
  event_level = "second"
)




#PRG GAIN AUC VSURF
prg_curve     <- prg::create_prg_curve(as.integer(results$outer_y == "Yes"), results$avg_final_outer_preds_vs)
auprg         <- prg::calc_auprg(prg_curve)
auprg

#Calibration slope and intercept
y01 <- as.integer(results$outer_y == "Yes")
logit_p <- qlogis(results$avg_final_outer_preds_vs)
fit_cal  <- stats::glm(y01 ~ logit_p, family = binomial())
cal_int  <- unname(coef(fit_cal)[1])
cal_slope<- unname(coef(fit_cal)[2])
fit_citl <- stats::glm(y01 ~ 1 + offset(logit_p), family = binomial())
citl     <- unname(coef(fit_citl)[1])
cal_int
cal_slope
citl
rocauc
prauc
logloss
brier
mcc

#Elastic Net metrics

rocauc_elas = yardstick::roc_auc_vec(
  truth = factor(results$outer_y, levels = c("No", "Yes")),
  estimate = results$avg_final_outer_preds_elas,
  event_level = "second"
)

prauc_elas = yardstick::pr_auc_vec(
  truth = factor(results$outer_y, levels = c("No", "Yes")),
  estimate = results$avg_final_outer_preds_elas,
  event_level = "second"
)

logloss_elas = yardstick::mn_log_loss_vec(
  truth = factor(results$outer_y, levels = c("No", "Yes")),
  estimate = results$avg_final_outer_preds_elas,
  event_level = "second"
)


ModelMetrics::prAUC(actual = results$outer_y, predicted = results$avg_final_outer_preds_elas)
MLmetrics::PRAUC(y_pred = results$avg_final_outer_preds_elas, y_true =  as.integer(results$outer_y == "Yes"))   # y in {0,1}

brier_elas = yardstick::brier_class(
  df,
  obs,
  Yes,
  event_level = "first"
)

brier = ModelMetrics::brier(actual = as.integer(results$outer_y == "Yes"), predicted = results$avg_final_outer_preds_elas)

mcc_elas = yardstick::mcc_vec(
  truth = factor(results$outer_y, levels = c("No", "Yes")),
  estimate = factor(ifelse(results$avg_final_outer_preds_vs >= 0.5, "Yes", "No"),
                    levels = c("No", "Yes")),
  event_level = "second"
)

#PRGGAIN AUC
prg_curve_elas     <- prg::create_prg_curve(as.integer(results$outer_y == "Yes"), results$avg_final_outer_preds_elas)
auprg_elas         <- prg::calc_auprg(prg_curve_elas)
auprg_elas

#Calibration slope and intercept

y01_elas <- as.integer(results$outer_y == "Yes")
logit_p_elas <- qlogis(results$avg_final_outer_preds_elas)
fit_cal_elas  <- stats::glm(y01_elas ~ logit_p_elas, family = binomial())
cal_int_elas  <- unname(coef(fit_cal_elas)[1])
cal_slope_elas<- unname(coef(fit_cal_elas)[2])
fit_citl_elas <- stats::glm(y01_elas ~ 1 + offset(logit_p_elas), family = binomial())
citl_elas     <- unname(coef(fit_citl_elas)[1])
cal_int_elas
cal_slope_elas
citl_elas
rocauc_elas
prauc_elas
logloss_elas
brier_elas
mcc_elas

# Calibration plot for VSURF
cal_plot_vs <- val.prob.ci.2(
  results$avg_final_outer_preds_vs,
  y01_elas,
  logistic = TRUE
) 

# Calibration plot for Elastic Net

cal_plot_elas <- val.prob.ci.2(
  results$avg_final_outer_preds_elas,
  y01_elas,
  logistic = TRUE,
  smooth = "loess"
)


val.prob(results$avg_final_outer_preds_elas, y01_elas)

# Diagnostic plot for VSURF
#df for dca
res_df = tibble(
  pred_elas = results$avg_final_outer_preds_elas,
  pred_vs = results$avg_final_outer_preds_vs,
  y01_elas = y01_elas,
  y = results$outer_y
)

df <- data.frame(
  Yes = results$avg_final_outer_preds_elas,           
  No  = 1-results$avg_final_outer_preds_elas,      
  obs = results$outer_y
)
brier = brier_score(df, positive = colnames(df)[1]) #Different RESULT!!!!!!!!!!!!!!!!!!!!!!!!!! in comparison to yardstick CHECK

df <- data.frame(
  Yes = results$avg_final_outer_preds_vs,           
  No  = 1-results$avg_final_outer_preds_vs,      
  obs = results$outer_y
)
eval = evalm(df) #similar - check

brier = brier_score(df, positive = colnames(df)[1]) #ANOTHER RESULT!!!!!!!!!!!!!!!!!!!!!!!!!! in comparison to yardstick CHECK
log_like = LL(df, positive = colnames(df)[1]) #similar - check

diag_plot_vs <- dcurves::dca(y01_elas ~ pred_elas + pred_vs, data = res_df, thresholds = seq(0, 0.5, by = 0.01),) 
diag_plot_vs

dca(y01_elas ~ pred_elas + pred_vs,
    data = res_df,
    thresholds = seq(0, 0.5, by = 0.01),
    label = list(pred_elas = "Variable selection with elastic net", pred_vs = "Variable selection with VSURF")
) %>%
  plot(smooth = TRUE)


#loocv prauc noncal ROC 0.98, PRG 0.71, PRC 0.79(MLeval+Mlmetrics), MCC 0.77 (std+yardstick) MCC 0.81 (threshold 0.98) not good calibr plot brier 0.061, LL -3.81 (MLeval) brier 0.06 (MLeval+yardstick) logloss 0.17 (yardstick+MLmetrics
#loocv mcc noncal ROC 0.96, PRG 0.65, PRC 0.72(MLeval+Mlmetrics), MCC 0.64 (std) MCC 0.81 (threshold 0.38) good calibr plot brier 0.078, LL -4.92 (MLeval) brier 0.77 logloss 0.22 (yardstick+MLmetrics

################################################################################Final model#######################################################################

# Final model with best hyperparameters
# Inner CV: optimize PR AUC (caret::prSummary returns PR AUC as "AUC")

pred_terms = c("if_ng_26", "mip_1b_67", "tn_fa_76", "il_12_p40_43", "tn_fb_77")

ctrl_inner <- trainControl(
  method = "LOOCV",
  classProbs = TRUE,
  savePredictions = "final",
  summaryFunction = mnLogLoss,   # PR metrics (AUC here = PR AUC)
  verboseIter = FALSE,
  allowParallel = FALSE
)

tg <- expand.grid(
  alpha  = seq(0, 1, by = 0.1),
  lambda = 10^seq(-4, 1, length.out = 50)
)

rec <- recipe(reformulate(pred_terms, response = "outcome"), data = df_imp_1) |>
  step_log(all_numeric_predictors(), offset = 1) |>
  step_center(all_predictors()) |>
  step_scale(all_predictors()) |>
  step_impute_knn(all_predictors(), neighbors = 5)

final_fit <- train(
  rec,
  data = df_imp_1,
  method = "glmnet",
  metric = "logLoss",            
  maximize = FALSE,
  tuneGrid = tg,
  trControl = ctrl_inner,
  standardize = FALSE,        # recipe handles it
  family = "binomial"
)

#Coef
coef_final <- coef(final_fit$finalModel, s = final_fit$bestTune$lambda)

# Predict on the full dataset

new_case = data.frame(if_ng_26 = 9, mip_1b_67 = 32, tn_fa_76 = 0, il_12_p40_43 = 0, tn_fb_77 = 0.6)
pred_probs <- predict(final_fit, newdata = new_case, type = "prob")[, positive_class] #Yes=0.78 No=0.22



#manual calculation of prob using coefs
df_imp_1_log <- df_imp_1 %>%
  select(all_of(pred_terms)) %>%
  mutate(across(everything(), log1p))


new_case_log = log1p(new_case) # log1p transformation

new_case_log_center <- new_case_log - c(mean(df_imp_1_log$if_ng_26), mean(df_imp_1_log$mip_1b_67),
                                        mean(df_imp_1_log$tn_fa_76),mean(df_imp_1_log$il_12_p40_43), mean(df_imp_1_log$tn_fb_77))

new_case_log_scale <- new_case_log_center / c(sd(df_imp_1_log$if_ng_26), sd(df_imp_1_log$mip_1b_67),
                                              sd(df_imp_1_log$tn_fa_76), sd(df_imp_1_log$il_12_p40_43), sd(df_imp_1_log$tn_fb_77))

log_odds <- sum(coef_final[-1] * as.numeric(new_case_log_scale)) + coef_final[1]

predicted_prob <- 1 / (1 + exp(-log_odds)) #0.86 - matched the caret prediction
plogis(log_odds) #0.86 - matched the caret prediction

#brier - mean error in prediction. Max 1 and min is 0. The case happened(1) and the model predicted 0.86, so the error is (1-0.86)^2 = 0.02 or 0.14 error in prob on average. Sum all the errors for all predictions and divide by the number of predictions = the mean error(loss)
#brier 0.27*(1-0.27)^2 + (1-0.27)*0.27^2= 0.1971 no better than  non-informative model. Depends on the prevalence. If the prevalence is 50%, the ineff model has 0.25 brier score


#Decision curves




#################More complicated code

# assumes you already trained `final_fit` exactly as before
pred_terms <- c("if_ng_26","mip_1b_67","tn_fa_76","il_12_p40_43","tn_fb_77")

# 1) Training means/SDs after log1p (publish these)
mu  <- sapply(pred_terms, function(v) mean(log1p(df_imp_1[[v]])))
sdv <- sapply(pred_terms, function(v)   sd(log1p(df_imp_1[[v]])))

# 2) Coefficients at the selected lambda (publish these)
cf      <- coef(final_fit$finalModel, s = final_fit$bestTune$lambda)
beta0   <- as.numeric(cf[1, , drop = TRUE])     # intercept
betas   <- as.numeric(cf[-1, , drop = TRUE])
names(betas) <- rownames(cf)[-1]

# (optional) ensure the same order across all vectors
ord <- names(betas)
mu_pub  <- mu[ord]; sd_pub <- sdv[ord]; beta_pub <- betas[ord]

# 3) Table to give to clinicians / include in paper
replication_table <- tibble::tibble(
  variable = ord,
  mean_log1p = as.numeric(mu_pub),
  sd_log1p   = as.numeric(sd_pub),
  beta       = as.numeric(beta_pub)
)

intercept <- beta0
replication_table
intercept

# 4) Sanity check: manual vs caret on an example case
new_case <- data.frame(
  if_ng_26 = 9, mip_1b_67 = 32, tn_fa_76 = 0,
  il_12_p40_43 = 0, tn_fb_77 = 0.6
)

# caret probability of Yes
p_yes_caret <- predict(final_fit, newdata = new_case, type = "prob")[, "Yes"]

# manual probability of Yes (using the publishable formula)
xlog <- log1p(unlist(new_case[ord]))
z    <- (xlog - mu_pub) / sd_pub
eta  <- intercept + sum(beta_pub * z)
p_yes_manual <- 1 / (1 + exp(-eta))  # logistic(-eta)

c(caret_Yes = p_yes_caret, manual_Yes = p_yes_manual, diff = p_yes_caret - p_yes_manual)
