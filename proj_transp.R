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
library(janitor)
library(report)
set.seed(123)


# Read xlsx and recode the outcome
df <- read_xlsx("D:/Packages/proj_transpl/Cytokines_heart_data_24.06 1 and 4 groups.xlsx",
                sheet = 1)

clean_df <- df %>%
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

###########################Missingness checks only
missingness_check <- clean_df %>%
  summarise(across(everything(), ~ sum(is.na(.)))) %>%
  pivot_longer(everything(), names_to = "variable", values_to = "n_missing") %>%
  mutate(percent_missing = n_missing / nrow(clean_df) * 100) %>%
  filter(n_missing > 0) %>%
  arrange(desc(n_missing))

print(missingness_check)

row_missingness_check <- clean_df %>%
  transmute(
    patient_id,
    n_missing = rowSums(is.na(across(everything()))),
    percent_missing = n_missing / ncol(clean_df) * 100
  ) %>%
  filter(n_missing > 0) %>%
  arrange(desc(n_missing))

print(row_missingness_check)

##################################################Near-zero-variance checks only
low_variance_check <- clean_df %>%
  select(-patient_id, -outcome) %>%
  caret::nearZeroVar(saveMetrics = TRUE) %>%
  tibble::rownames_to_column("variable") %>%
  filter(zeroVar | nzv | percentUnique < 50) %>%
  arrange(percentUnique, desc(freqRatio))

print(low_variance_check)

# No imputation or low-variance filtering is applied here. The modeling
# functions handle imputation inside resampling recipes; patient_id is removed
# only because it is an identifier, not a biological predictor.
model_df <- clean_df %>%
  select(-patient_id)

#dont delete observation 2 because this is the only one female among cases and you will get the separation issue,
#however this observation is an outlier(low value of if while being rejected), so you can keep it for sensitivity analysis

boxplot(model_df %>% select(where(is.numeric)),outline=FALSE,col="cornflowerblue")
##################################################Correlation analysis

numeric_predictors <- model_df %>%
  select(where(is.numeric))

cor_mat=cor(numeric_predictors, method = "pearson", use = "pairwise.complete.obs")

#' Find correlated variable pairs in a correlation matrix

find_correlated <- function(cor_mat, cutoff = 0.70) {

  if (!is.matrix(cor_mat) || !is.numeric(cor_mat))
    stop("cor_mat must be a numeric matrix (e.g. the output of cor()).")

  if (!identical(rownames(cor_mat), colnames(cor_mat)))
    stop("Row and column names must match (a square correlation matrix).")

  ## Keep only the upper triangle so we do not report duplicates or self-correlations.
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

vsurf_exploration_df <- model_df %>%
  mutate(across(where(is.numeric), ~ log1p(.))) %>%
  tidyr::drop_na()

x_train = vsurf_exploration_df %>%
  select(-outcome)

y_train <- vsurf_exploration_df %>%
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

ggplot(model_df, aes(x=if_ng_26)) +
  geom_histogram()

ggplot(model_df, aes(x=il_12_p40_43)) +
  geom_histogram()

ggplot(model_df, aes(x=tn_fb_77)) +
  geom_histogram() #different scale

ggplot(model_df, aes(x=mip_1b_67)) +
  geom_histogram()

ggplot(model_df, aes(x=tn_fa_76)) +
  geom_histogram()

selected_marker_cor <- model_df %>%
  select(if_ng_26, il_12_p40_43, tn_fb_77, mip_1b_67, tn_fa_76) %>%
  mutate(across(where(is.numeric), ~ log1p(.))) %>%
  cor(use = "pairwise.complete.obs")  #transformation is desirable

selected_marker_cor


rf = cforest(outcome ~ if_ng_26 , data = model_df)
estimates(rf)
#############################################Algorithm run

results = nested_elastic_binary_outcome(
  model_df,
  outcome_var = "outcome",
  positive_class = "Yes",
  negative_class = "No",
  cv_outer_folds = 5,
  cv_outer_repeats = 20,
  seed = 1,
  alpha_grid  = seq(0, 1, by = 0.1),
  lambda_grid = 10^seq(-4, 1, length.out = 50),
  hyperparam_search = "grid",
  random_search_size = 25,
  ntree = 1000,
  nforests = 20,
  inner_cv_method = "LOOCV",
  inner_cv_repeats = 20,
  inner_cv_folds = 5,
  selection_rule = "best",
  numeric_imputation = "median",
  nominal_imputation = "mode",
  add_missing_indicators = "no",
  remove_zero_variance = "no",
  add_interactions = "no",
  add_polynomials = "no")

# ---- Outer resample summaries ----
outer_res <- outer_perf_nested_binary(trained_object = results,
                                      positive_class = "Yes",
                                      negative_class = "No")

outer_predictions <- results$outer_predictions

#AUC ROC-the probability is that a randomly chosen cancer patient is ranked higher (given higher prob) than a randomly chosen healthy patient
#the AUC judges a correctly assigned biopsy (TP) and an unnecessary biopsy (FP) as equally important.
# ---- Calibration plots (rms) ----
y01 <- as.integer(outer_predictions$outcome == "Yes")  # or pass in positive_class if not "Yes"

cal_plot_vsurf_inter <- val.prob.ci.2(
  outer_predictions$pred_vsurf_interpretation,
  y01,
  logistic = TRUE,
  col.log = "blue",
  CL.smooth = FALSE
)

cal_plot_vsurf_pred <- val.prob.ci.2(
  outer_predictions$pred_vsurf_prediction,
  y01,
  logistic = TRUE,
  col.log = "blue",
  CL.smooth = FALSE
)

cal_plot_elas <- val.prob.ci.2(
  outer_predictions$pred_elastic_net,
  y01,
  logistic = TRUE,
  col.log = "blue",
  CL.smooth = FALSE
)

cal_plot_lasso <- val.prob.ci.2(
  outer_predictions$pred_lasso,
  y01,
  logistic = TRUE,
  col.log = "blue",
  CL.smooth = FALSE
)

cal_plot_prev <- val.prob.ci.2(
  outer_predictions$pred_prevalence,
  y01,
  logistic = TRUE,
  col.log = "blue",
  CL.smooth = FALSE
)

# Optional extra check
val.prob(outer_predictions$pred_elastic_net,          y01)
val.prob(outer_predictions$pred_lasso,                y01)
val.prob(outer_predictions$pred_vsurf_interpretation, y01)
val.prob(outer_predictions$pred_vsurf_prediction,     y01)
val.prob(outer_predictions$pred_prevalence,           y01)


# ---- Decision curve analysis (dcurves) ----
res_df <- tibble::tibble(
  pred_elas = outer_predictions$pred_elastic_net,
  pred_lasso = outer_predictions$pred_lasso,
  pred_vsurf_inter = outer_predictions$pred_vsurf_interpretation,
  pred_vsurf_pred = outer_predictions$pred_vsurf_prediction,
  pred_baseline = outer_predictions$pred_prevalence,
  y01       = y01,
  if_ng_26 = model_df$if_ng_26[outer_predictions$row_id],
  sex = model_df$sex[outer_predictions$row_id]
)

dcurves::dca(
  y01 ~ pred_elas + pred_lasso + pred_vsurf_inter + pred_vsurf_pred + pred_baseline,
  data = res_df,
  thresholds = seq(0, 0.6, by = 0.01)
) %>%
  plot(smooth = T)

#The value of 0.16 at a threshold probability of 20% can be interpreted as: compared with no intervention,
#intervention on the basis of the elastic model is the equivalent of a strategy that found 16 rejections per hundred patients without conducting any unnecessary interventions.

dca(y01 ~ pred_elas + pred_lasso + pred_vsurf_inter + pred_vsurf_pred + pred_baseline,
    data = res_df,
    thresholds = seq(0.05, 0.4, 0.01)
) %>%
  net_intervention_avoided() %>%
  plot(smooth = T)

#At a probability threshold of 15-30%, the net reduction in interventions is about 0.25.
#In other words, at this probability threshold, biopsying patients on the basis of the model is the equivalent of
#a strategy that led to an absolute 25% reduction in the number of biopsies without missing any heart rejections (25 unnecessary biopsies avoided per 100 patients).

dca(y01 ~ pred_elas + pred_lasso + pred_vsurf_inter + pred_vsurf_pred + pred_baseline,
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

fin_mpd = final_model_with_coefs(model_df,
                                   outcome_var = "outcome",
                                   positive_class = "Yes",
                                   negative_class = "No",
                                   family = "binomial",
                                   cv_method = "LOOCV",
                                   cv_folds = 5,
                                   cv_repeats = 20,
                                   alpha_grid  = seq(0, 1, by = 0.1),
                                   lambda_grid = 10^seq(-4, 1, length.out = 100),
                                   hyperparam_search = "grid",
                                   random_search_size = 25,
                                   ntree = 1000,
                                   nforests = 20,
                                   selection_rule = "best",
                                   numeric_imputation = "median",
                                   nominal_imputation = "mode",
                                   add_missing_indicators = "no",
                                   remove_zero_variance = "no",
                                   add_interactions = "no",
                                   add_polynomials = "no",
                                   post_selection_refit = "vsurf_only")
##############################################################################################Example prediction using the VSURF post-selection model
new_case = model_df %>%
  slice(1) %>%
  mutate(if_ng_26 = 9, mip_1b_67 = 32, tn_fb_77 = 0.6, tn_fa_76 = 11, il_12_p40_43 = 0)

pred_probs <- predict(fin_mpd$post_selection_models$VSURF, newdata = new_case, type = "prob")[, "Yes"] #Yes=0.96
pred_probs

model_df_log <- model_df %>%
  select(where(is.numeric)) %>%
  mutate(across(everything(), log1p))

new_case_log = new_case %>%
  select(where(is.numeric)) %>%
  mutate(across(where(is.numeric), log1p)) %>%
  select(if_ng_26, mip_1b_67, tn_fa_76, il_12_p40_43, tn_fb_77)

new_case_log_center <- new_case_log - c(mean(model_df_log$if_ng_26, na.rm = TRUE),
                                        mean(model_df_log$mip_1b_67, na.rm = TRUE),
                                        mean(model_df_log$tn_fa_76, na.rm = TRUE),
                                        mean(model_df_log$il_12_p40_43, na.rm = TRUE),
                                        mean(model_df_log$tn_fb_77, na.rm = TRUE))

new_case_log_scale <- new_case_log_center / c(sd(model_df_log$if_ng_26, na.rm = TRUE),
                                              sd(model_df_log$mip_1b_67, na.rm = TRUE),
                                              sd(model_df_log$tn_fa_76, na.rm = TRUE),
                                              sd(model_df_log$il_12_p40_43, na.rm = TRUE),
                                              sd(model_df_log$tn_fb_77, na.rm = TRUE))

vsurf_post_coefs <- fin_mpd$coefficient_estimates$post_selection[[1]]$VSURF
intercept = vsurf_post_coefs["(Intercept)"]
beta_vector = vsurf_post_coefs[setdiff(names(vsurf_post_coefs), "(Intercept)")]
matching_terms <- intersect(names(beta_vector), names(new_case_log_scale))
log_odds = intercept + sum(beta_vector[matching_terms] * as.numeric(new_case_log_scale[matching_terms]))

predicted_prob <- 1 / (1 + exp(-log_odds)) #0.96 - matched the caret prediction
plogis(log_odds) #0.96 - matched the caret prediction

eval_elastic <- data.frame(
  Yes = outer_predictions$pred_elastic_net,
  No  = 1 - outer_predictions$pred_elastic_net,
  obs = outer_predictions$outcome
)
evalm(eval_elastic, positive = "Yes", optimise = "MCC") #similar - check

df_lasso <- data.frame(
  Yes = outer_predictions$pred_lasso,
  No  = 1 - outer_predictions$pred_lasso,
  obs = outer_predictions$outcome
)

evalm(df_lasso, positive = "Yes", optimise = "MCC") #similar - check


eval_vsurf_inter <- data.frame(
  Yes = outer_predictions$pred_vsurf_interpretation,
  No  = 1 - outer_predictions$pred_vsurf_interpretation,
  obs = outer_predictions$outcome
)

evalm(eval_vsurf_inter, positive = "Yes", optimise = "MCC") #similar - check

eval_vsurf_pred <- data.frame(
  Yes = outer_predictions$pred_vsurf_prediction,
  No  = 1 - outer_predictions$pred_vsurf_prediction,
  obs = outer_predictions$outcome
)

evalm(eval_vsurf_pred, positive = "Yes", optimise = "MCC") #similar - check
