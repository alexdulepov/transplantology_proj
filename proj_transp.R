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
#############################################Algorithm run

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
inner_res <- inner_model_perf_nested_binary(trained_object = results)

# ---- Outer resample summaries ----
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
    thresholds = seq(0.05, 0.4, 0.01)
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

fin_mpd = final_model_with_coefs(df_imp_1,
                                   outcome_var = "outcome",
                                   positive_class = "Yes",
                                   negative_class = "No",
                                   family = "gaussian",
                                   cv_method = "repeatedcv",
                                   cv_folds = 5,
                                   cv_repeats = 20,
                                   alpha_grid  = seq(0, 1, by = 0.1),
                                   lambda_grid = 10^seq(-4, 1, length.out = 20),
                                   ntree = 1000,
                                   nfor = 20,
                                   selection_rule = "oneSE")
##############################################################################################Manual calculation of probs using coefs
new_case = df_imp_1 %>%
  slice(1) %>%
  mutate(if_ng_26 = 9, mip_1b_67 = 32, tn_fb_77 = 0.6, tn_fa_76 = 11, il_12_p40_43 = 0)

pred_probs <- predict(final_fit, newdata = new_case, type = "prob")[, positive_class] #Yes=0.96
pred_probs

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

