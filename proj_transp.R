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
set.seed(1234)

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
  select(-c(egf_13, il_17e_il_25_51, il_17f_53, il_22_55))#%>%
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
  select(3:length(df_imp_1))

y_train <- df_imp_1 %>% 
  pull(outcome) 

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

############################################################################################## Elastic Net Regression with caret

# y: factor with the *positive class first* (needed for ROC in caret)
# Example: y <- factor(y, levels = c("case","control"))
ctrl <- trainControl(
  method = "repeatedcv",
  number = 5,
  repeats = 100,
  classProbs = TRUE,
  summaryFunction = twoClassSummary,   # gives ROC, Sens, Spec
  savePredictions = "final",
  returnResamp = "all",
  allowParallel = T
  # sampling = "down"  # or "up"/"smote" if you have class imbalance
)

ctrl <- trainControl(
  method          = "LOOCV",
  classProbs      = TRUE,
  summaryFunction = twoClassSummary,
  savePredictions = "final",
  returnResamp = "all")


ctrl <- trainControl(
  method          = "boot632",
  number = 1000,
  classProbs      = TRUE,
  summaryFunction = twoClassSummary,
  savePredictions = "final",
  returnResamp = "all"
)


fit_1 <- train(
  outcome ~ if_ng_26,
  method = "glm",
  family = "binomial",  
  preProcess = c("expoTrans","center","scale"),
  metric = "ROC",
  trControl = ctrl,
  data = df_imp_1
)

fit_2 <- train(
  outcome ~ if_ng_26 + mip_1b_67,
  method = "glmnet",
  preProcess = c("expoTrans","center","scale"),  
  metric = "ROC",
  tuneGrid = expand.grid(
    alpha  = seq(0, 1, by = 0.1),
    lambda = 10^seq(-4, 1, length.out = 50)
  ),
  trControl = ctrl,
  standardize = FALSE,
  data = df_imp_1
)

fit_3 <- train(
  outcome ~ if_ng_26 + mip_1b_67+tn_fa_76,
  method = "glmnet",
  preProcess = c("expoTrans","center","scale"), 
  metric = "ROC",
  tuneGrid = expand.grid(
    alpha  = seq(0, 1, by = 0.1),
    lambda = 10^seq(-4, 1, length.out = 50)
  ),
  trControl = ctrl,
  standardize = FALSE,
  data = df_imp_1
)

fit_4 <- train(
  outcome ~ if_ng_26 + mip_1b_67+tn_fa_76 + il_12_p40_43,
  method = "glmnet",
  preProcess = c("expoTrans","center","scale"), 
  metric = "ROC",
  tuneGrid = expand.grid(
    alpha  = seq(0, 1, by = 0.1),
    lambda = 10^seq(-4, 1, length.out = 50)
  ),
  trControl = ctrl,
  standardize = FALSE,
  data = df_imp_1
)

fit_5 <- train(
  outcome ~ if_ng_26 + mip_1b_67+tn_fa_76 + il_12_p40_43 + tn_fb_77,
  method = "glmnet",
  preProcess = c("expoTrans","center","scale"),  
  metric = "ROC",
  tuneGrid = expand.grid(
    alpha  = seq(0, 1, by = 0.1),
    lambda = 10^seq(-4, 1, length.out = 50)
  ),
  trControl = ctrl,
  standardize = FALSE,
  data = df_imp_1
)

# Coef from regula glm with 1 var
summary(fit_1$finalModel)$coefficients  

#Coef and vars from elastic 
model_names <- sprintf("fit_%d", 1:5)
fits <- mget(model_names, inherits = TRUE)  

# Helper to extract best alpha/lambda and coefficients at that lambda
extract_glmnet <- function(fit, name) {
  best_alpha  <- fit$bestTune$alpha[1]
  best_lambda <- fit$bestTune$lambda[1]
  
  # Coefficients at best lambda (includes "(Intercept)")
  cm <- coef(fit$finalModel, s = best_lambda)
  v  <- as.numeric(cm)
  names(v) <- rownames(cm)
  
  # Keep nonzero coefficients (optional)
  v <- v[v != 0]
  
  tibble(
    fit         = name,
    best_alpha  = best_alpha,
    best_lambda = best_lambda,
    coef        = list(v)#,     # list-column
    #auc        = fit$results$ROC[fit$bestTune$alpha == best_alpha & fit$bestTune$lambda == best_lambda],
    #auc_sd = fit$results$ROCSD[fit$bestTune$alpha == best_alpha & fit$bestTune$lambda == best_lambda]
  )
}

elastic_df <- bind_rows(
  lapply(names(fits), function(nm) extract_glmnet(fits[[nm]], nm))
)

elastic_df$coef #tn_fa and il excluded from fit_5 => unstable + tn_fb has multiple null values

######################################Performance evaluation

res <- evalm(list(fit_1, fit_2, fit_3, fit_4, fit_5), gnames = c("if","if+mip","if+mip+tn_fa","if+mip+tn_fa+il","if+mip+tn_fa+il+tn_fb"), 
             rlinethick = 0.8, fsize = 8, plots = "cc", positive = "Yes", optimise = "MCC", bins=4)

res <- evalm(list(fit_4),  
             rlinethick = 0.8, fsize = 8, plots = "cc", positive = "Yes", optimise = "MCC", bins=4)
#calibration curve is much better for if+mip

res$optres 
res$stdres 
res[["probs"]]$"if+mip"[which.max(res[["probs"]]$"if+mip"$MCC),]

roc_obj <- roc(response = res[["probs"]][["if+mip"]][["obs"]], predictor =res[["probs"]][["if+mip"]][["Yes"]], quiet = TRUE)
opt_cut <- coords(roc_obj, "best", ret = c("threshold","sensitivity","specificity"), best.method="youden")

#fit_2
#ggplot(fit_2)


#####################################Calibration
# Intercept and slope

cal_df = data.frame(
  fit = c(),
  intercept = c(),
  slope = c() 
)

for (i in names(res$optres)) {
  p =  res[["probs"]][[i]][["Yes"]]
  y = res[["probs"]][[i]][["obs"]]
  glm_cal <- glm(y ~ qlogis(p), family = binomial)
  inter = glm_cal$coefficients[1]
  slope = glm_cal$coefficients[2]
  name = i
  cal_df = rbind(cal_df, data.frame(fit = name, intercept = inter, slope = slope))
}
  
cal_df

#Brier score
brier_ll__df = data.frame(
  fit = c(),
  brier_score = c(),
  log_likelihood = c()
)

for (i in names(res$optres)) {
  df <- data.frame(
    Yes = res[["probs"]][[i]][["Yes"]],           
    No  = res[["probs"]][[i]][["No"]],      
    obs = res[["probs"]][[i]][["obs"]]
  )
  brier = brier_score(df, positive = colnames(df)[1])
  log_like = LL(df, positive = colnames(df)[1])
  name = i
  brier_ll__df = rbind(brier_ll__df, data.frame(fit = name, brier_score = brier, log_likelihood = log_like))
}

brier_ll__df

# caret calibration plot
preds <- res[["probs"]][[1]]
cal_raw <- caret::calibration(obs ~ Yes, data = preds, class = "Yes", cuts = 4)
xyplot(cal_raw, auto.key = list(columns = 1))

############## NESTED LOOCV???????
#Model assumptions?
