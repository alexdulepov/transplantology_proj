library(tidyverse)
library(readxl)
library(statip)
library(pheatmap)
library(VIM)
library(writexl)
library(glmnet)
library(caret)
library(party)
library(flexplot)
library(VSURF)
library(janitor)
library(car)
library(tidymodels)
library(olsrr)
library(pROC)
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
  select(-c(egf_13, il_17e_il_25_51, il_17f_53, il_22_55))#%>%
  #slice(-c(6))  
#dont delete observation 2 because this is the only one female among cases and you will get the separation issue,
#however this observation is an outlier(low value of if while being rejected), so you can keep it for sensitivity analysis

x_train = df_imp_1 %>%
  select(3:length(df_imp)) %>%
  mutate(
    across(                                                    
      .cols = where(is.numeric),                  
      .fns  = ~ log1p(.x),                                 
      .names = "{.col}"                                  
    ),
    across(                                                    
      .cols = where(is.numeric),                  
      .fns  = ~ scale(.x),                                 
      .names = "{.col}"                                  
    )) 

y_train <- df_imp_1 %>% 
  pull(outcome) %>% 
  factor(levels = c(0, 1),
         labels = c("No", "Yes"))

x1 = df_imp_1 %>%
  select(age,time_tr,if_ng_26, il_12_p40_43, tn_fb_77, mip_1b_67, tn_fa_76) %>%
  cor()

x2 = x_train %>%
  select(age,time_tr,if_ng_26, il_12_p40_43, tn_fb_77, mip_1b_67, tn_fa_76) %>%
  cor()
  
boxplot(df_imp,outline=FALSE,col="cornflowerblue")
##################################################Correlation analysis

cor_mat=cor(df_imp[,4:length(df_imp)], method = "pearson")

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

############ ELASTIC NET REGRESSION
set.seed(42)

x <- x_train  # design matrix
y <- as.numeric(y_train)-1  # response vector

# Cross‑validate over both λ and α
alphas <- seq(0, 1, length = 100)                  # 0, 0.1, …, 1
cv_list <- lapply(alphas, function(a)
  cv.glmnet(x, y, alpha = a, nfolds = 22, standardize = F, type.measure = "class"))

# Pick the alpha with the lowest CV error
cv_means <- sapply(cv_list, function(z) min(z$cvm))
best_a   <- alphas[ which.min(cv_means) ]
best_fit <- cv_list[[ which.min(cv_means) ]]

coef(best_fit, s = "lambda.min")  # coefficients at the best λ


#mip_1b_67, if_ng_26, tn_fb_77           


###
    # only needed for ggplot variant, optional

## ---------------- 1.  coefficient paths ----------------
## Re‑fit one glmnet model at the chosen alpha so we get
## every λ along the path (cv.glmnet keeps it too, but this
## makes the call explicit and lets us plot all variables).

fit_path <- glmnet(
  x, y,
  alpha       = best_a,          # <‑‑ winner from your CV loop
  family      = "binomial",
  standardize = TRUE
)

# Base‑graphics quick plot (built‑in to glmnet):
plot(fit_path, xvar = "lambda", label = TRUE)   # label names on right margin
abline(v = log(best_fit$lambda.min), lty = 2, col = "red")  # CV‑chosen λ
title(main = sprintf("Elastic‑net coefficient paths (alpha = %.2f)", best_a))

set.seed(1234)

###################################VSURF

vsurf_model = VSURF(y=y_train, x=x_train, ntree.thres = 10000,nfor.thres = 100, 
                    ntree.interp = 10000, nfor.interp=100, 
                    ntree.pred = 10000, nfor.pred = 100,
                    RFimplem = "ranger")
summary(vsurf_model)
variables=attr(vsurf_model$terms, "term.labels") #the variables order used in VSURF

plot(vsurf_model, step="thres", imp.mean=FALSE, main="Variable importance plot") # variable importance plot
#green line - a smoothing function using CART
#red line - threshold which is a minimum predicted value from green line (everything below threshold is rejected for further steps)
variables[vsurf_model$varselect.thres] #the selected variables after threshold step]
colnames(x_train[vsurf_model[["varselect.thres"]]])

plot(vsurf_model, step="interp", imp.mean=FALSE, main="Variable importance plot")
#red line - threshold, the smallest model with OOB less than worst OOB+ 1 sd
variables[vsurf_model$varselect.interp]
colnames(x_train[vsurf_model[["varselect.interp"]]])

plot(vsurf_model, step="pred", imp.mean=FALSE, main="Variable importance plot")
#red line - threshold, the smallest model with OOB less than worst OOB+ 1 sd
variables[vsurf_model$varselect.pred] #if_ng_26 + il_12_p40_43 +(tn_fb_77?)
colnames(x_train[vsurf_model[["varselect.pred"]]])

mod_rfor = cforest(outcome ~ if_ng_26+il_12_p40_43+mip_1b_67+tn_fa_76, data = df_imp[,-1])
compare.fits(outcome ~ if_ng_26, data = df_imp[,-1], mod_rfor)

data("avengers")
mod_rfor = cforest(ptsd ~ ., data = avengers)
compare.fits(ptsd ~ strength , data =avengers, mod_rfor)


glm_mod = glm(outcome~if_ng_26+il_12_p40_43, data = df_imp_1[,c(2:length(df_imp_1))],family = "binomial")
summary(glm_mod)
exp(coef(glm_mod))
vif(glm_mod) # Variance inflation factor, no multicollinearity

ggplot(df_imp, aes(x=if_ng_26, y=outcome)) +
  geom_point()

ggplot(df_imp, aes(x=if_ng_26)) +
  geom_histogram()

stats::cooks.distance(glm_mod) 
# Plot Cook's distance
ols_plot_dfbetas(glm_mod)
ols_plot_cooksd_chart(glm_mod, type = 1, threshold = NULL, print_plot = TRUE)

#6, 22 obs are influential - keep track of outliers. better to check all the variables for outliers. sens anaylsis with and without outliers.
#plot logistic curve


############Confusion matrix
# Get the actual responses from churn
actual_response <- df_imp_1$outcome

# Get the predicted responses from the model
predicted_response <- round(fitted(glm_mod))

# Get a table of these values
outcomes <- table(predicted_response, actual_response)

# Convert the table to a conf_mat object
confusion <- conf_mat(outcomes)

# "Automatically" plot the confusion matrix
autoplot(confusion)

# Get summary metrics
summary(confusion, event_level = "second")


####LOOCV
library(caret)

df_imp_0$outcome <- factor(df_imp_0$outcome,
                           levels = c(0, 1),
                           labels = c("No", "Yes"))

ctrl <- trainControl(
  method          = "LOOCV",
  savePredictions = "final",   # <-- keeps every hold‑out prediction
  classProbs      = TRUE,
  summaryFunction = twoClassSummary 
)

trctrl <- trainControl(method = "boot632",number=10000,
                       returnResamp="all", savePredictions = "final",   # <-- keeps every hold‑out prediction
                       classProbs      = TRUE,
                       summaryFunction = twoClassSummary)

set.seed(42)
loocv_fit <- train(
  y= y_train, 
  x= x_train[,c("if_ng_26","il_12_p40_43"), drop = FALSE],
  method      = "glm",
  family      = binomial,
  trControl   = trctrl,
  metric     = "ROC"
)

#sens 0.833 with 2 var
#1 var - more stable, but less precise
#2 var - les stable , but more precise

#bootstrap with 2 var : sens = 0.71, spec =0.897, ROC = 0.798, PPV=0.72,NPV=0.89, bal.Ac = 0.81 (10000 samples enough)
#bootstrap with 1 var : sens = 0.67, spec =0.91, ROC = 0.790, PPV=0.737,NPV=0.88, bal.Ac = 0.79

loocv_fit$results
outcomes <- table(loocv_fit$pred$pred, loocv_fit$pred$obs)
confusionMatrix(outcomes, positive = "Yes")

# Convert the table to a conf_mat object
confusion <- conf_mat(outcomes)

# "Automatically" plot the confusion matrix
autoplot(confusion)

# Get summary metrics
summary(confusion, event_level = "second")

rocCurve <- pROC::roc(
  response  = loocv_fit$pred$obs,   # factor: "No"/"Yes"
  predictor = loocv_fit$pred$Yes    # numeric: P(class = "Yes")
)
plot(rocCurve, print.auc = TRUE)
plot(x1, print.auc = TRUE) #1 var
plot(x2, print.auc = TRUE) #2 var
pROC::coords(rocCurve, "best", ret = c("threshold", "sensitivity", "specificity"))

###
# choose the cut-off threshold based on the ROC curve

# 1) Build the ROC (make sure the positive class is correct)
rocCurve <- roc(
  response  = loocv_fit$pred$obs,     # factor with levels c("No","Yes")
  predictor = loocv_fit$pred$Yes,     # P(class == "Yes")
  levels = c("No","Yes"),             # ensure "Yes" is the positive class
  direction = "<"                     # higher prob -> "Yes"
)

# 2) Get the best threshold (Youden’s J or closest-to-(0,1))
best <- coords(
  rocCurve,
  x = "best",
  best.method = "youden",             # or "closest.topleft"
  ret = c("threshold","sensitivity","specificity"),
  transpose = FALSE                   # returns a named vector
)

best_threshold <- best["threshold"]
best_threshold #0.22
thr90sens <- pROC::coords(rocCurve, x = 0.80, input = "sensitivity", ret = "threshold")
thr95spec <- pROC::coords(rocCurve, x = 0.85, input = "specificity", ret = "threshold")

### pr-recall
library(precrec)

scores <- loocv_fit$pred$Yes
labels <- loocv_fit$pred$obs  # factor with "No"/"Yes"

mmod <- evalmod(scores = scores,
                labels = labels,
                posclass = "Yes")   # <-- not 'positive'
# includes AUPRC and AUROC
autoplot(mmod, "PRC") + ggtitle("Precision–Recall")
autoplot(mmod, "ROC") + ggtitle("ROC")
ap_df <- auc(mmod)                        # AUCs for both ROC and PR
auprc <- subset(ap_df, curvetypes == "PRC")$aucs
auprc

# using your objects
scores <- loocv_fit$pred$Yes
labels <- factor(loocv_fit$pred$obs, levels = c("No","Yes"))

prevalence <- mean(labels == "Yes")   # PR baseline
auprc <- 0.61                          # your PR AUC
roc_auc <- 0.79                        # your ROC AUC

# "Lift" over baseline, scaled to [0,1]
pr_lift <- (auprc - prevalence) / (1 - prevalence)
prevalence; pr_lift

#Prevalence ≈ 0.271 → a random classifier’s baseline AUPRC ≈ 0.271.

#Your AUPRC = 0.61 → +0.339 absolute above baseline, or about 2.25× the baseline (0.61 / 0.271).

#pr_lift = 0.465 → you’re 46.5% of the way from random to perfect in PR space (where 0 = baseline, 1 = perfect).

#ROC AUC = 0.79 looks higher because ROC is less sensitive to class imbalance; with ~27% positives, it’s normal for AUPRC to be notably lower than ROC AUC.

library(dplyr); library(yardstick); library(ggplot2)

df <- loocv_fit$pred %>%
  transmute(truth = factor(obs, levels = c("No","Yes")),
            score = Yes)

pr <- pr_curve(df, truth, score, event_level = "second")

best <- pr %>%
  mutate(F1 = 2*precision*recall/(precision+recall)) %>%
  filter(!is.na(F1)) %>% slice_max(F1, n = 1)

best$.threshold; best$precision; best$recall #best threshold is 0.42 based on F1 score

#Plot the PR curve + annotate baseline and the chosen point
#The dashed horizontal line is the prevalence (= baseline precision of a random classifier).

#The point marks your best‑F1 operating point (recall, precision) on the curve.

autoplot(pr) +
  geom_hline(yintercept = mean(df$truth == "Yes"), linetype = 2) +
  geom_point(aes(x = best$recall, y = best$precision))

#Uncertainty
library(precrec)

set.seed(1)
B <- 2000
n <- nrow(df)
ap <- roc <- numeric(B)

for (b in 1:B) {
  i <- sample.int(n, n, replace = TRUE)
  m <- evalmod(scores = df$score[i], labels = df$truth[i], posclass = "Yes")
  ap[b]  <- subset(auc(m),  curvetypes == "PRC")$aucs
  roc[b] <- subset(auc(m),  curvetypes == "ROC")$aucs
}

quantile(ap,  c(.025, .5, .975))   # AUPRC CI
quantile(roc, c(.025, .5, .975))   # ROC AUC CI



#roc for if_ng_26 + il_12_p40_43 = 0.781 , but sens is higher, however 2 times were convergence issue
#roc for if_ng_26 = 0.823
# The further from 1 that the CPR is, the more misleading ROC analysis can be. In the cases where CPR is not close to 1, consider 
#using precision-recall analysis or evaluating Matthew’s correlation coefficient (MCC) (e.g. plotting MCC versus sensitivity), 
#especially when the algorithm’s positive predictions will be used for clinical decision making.+#ROC + IMBALANce robust metrics

################
library(caret)
library(glmnet)

## 1 ───────────────────────────────────────────────────────────────
##  Control object: leave‑one‑out CV + ROC
ctrl <- trainControl(
  method          = "LOOCV",
  savePredictions = "final",
  classProbs      = TRUE,
  summaryFunction = twoClassSummary          # gives ROC, Sens, Spec
)

## 2 ───────────────────────────────────────────────────────────────
##  Tuning grid for Elastic‑Net
grid <- expand.grid(
  alpha  = seq(0, 1,  by = 0.10),            # 0 = ridge, 1 = lasso
  lambda = 10^seq(-4, 1, length = 25)        # shrinkage values
)

## 3 ───────────────────────────────────────────────────────────────
##  Predictor matrix: *all* 60 columns
#  (x_train is your pre‑processed data‑frame/tibble:  n × 60)
x_mat <- as.matrix(x_train)                  # keeps column names

## 4 ───────────────────────────────────────────────────────────────
##  LOOCV Elastic‑Net, optimised on ROC AUC
set.seed(42)
loocv_fit_en <- train(
  x          = x_mat,
  y          = y_train,                      # factor with two levels
  method     = "glmnet",
  family     = "binomial",
  tuneGrid   = grid,                         # or tuneLength = 25
  metric     = "Spec",                        # must match summaryFunction
  trControl  = ctrl
)

## 5 ───────────────────────────────────────────────────────────────
##  Quick checks
loocv_fit_en$bestTune                  # α & λ chosen on LOOCV ROC
getTrainPerf(loocv_fit_en)            # ROC, Sens, Spec of best model
head(loocv_fit_en$results)            # full resampling table

## 6 ───────────────────────────────────────────────────────────────
##  Confusion matrix for the best α‑λ pair
library(dplyr)
best <- loocv_fit_en$bestTune

pred_best <- loocv_fit_en$pred %>% 
  filter(alpha  == best$alpha,
         lambda == best$lambda)

confusionMatrix(pred_best$pred,
                pred_best$obs,
                positive = "Yes")      # adjust if your “positive” ≠ "Yes"

# Pull the best lambda chosen by caret
best_lambda <- loocv_fit_en$bestTune$lambda

# Get the coefficient vector at that λ
beta <- coef(loocv_fit_en$finalModel, s = best_lambda)

# Drop the intercept and keep non‑zero terms
keep_idx <- which(beta != 0)[-1]        # -1 removes the intercept row
selected_vars <- rownames(beta)[keep_idx]

selected_vars
