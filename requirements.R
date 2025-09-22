# R Package Requirements for Nested LOOCV with VSURF and Elastic Net
# Install these packages before running the main script

# Core + extras
install.packages(c(
  "VSURF","caret","glmnet","tidyverse","recipes","ModelMetrics","yardstick",
  "MLeval","CalibrationCurves","dcurves","readxl","statip","pheatmap","janitor",
  "doParallel","devtools","ranger","kknn","foreach","iterators","VIM"
))

library(devtools)
devtools::install_github("meeliskull/prg", subdir = "R_package/prg")

# Load
library(caret)
library(tidyverse)
library(ModelMetrics)
library(recipes)
library(glmnet)
library(yardstick)
library(MLeval)
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


cat("All packages loaded successfully!\n")