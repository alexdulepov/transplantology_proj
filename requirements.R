# R Package Requirements for Nested LOOCV with VSURF and Elastic Net
# Install these packages before running the main script

# Core packages
install.packages("VSURF")
install.packages("caret")
install.packages("pROC")
install.packages("glmnet")
install.packages("ggplot2")
install.packages("dplyr")
install.packages("PRROC")

# Additional dependencies that might be needed
install.packages("randomForest")  # Required by VSURF
install.packages("e1071")        # Required by caret
install.packages("foreach")       # Required by caret
install.packages("iterators")     # Required by caret
install.packages("parallel")      # Required by caret

# If you prefer using BiocManager for some packages:
# if (!requireNamespace("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
# BiocManager::install("PRROC")

# Load all packages to check if they're working
library(VSURF)
library(caret)
library(pROC)
library(glmnet)
library(ggplot2)
library(dplyr)
library(PRROC)
library(randomForest)
library(e1071)

cat("All packages loaded successfully!\n")