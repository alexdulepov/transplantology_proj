# Demo Usage of Nested LOOCV with VSURF and Elastic Net
# This script demonstrates how to use the main functions with sample data

# Source the main functions
source("Functions_nested_CV(VSURF+elastic_net).R")

# Generate sample data for demonstration
set.seed(123)
n_samples <- 100
n_features <- 20

# Generate features (some informative, some noise)
X <- matrix(rnorm(n_samples * n_features), nrow = n_samples, ncol = n_features)

# Create informative features (first 5 features are important)
informative_features <- 1:5
X[, informative_features] <- X[, informative_features] + rnorm(n_samples, 0, 0.5)

# Generate target variable based on informative features
logit_prob <- 0.5 + 0.8 * X[, 1] - 0.6 * X[, 2] + 0.4 * X[, 3] - 0.3 * X[, 4] + 0.2 * X[, 5]
prob <- 1 / (1 + exp(-logit_prob))
y <- rbinom(n_samples, 1, prob)

# Convert to factor for classification
y <- as.factor(y)

# Check class balance
cat("Class distribution:\n")
print(table(y))

# Run nested LOOCV (this will take some time)
cat("\nStarting nested LOOCV with VSURF and Elastic Net...\n")
cat("This may take several minutes depending on your data size...\n\n")

# For demonstration, you can also use fewer folds to speed up the process
# results <- nested_loocv_vsurf_elasticnet(X, y, n_outer_folds = 10)

# Full LOOCV
results <- nested_elastic_binary_outcome(X, y)

# Print comprehensive summary
print_summary(results)

# Plot calibration curve
cal_plot <- plot_calibration(results$calibration)
print(cal_plot)

# Additional analysis and visualization

# 1. Variable selection frequency analysis
cat("\nVariable Selection Analysis:\n")
var_freq <- results$selected_variables_summary
var_freq_df <- data.frame(
  variable = as.numeric(names(var_freq)),
  frequency = as.numeric(var_freq),
  is_informative = as.numeric(names(var_freq)) %in% informative_features
)

print(var_freq_df)

# 2. Performance across folds
fold_performance <- data.frame(
  fold = 1:length(results$outer_results),
  prauc = sapply(results$outer_results, function(x) x$test_prauc),
  n_variables = sapply(results$outer_results, function(x) length(x$selected_variables))
)

cat("\nPerformance Summary Across Folds:\n")
print(summary(fold_performance$prauc))

# 3. Create performance plots
library(ggplot2)

# PRAUC distribution across folds
prauc_plot <- ggplot(fold_performance, aes(x = prauc)) +
  geom_histogram(bins = 10, fill = "steelblue", alpha = 0.7) +
  geom_vline(xintercept = results$overall_prauc, color = "red", linetype = "dashed") +
  labs(title = "Distribution of PRAUC Across Folds",
       subtitle = paste("Overall PRAUC:", round(results$overall_prauc, 4)),
       x = "PRAUC", y = "Frequency") +
  theme_minimal()

print(prauc_plot)

# Number of variables selected across folds
var_count_plot <- ggplot(fold_performance, aes(x = n_variables)) +
  geom_histogram(bins = 10, fill = "darkgreen", alpha = 0.7) +
  labs(title = "Distribution of Number of Variables Selected Across Folds",
       x = "Number of Variables", y = "Frequency") +
  theme_minimal()

print(var_count_plot)

# 4. Save results
saveRDS(results, "demo_nested_loocv_results.rds")
cat("\nResults saved to 'demo_nested_loocv_results.rds'\n")

# 5. Example of how to load and use saved results
cat("\nExample of loading saved results:\n")
loaded_results <- readRDS("demo_nested_loocv_results.rds")
cat("Loaded results contain", length(loaded_results$outer_results), "folds\n")

# 6. Extract specific information from results
cat("\nExtracting specific information:\n")
cat("Best performing fold:", which.max(sapply(loaded_results$outer_results, function(x) x$test_prauc)), "\n")
cat("Worst performing fold:", which.min(sapply(loaded_results$outer_results, function(x) x$test_prauc)), "\n")

# 7. Create a summary report
cat("\n" + "="*60 + "\n")
cat("DEMO SUMMARY REPORT\n")
cat("="*60 + "\n")
cat("Dataset: Generated sample data with", n_samples, "samples and", n_features, "features\n")
cat("Informative features:", paste(informative_features, collapse = ", "), "\n")
cat("Cross-validation: Nested LOOCV with VSURF variable selection\n")
cat("Model: Elastic Net with PRAUC optimization\n")
cat("Overall PRAUC:", round(results$overall_prauc, 4), "\n")
cat("Calibration ECE:", round(results$calibration$ece, 4), "\n")
cat("Average variables selected:", round(mean(fold_performance$n_variables), 2), "\n")
cat("="*60 + "\n")