# Nested LOOCV with VSURF for Variable Selection and Elastic Net with Calibration
# Author: AI Assistant
# Date: 2024

# Load required libraries
library(VSURF)
library(caret)
library(pROC)
library(glmnet)
library(ggplot2)
library(dplyr)
library(PRROC)

# Function to calculate PRAUC
calculate_prauc <- function(actual, predicted_probs) {
  if (length(unique(actual)) < 2) {
    return(0)
  }
  
  # Handle edge cases
  if (all(predicted_probs == 0)) {
    return(0)
  }
  
  tryCatch({
    pr_curve <- pr.curve(scores.class0 = predicted_probs[actual == 1],
                         scores.class1 = predicted_probs[actual == 0],
                         curve = TRUE)
    return(pr_curve$auc.integral)
  }, error = function(e) {
    return(0)
  })
}

# Function to perform VSURF variable selection
perform_vsurf_selection <- function(X_train, y_train) {
  # Ensure data types are correct
  X_train <- as.data.frame(X_train)
  y_train <- as.factor(y_train)
  
  # Perform VSURF
  vsurf_result <- VSURF(x = X_train, y = y_train, 
                        ntree = 100, 
                        mtry = max(1, floor(sqrt(ncol(X_train)))),
                        parallel = FALSE, 
                        ncores = 1)
  
  # Return interpretation set variables
  return(vsurf_result$varselect.interp)
}

# Function to train elastic net with caret
train_elastic_net <- function(X_train, y_train, X_val, y_val) {
  # Create training control for caret
  ctrl <- trainControl(
    method = "none",  # We'll handle CV manually
    classProbs = TRUE,
    summaryFunction = twoClassSummary
  )
  
  # Grid for alpha and lambda
  grid <- expand.grid(
    alpha = seq(0, 1, by = 0.1),
    lambda = 10^seq(-4, 1, length.out = 20)
  )
  
  best_prauc <- 0
  best_model <- NULL
  best_params <- NULL
  
  # Grid search for best parameters
  for (i in 1:nrow(grid)) {
    alpha_val <- grid$alpha[i]
    lambda_val <- grid$lambda[i]
    
    # Train model with current parameters
    model <- train(
      x = X_train, 
      y = y_train,
      method = "glmnet",
      trControl = ctrl,
      tuneGrid = data.frame(alpha = alpha_val, lambda = lambda_val),
      metric = "ROC",  # Caret doesn't have PRAUC built-in
      preProcess = c("center", "scale")
    )
    
    # Get predictions on validation set
    pred_probs <- predict(model, X_val, type = "prob")[, 2]
    
    # Calculate PRAUC
    prauc <- calculate_prauc(y_val, pred_probs)
    
    # Update best model if current PRAUC is better
    if (prauc > best_prauc) {
      best_prauc <- prauc
      best_model <- model
      best_params <- c(alpha = alpha_val, lambda = lambda_val)
    }
  }
  
  return(list(model = best_model, prauc = best_prauc, params = best_params))
}

# Function to assess calibration performance
assess_calibration <- function(actual, predicted_probs, n_bins = 10) {
  # Create bins
  bin_edges <- seq(0, 1, length.out = n_bins + 1)
  bin_centers <- (bin_edges[-1] + bin_edges[-(n_bins + 1)]) / 2
  
  # Assign predictions to bins
  bin_indices <- cut(predicted_probs, breaks = bin_edges, labels = FALSE, include.lowest = TRUE)
  
  # Calculate observed vs expected probabilities for each bin
  calibration_data <- data.frame(
    bin = 1:n_bins,
    bin_center = bin_centers,
    predicted_prob = bin_centers,
    observed_prob = numeric(n_bins),
    n_observations = numeric(n_bins)
  )
  
  for (i in 1:n_bins) {
    bin_mask <- bin_indices == i
    if (sum(bin_mask) > 0) {
      calibration_data$observed_prob[i] <- mean(actual[bin_mask] == 1)
      calibration_data$n_observations[i] <- sum(bin_mask)
    }
  }
  
  # Calculate calibration metrics
  # Hosmer-Lemeshow test statistic
  hl_statistic <- sum((calibration_data$observed_prob - calibration_data$predicted_prob)^2 * 
                      calibration_data$n_observations / 
                      (calibration_data$predicted_prob * (1 - calibration_data$predicted_prob)))
  
  # ECE (Expected Calibration Error)
  ece <- sum(abs(calibration_data$observed_prob - calibration_data$predicted_prob) * 
             calibration_data$n_observations) / sum(calibration_data$n_observations)
  
  return(list(
    calibration_data = calibration_data,
    hl_statistic = hl_statistic,
    ece = ece
  ))
}

# Main nested LOOCV function
nested_loocv_vsurf_elasticnet <- function(X, y, n_outer_folds = NULL) {
  n_samples <- nrow(X)
  
  # If n_outer_folds is NULL, use LOOCV
  if (is.null(n_outer_folds)) {
    n_outer_folds <- n_samples
  }
  
  # Initialize results storage
  outer_results <- list()
  all_predictions <- numeric(n_samples)
  all_actuals <- y
  
  # Outer loop
  for (outer_fold in 1:n_outer_folds) {
    if (n_outer_folds == n_samples) {
      # LOOCV: leave one sample out
      test_indices <- outer_fold
      train_indices <- setdiff(1:n_samples, test_indices)
    } else {
      # K-fold CV
      fold_size <- floor(n_samples / n_outer_folds)
      start_idx <- (outer_fold - 1) * fold_size + 1
      end_idx <- ifelse(outer_fold == n_outer_folds, n_samples, outer_fold * fold_size)
      test_indices <- start_idx:end_idx
      train_indices <- setdiff(1:n_samples, test_indices)
    }
    
    cat("Processing outer fold", outer_fold, "of", n_outer_folds, "\n")
    
    # Split data
    X_train_outer <- X[train_indices, , drop = FALSE]
    y_train_outer <- y[train_indices]
    X_test_outer <- X[test_indices, , drop = FALSE]
    y_test_outer <- y[test_indices]
    
    # Inner loop for VSURF variable selection
    cat("  Performing VSURF variable selection...\n")
    selected_vars <- perform_vsurf_selection(X_train_outer, y_train_outer)
    
    if (length(selected_vars) == 0) {
      cat("  Warning: No variables selected by VSURF, using all variables\n")
      selected_vars <- 1:ncol(X_train_outer)
    }
    
    cat("  Selected", length(selected_vars), "variables\n")
    
    # Use selected variables for training
    X_train_selected <- X_train_outer[, selected_vars, drop = FALSE]
    X_test_selected <- X_test_outer[, selected_vars, drop = FALSE]
    
    # Train elastic net with best parameters
    cat("  Training elastic net...\n")
    elastic_result <- train_elastic_net(X_train_selected, y_train_outer, 
                                       X_test_selected, y_test_outer)
    
    # Get predictions on test set
    test_predictions <- predict(elastic_result$model, X_test_selected, type = "prob")[, 2]
    
    # Store results
    outer_results[[outer_fold]] <- list(
      selected_variables = selected_vars,
      elastic_params = elastic_result$params,
      test_prauc = elastic_result$prauc,
      test_predictions = test_predictions,
      test_actuals = y_test_outer
    )
    
    # Store predictions for overall assessment
    all_predictions[test_indices] <- test_predictions
    
    cat("  Test PRAUC:", round(elastic_result$prauc, 4), "\n")
  }
  
  # Overall performance assessment
  overall_prauc <- calculate_prauc(all_actuals, all_predictions)
  
  # Calibration assessment
  calibration_result <- assess_calibration(all_actuals, all_predictions)
  
  # Compile final results
  final_results <- list(
    outer_results = outer_results,
    overall_predictions = all_predictions,
    overall_actuals = all_actuals,
    overall_prauc = overall_prauc,
    calibration = calibration_result,
    selected_variables_summary = table(unlist(lapply(outer_results, function(x) x$selected_variables)))
  )
  
  return(final_results)
}

# Function to plot calibration curve
plot_calibration <- function(calibration_result) {
  cal_data <- calibration_result$calibration_data
  
  p <- ggplot(cal_data, aes(x = predicted_prob, y = observed_prob)) +
    geom_point(size = 3, color = "blue") +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
    geom_line(color = "blue") +
    labs(title = "Calibration Plot",
         subtitle = paste("ECE =", round(calibration_result$ece, 4),
                         "| HL Statistic =", round(calibration_result$hl_statistic, 4)),
         x = "Predicted Probability",
         y = "Observed Probability") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5))
  
  return(p)
}

# Function to print summary results
print_summary <- function(results) {
  cat("\n" + "="*50 + "\n")
  cat("NESTED LOOCV WITH VSURF AND ELASTIC NET RESULTS\n")
  cat("="*50 + "\n")
  
  cat("\nOverall Performance:\n")
  cat("Overall PRAUC:", round(results$overall_prauc, 4), "\n")
  
  cat("\nCalibration Performance:\n")
  cat("Expected Calibration Error (ECE):", round(results$calibration$ece, 4), "\n")
  cat("Hosmer-Lemeshow Statistic:", round(results$calibration$hl_statistic, 4), "\n")
  
  cat("\nVariable Selection Summary:\n")
  var_summary <- results$selected_variables_summary
  for (i in 1:length(var_summary)) {
    cat("Variable", names(var_summary)[i], "selected", var_summary[i], "times\n")
  }
  
  cat("\nOuter Fold Results:\n")
  for (i in 1:length(results$outer_results)) {
    fold_result <- results$outer_results[[i]]
    cat("Fold", i, ": PRAUC =", round(fold_result$test_prauc, 4),
        "| Variables =", length(fold_result$selected_variables), "\n")
  }
}

# Example usage
# Uncomment and modify the following code to run the analysis:

# # Load your data
# # data <- read.csv("your_data.csv")
# # X <- data[, -ncol(data)]  # Features
# # y <- data[, ncol(data)]   # Target variable
# 
# # Ensure target is binary (0/1 or factor)
# y <- as.factor(y)
# 
# # Run nested LOOCV
# results <- nested_loocv_vsurf_elasticnet(X, y)
# 
# # Print summary
# print_summary(results)
# 
# # Plot calibration
# cal_plot <- plot_calibration(results$calibration)
# print(cal_plot)
# 
# # Save results
# saveRDS(results, "nested_loocv_results.rds")