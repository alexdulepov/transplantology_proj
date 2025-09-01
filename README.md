# Nested LOOCV with VSURF and Elastic Net

This repository contains R code for implementing nested Leave-One-Out Cross-Validation (LOOCV) with VSURF for variable selection and Elastic Net for classification, optimized for PRAUC (Precision-Recall Area Under Curve) and including calibration performance assessment.

## Features

- **Nested LOOCV**: Outer loop for model evaluation, inner loop for variable selection and parameters optimization
- **VSURF Variable Selection**: Uses VSURF (Variable Selection Using Random Forests) and elastic net for variable selection (comparison of the 2 methods)
- **Elastic Net Training**: Elastic net models trained using caret package
- **LogLoss Optimization**: Minimizes logloss (cross-entropy)
- **Calibration Assessment**: Comprehensive calibration performance metrics including calibration-in-the large, calibraion slope and calibraion plot
- **Decision analysis**: Implementation of the decision curve analysis for the assessment of the model clinical utility

## Files

- `nested_loocv_vsurf_elasticnet.R` - Main implementation with all functions
- `requirements.R` - Package installation requirements
- `demo_usage.R` - Example usage with sample data
- `README.md` - This documentation file

## Installation

1. Install required R packages:

```r
source("requirements.R")
```

Or install manually:

```r
install.packages(c("VSURF", "caret", "pROC", "glmnet", "ggplot2", "dplyr", "PRROC", "randomForest", "e1071"))
```

## Usage

### Basic Usage

```r
# Source the main functions
source("nested_loocv_vsurf_elasticnet.R")

# Load your data
data <- read.csv("your_data.csv")
X <- data[, -ncol(data)]  # Features
y <- data[, ncol(data)]   # Target variable

# Ensure target is binary (0/1 or factor)
y <- as.factor(y)

# Run nested LOOCV
results <- nested_loocv_vsurf_elasticnet(X, y)

# Print summary
print_summary(results)

# Plot calibration
cal_plot <- plot_calibration(results$calibration)
print(cal_plot)
```

### Demo with Sample Data

```r
# Run the demo to see the system in action
source("demo_usage.R")
```

## Function Details

### Main Function: `nested_loocv_vsurf_elasticnet()`

**Parameters:**
- `X`: Feature matrix (n_samples × n_features)
- `y`: Target variable vector (binary)
- `n_outer_folds`: Number of outer CV folds (NULL for LOOCV)

**Returns:**
- `outer_results`: Results for each outer fold
- `overall_predictions`: Combined predictions across all folds
- `overall_actuals`: Actual target values
- `overall_prauc`: Overall PRAUC performance
- `calibration`: Calibration performance metrics
- `selected_variables_summary`: Summary of variable selection frequency

### Variable Selection: `perform_vsurf_selection()`

Uses VSURF algorithm to select interpretable variables:
- Random forest-based variable importance
- Three-step selection process (thresholding, interpretation, prediction)
- Returns interpretation set variables

### Elastic Net Training: `train_elastic_net()`

Grid search over alpha and lambda parameters:
- Alpha: Mixing parameter (0 = Ridge, 1 = Lasso)
- Lambda: Regularization strength
- Optimizes for PRAUC on validation set

### Calibration Assessment: `assess_calibration()`

Provides comprehensive calibration metrics:
- **ECE (Expected Calibration Error)**: Average absolute difference between predicted and observed probabilities
- **Hosmer-Lemeshow Statistic**: Chi-squared test for calibration
- **Calibration Plot**: Visual assessment of probability calibration

## Output Interpretation

### Performance Metrics

- **PRAUC**: Precision-Recall AUC (higher is better, range 0-1)
- **Calibration ECE**: Expected Calibration Error (lower is better, range 0-1)
- **Hosmer-Lemeshow**: Chi-squared statistic for calibration (lower p-value indicates better calibration)

### Variable Selection

- **Selection Frequency**: How often each variable is selected across folds
- **Interpretation Set**: Variables selected by VSURF for interpretability
- **Stability**: Consistency of variable selection across folds

## Customization Options

### VSURF Parameters

```r
# Modify VSURF parameters in perform_vsurf_selection()
vsurf_result <- VSURF(x = X_train, y = y_train, 
                      ntree = 200,           # Number of trees
                      mtry = sqrt(ncol(X)),  # Variables per split
                      parallel = TRUE,       # Enable parallel processing
                      ncores = 4)            # Number of cores
```

### Elastic Net Grid

```r
# Modify parameter grid in train_elastic_net()
grid <- expand.grid(
  alpha = seq(0, 1, by = 0.05),           # Finer alpha grid
  lambda = 10^seq(-6, 2, length.out = 50) # Extended lambda range
)
```

### Calibration Bins

```r
# Modify number of calibration bins
calibration_result <- assess_calibration(actual, predicted_probs, n_bins = 20)
```

## Performance Considerations

- **LOOCV**: Computationally expensive for large datasets
- **VSURF**: Can be slow with many features or samples
- **Grid Search**: Parameter grid size affects training time
- **Memory**: Large datasets may require significant memory

### Optimization Tips

1. **Reduce VSURF trees**: Use fewer trees for faster execution
2. **Coarse parameter grid**: Start with coarse grid, refine later
3. **Parallel processing**: Enable VSURF parallel processing
4. **K-fold CV**: Use K-fold instead of LOOCV for large datasets

## Example Output

```
==================================================
NESTED LOOCV WITH VSURF AND ELASTIC NET RESULTS
==================================================

Overall Performance:
Overall PRAUC: 0.8234

Calibration Performance:
Expected Calibration Error (ECE): 0.0456
Hosmer-Lemeshow Statistic: 8.2341

Variable Selection Summary:
Variable 1 selected 95 times
Variable 2 selected 87 times
Variable 3 selected 92 times
...

Outer Fold Results:
Fold 1 : PRAUC = 0.8123 | Variables = 8
Fold 2 : PRAUC = 0.8345 | Variables = 7
...
```

## Troubleshooting

### Common Issues

1. **VSURF errors**: Ensure target variable is factor and features are numeric
2. **Memory issues**: Reduce VSURF trees or use K-fold CV
3. **Package conflicts**: Check package versions and dependencies
4. **Slow execution**: Enable parallel processing where possible

### Error Handling

The code includes error handling for:
- Edge cases in PRAUC calculation
- VSURF failures
- Invalid predictions
- Missing data

## Citation

If you use this code in your research, please cite:

- VSURF: Genuer, R., Poggi, J. M., & Tuleau-Malot, C. (2010). Variable selection using random forests. Pattern Recognition Letters, 31(14), 2225-2236.
- Caret: Kuhn, M. (2008). Building predictive models in R using the caret package. Journal of Statistical Software, 28(5), 1-26.

## License

This code is provided for educational and research purposes. Please ensure compliance with the licenses of the included packages.

## Support

For issues or questions:
1. Check the troubleshooting section
2. Verify package installations
3. Review error messages for specific issues
4. Ensure data format requirements are met
