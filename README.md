# resilience
Forest Resilience Analysis

## Overview

This project contains Python scripts for time series analysis and machine learning model training, including 1-lag autocorrelation (AR1) analysis, STL (Seasonal-Trend decomposition using LOESS) decomposition, XGBoost model training with hyperparameter tuning, and SHAP (SHapley Additive exPlanations) model interpretation.

## Project Structure
├── XGBoost_SHAP_ADs.py # XGBoost training with SHAP feature importance analysis
├── AR1_Var_calculating.py # AR1 autoregressive analysis and variance calculation
├── STL_Decomposition.py # STL decomposition for time series
├── shap_values.csv # SHAP values output
├── shap_importance.csv # SHAP feature importance scores
├── data/ # Input data directory
└── output/ # Output results directory


## Dependencies

- Python >= 3.7
- pandas >= 1.3.0
- numpy >= 1.21.0
- scipy >= 1.7.0
- statsmodels >= 0.13.0
- xgboost >= 1.5.0
- scikit-learn >= 0.24.0
- shap >= 0.39.0

## Install dependencies using:
```bash
pip install pandas numpy scipy statsmodels xgboost scikit-learn shap


##Script Descriptions

1. XGBoost_SHAP_ADs.py

This script trains an XGBoost regression model using RandomizedSearchCV for hyperparameter optimization, evaluates model performance, and generates SHAP explanations.

Key Features:

Automated hyperparameter tuning using RandomizedSearchCV

Model performance evaluation (RMSE, R²)

SHAP value calculation for model interpretation

Feature importance ranking based on SHAP values

Usage:

python XGBoost_SHAP_ADs.py

Input:

CSV file containing features (X1, X2, X3, ...) and target variable (Y)

Output:

shap_values.csv: SHAP values for each sample

shap_importance.csv: Feature importance scores

SHAP summary plot

Key Parameters:

n_iter: Number of parameter combinations to try (default: 50)

test_size: Proportion of test set (default: 0.2)

random_state: Random seed for reproducibility (default: 4016)


2. AR1_Var_calculating.py

This script calculates AR1 (lag-1 autocorrelation) coefficients and variance using a sliding window approach for ecosystem resilience assessment.

Key Features:

Sliding temporal window analysis for AR1 coefficient calculation

Variance computation within each window

Support for multiple window sizes

Parallel processing for efficiency

Usage:

python AR1_Var_calculating.py


Input:

Time series residual data from STL decomposition

Output:

AR1 coefficient time series

Variance time series

Key Parameters:

window_sizes: List of sliding window sizes (e.g., [36, 48, 60, 72, 84])

n_processes: Number of parallel processes (default: 50)

date_start, date_end: Analysis period

3. STL_Decomposition.py

This script performs Seasonal-Trend decomposition using LOESS on vegetation index time series data.

Key Features:

Robust STL decomposition

Automatic calculation of trend and low-pass filter lengths

Parallel processing for multiple time series

Separation of trend, seasonal, and residual components

Usage:

python STL_Decomposition.py

Input:

Time series data in CSV format

Output:

trend/: Trend component files

season/: Seasonal component files

resid/: Residual component files

Key Parameters:

period: Period of seasonal component (default: 12)

smooth_length: Length of seasonal smoother (default: 7)

n_processes: Number of parallel processes (default: 50)

## Workflow

1. Data Preprocessing

Prepare time series data in CSV format

Ensure proper date indexing

2. STL Decomposition

python STL_Decomposition.py

-Separates time series into trend, seasonal, and residual components

3. AR1 and Variance Calculation

python AR1_Var_calculating.py

-Calculates resilience indicators from residuals

4. XGBoost Model Training

python XGBoost_SHAP_ADs.py

-Trains model with optimal hyperparameters

-Generates SHAP interpretations

## Configuration

Before running the scripts, update the following paths in each file:

# Example configuration

input_path = "./path/to/input/data/"

output_path = "./path/to/output/"

## Contact

Author: Yiling Cai

Institution: School of Geography and Remote Sensing, Guangzhou University, Guangzhou 510006, China;
            Guangdong Guodi Science Technology Co., Ltd, Guangzhou 510075, China
			
Email: caiyling6@alumni.sysu.edu.cn
