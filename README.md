# Overview

This project contains Python scripts for time series analysis and machine learning model training, including 1-lag autocorrelation (AR1) analysis, STL (Seasonal-Trend decomposition using LOESS) decomposition, XGBoost model training with hyperparameter tuning, and SHAP (SHapley Additive exPlanations) model interpretation.

# Project Structure
├── STL_Decomposition.py # STL decomposition for time series  
├── AR1_Var_calculating.py # AR1 autoregressive analysis and variance calculation  
├── XGBoost_SHAP_ADs.py # XGBoost training with SHAP feature importance analysis  
├── input/ # Input data directory  
└── output/ # Output results directory  


# Requirements

- Python >= 3.7
- pandas >= 1.3.0
- numpy >= 1.21.0
- scipy >= 1.7.0
- statsmodels >= 0.13.0
- xgboost >= 1.5.0
- scikit-learn >= 0.24.0
- shap >= 0.39.0

# Install dependencies using:
pip install pandas numpy scipy statsmodels xgboost scikit-learn shap  

# Script Description  

## 1. STL_Decomposition.py

This script performs Seasonal-Trend decomposition using LOESS on vegetation index time series data.

### Key Features:

Robust STL decomposition

Automatic calculation of trend and low-pass filter lengths

Parallel processing for multiple time series

Separation of trend, seasonal, and residual components

### Usage:

python STL_Decomposition.py

### Input:

Time series data in CSV format

### Output:

trend/: Trend component files

season/: Seasonal component files

resid/: Residual component files

### Key Parameters:

period: Period of seasonal component (default: 12)

smooth_length: Length of seasonal smoother (default: 7)

n_processes: Number of parallel processes (default: 50)

date_start, date_end: Analysis period (default: 1982, 2022)

vegetation_index: vegetation index (default: kNDVI)

## 2. AR1_Var_calculating.py

This script calculates AR1 (lag-1 autocorrelation) coefficients and variance using a sliding window approach for ecosystem resilience assessment.

### Key Features:

Sliding temporal window analysis for AR1 coefficient calculationC

Variance computation within each window

Support for multiple window sizes

Parallel processing for efficiency

### Usage:

python AR1_Var_calculating.py

### Input:

Time series residual data from STL decomposition

### Output:

AR1 coefficient time series

Variance time series

### Key Parameters:

window_sizes: List of sliding window sizes (e.g., [36, 48, 60, 72, 84])

n_processes: Number of parallel processes (default: 50)

date_start, date_end: Analysis period (default: 1982, 2022)

vegetation_index: vegetation index (default: kNDVI)

## 3. XGBoost_SHAP_ADs.py

This script trains an XGBoost regression model using RandomizedSearchCV for hyperparameter optimization, evaluates model performance, and generates SHAP explanations.

### Key Features:

Automated hyperparameter tuning using RandomizedSearchCV

Model performance evaluation (RMSE, R²)

SHAP value calculation for model interpretation

Feature importance ranking based on SHAP values

### Usage:

python XGBoost_SHAP_ADs.py

### Input:

CSV file containing features (X1, X2, X3, ...) and target variable (Y)

### Output:

shap_values.csv: SHAP values for each sample

shap_importance.csv: Feature importance scores

SHAP summary plot

### Key Parameters:

n_iter: Number of parameter combinations to try (default: 50)

test_size: Proportion of test set (default: 0.2)

random_state: Random seed for reproducibility (default: 4016)

# Workflow

## 1. Data Preprocessing

Prepare time series data in CSV format

## 2. STL Decomposition

python STL_Decomposition.py

-Separates time series into trend, seasonal, and residual components

## 3. AR1 and Variance Calculation

python AR1_Var_calculating.py

-Calculates resilience indicators from residuals

## 4. XGBoost Model Training

python XGBoost_SHAP_ADs.py

-Trains model with optimal hyperparameters

-Generates SHAP interpretations

# Configuration

Before running the scripts, update the following paths in each file:

# Example configuration

input_path = "./path/to/input/"

output_path = "./path/to/output/"

# Data Format

## Input CSV Requirements:

First column: series_id (example: PID = xxxxxxx)

Subsequent columns: features (X1, X2, X3, ...) and target variable (Y), or time series data (format:YYYY-MM-DD)

No missing values in critical analysis periods

### Example1:

PID	1982-01-01	1982-02-01	1982-03-01	1982-04-01 ...

3	0.0034 	0.0053 	0.0132 	0.0585 

1	0.0023 	0.0106 	0.0110 	0.0851 

4	0.0023 	0.0108 	0.0112 	0.1008 

5	0.0032 	0.0106 	0.0119 	0.1098 

2	0.1118 	0.1929 	0.0193 	0.0816 

...

### Example2:

PID	TAC_AD	TAC_mean	FC	FL	FG ...

516142	0.3546 	0.1514 	57.8624 	0.0288 	2.0814 

387407	0.2421 	0.0808 	87.4868 	0.0306 	0.4815 

370881	0.2758 	0.1918 	78.5291 	0.0627 	1.0123 

406671	0.0851 	0.1240 	75.0423 	0.0378 	1.6677 

354255	0.1316 	0.0987 	61.6825 	13.3503 	0.1026 

...


## Output Files

### STL Decomposition:

trend/{vegetation_index}_STL_trend_{series_id}.csv: Trend component

season/{vegetation_index}_STL_season_{series_id}.csv: Seasonal component

resid/{vegetation_index}_STL_resid_{series_id}.csv: Residual component

### AR1 Analysis:

{vegetation_index}_AR1_{window_size}/{vegetation_index}_AR1_{window_size}_STL_{series_id}.csv: AR1 coefficients

{vegetation_index}_Var_{window_size}/{vegetation_index}_Var_{window_size}_STL_{series_id}.csv: Variance values


### SHAP Analysis:

shap_values.csv: Individual SHAP values for each prediction

shap_importance.csv: Aggregated feature importance scores

# Performance Considerations

Use parallel processing (n_processes) for large datasets

Adjust n_iter in RandomizedSearchCV based on computational resources

Consider memory usage when processing multiple time series simultaneously

# Contact

Author: Yiling Cai

Institution: School of Geography and Remote Sensing, Guangzhou University, Guangzhou 510006, China;
            Guangdong Guodi Science Technology Co., Ltd, Guangzhou 510075, China

Email: caiyling6@alumni.sysu.edu.cn

# License

This project is licensed under the GNU GPLv3 License - see the LICENSE file for details.

# Version History

v1.0.0 (2026-01-16): Initial release
