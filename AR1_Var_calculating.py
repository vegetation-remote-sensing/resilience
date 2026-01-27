#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Autoregressive Model (AR1) Analysis for Time Series Resilience Indicators

This script calculates autoregressive lag-1 (AR1) coefficients and variance
using sliding window approach for ecosystem resilience assessment.

Author: Yiling Cai
Institution: School of Geography and Remote Sensing, Guangzhou University, Guangzhou 510006, China;
             Guangdong Guodi Science Technology Co., Ltd, Guangzhou 510075, China;
Date: 2026-01-16

Requirements:
    - Python >= 3.7
    - pandas >= 1.3.0
    - numpy >= 1.21.0
    - scipy >= 1.7.0
    - statsmodels >= 0.13.0

"""

import datetime
import gc
import multiprocessing
import os
import shutil
import time
import statsmodels.api as sm
import numpy as np
import scipy.stats, scipy.signal
import pandas as pd
from scipy.optimize import curve_fit
from pathlib import Path

### Code for Trend Estimation ###
def calculate_ar1_coefficient(x):
    """
        Calculate 1-lag autocorrelation coefficient (AR1) using OLS.

        Fits the model: X(t) = α + β*X(t-1) + ε

        Args:
            x: Time series array

        Returns:
            AR1 coefficient (β) or NaN if calculation fails
    """
    try:
        # Prepare lagged data
        x_lag = x[:-1]  # X(t-1)
        x_current = x[1:]  # X(t)

        # Add constant term for intercept
        x_lag_const = sm.add_constant(x_lag)

        # Fit OLS regression
        model = sm.OLS(x_current, x_lag_const, missing="drop")
        results = model.fit()

        # Return AR1 coefficient (slope parameter)
        return results.params[1]

    except Exception as e:
        print(f"Failed to calculate AR1 coefficient: {e}")
        return np.nan

def sliding_window_analysis(timeseries, series_id, window_size):
    """
        Apply sliding window to calculate AR1 and variance.

        Args:
            timeseries: Input time series array
            series_id: Identifier for the time series
            window_size: Size of the sliding window

        Returns:
            Tuple of (ar1_array, variance_array) with same length as input

        Note:
            Values at the edges (< half_window from start/end) are set to NaN
    """
    try:
        ln = timeseries.shape[0]
        half_window = int(window_size / 2)

        # Initialize output arrays with NaN
        ar1 = np.empty(ln)
        ar1.fill(np.nan)
        variance = np.empty(ln)
        variance.fill(np.nan)

        # Calculate for valid window positions
        for i in range(half_window, ln - half_window):
            # Extract window
            window_data  = timeseries[i - half_window: i + half_window]

            # Calculate metrics
            ar1_ = calculate_ar1_coefficient(window_data)
            variance_ = np.nanvar(window_data , ddof=1)
            ar1[i] = ar1_
            variance[i] = variance_

            del ar1_, variance_
            gc.collect()

        ar1 = np.array(ar1)
        variance = np.array(variance)
        return ar1, variance
    except Exception as e:
        print(f"Failed to sliding window analysis: {series_id} {e}")
        return None, None


def process_series(input_file, output_dir, vegetation_index, series_id, window_size, date_start, date_end):
    """
        Process a single time series to calculate AR1 and variance.

        Args:
            input_file: Path to input CSV file with residuals
            output_dir: Base output directory
            series_id: Identifier for the time series
            window_size: Size of the sliding window
            date_start: Starting date of analysis
            date_end: Ending date of analysis

        Returns:
            True if successful, False otherwise
    """
    try:
        # Define output paths
        ar1_dir = output_dir + vegetation_index + "_AR1_" + str(window_size) + "/"
        Path(ar1_dir).mkdir(parents=True, exist_ok=True)
        ar1_file = vegetation_index + "_AR1_" + str(window_size) + "_STL_" + str(series_id).zfill(7) + ".csv"

        var_dir = output_dir + vegetation_index + "_Var_" + str(window_size) + "/"
        Path(var_dir).mkdir(parents=True, exist_ok=True)
        var_file = vegetation_index + "_Var_" + str(window_size) + "_STL_" + str(series_id).zfill(7) + ".csv"

        # Skip if already processed
        if Path(ar1_file).exists() and Path(var_file).exists():
            print(f"Series {series_id} already processed, skipping")
            return True

        # Load data
        df = pd.read_csv(input_file, index_col=0)

        # Extract time series for analysis period
        df_data = df.loc[date_start:date_end, series_id].copy(deep=True)
        series_data = np.array(df_data)

        # Calculate AR1 and variance
        ar1, var = sliding_window_analysis(series_data, series_id, window_size)

        # Save results
        df_ar1 = df.copy(deep=True)
        df_ar1.loc[date_start:date_end, series_id] = ar1
        df_ar1.to_csv(ar1_file)

        df_var = df.copy(deep=True)
        df_var.loc[date_start:date_end, series_id] = var
        df_var.to_csv(var_file)

        print(f"Successfully processed series {series_id}")

        del ar1_dir, ar1_file, var_dir, var_file, df, df_data
        del series_data, ar1, var, df_ar1, df_var
        gc.collect()

        return True
    except Exception as e:
        print(f"Failed to process series {series_id}: {e}")
        return False


if __name__ == '__main__':
    # Initialize
    n_processes = 50
    window_sizes = [36, 48, 60, 72, 84]
    vegetation_index = "kNDVI"
    y_start = 1982
    y_end = 2022
    date_start = "-".join([str(y_start), "01", "01"])
    date_end = "-".join([str(y_end), "12", "01"])

    catalog_file = "./path/to/input/kndvi_data.csv"
    input_dir = "./path/to/input/residual/files/"
    output_dir = "./path/to/output/"

    # Load list of series identifiers from catalog file.
    df = pd.read_csv(catalog_file, index_col=0) #index_col = "series_ids"
    df_transposed = df.T
    series_list = df_transposed.columns.tolist()

    # Process entire dataset with multiple window sizes.
    for window_size in window_sizes:

        # Process in parallel
        MyPool = multiprocessing.Pool(processes=pp)
        ResultsList = []
        for series_id in series_list:
            input_file = input_dir + vegetation_index + "_STL_resid_" + str(series_id).zfill(7) + ".csv"
            if not Path(input_file).exists():
                print(f"Input file not found: {input_file}")
                continue

            u = MyPool.apply_async(process_series, (input_file, output_dir, vegetation_index, str(series_id), window_size,
                                                    date_start, date_end,))
            ResultsList.append(u)
        MyPool.close()
        MyPool.join()



