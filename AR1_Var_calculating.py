#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Lag-1 Autocorrelation (AR1) and Variance for Time Series Resilience Indicators

This script calculates the lag-1 autoregressive coefficient (AR1) and the
variance of STL residuals with a sliding-window approach. Rising AR1 and
variance indicate critical slowing down, i.e. declining ecosystem resilience
(Scheffer et al., 2009).

Workflow:
    Step 1. Set the analysis parameters and file paths.
    Step 2. Load the list of series identifiers.
    Step 3. For every window size, process every series in parallel:
            3.1 define the output files and skip finished series;
            3.2 load the STL residuals for the analysis period;
            3.3 slide the window along the series and calculate AR1 and
                variance in every window;
            3.4 save the AR1 and variance series.
    Step 4. Summarise how many series succeeded or failed.

Window definition:
    With half_window = window_size // 2, the value stored at time step i is
    calculated from the observations i - half_window, ..., i + half_window - 1.
    The first and the last half_window time steps are set to NaN.

Usage:
    python AR1_Var_calculating.py     (run STL_Decomposition.py first)

Requirements:
    - Python >= 3.7
    - pandas >= 1.3.0
    - numpy >= 1.21.0

Reference:
    Scheffer, M. et al. Early-warning signals for critical transitions.
    Nature 461, 53-59 (2009).
"""

import multiprocessing
from pathlib import Path

import numpy as np
import pandas as pd


def calculate_ar1_coefficient(x):
    """
    Calculate the lag-1 autoregressive coefficient (AR1) using OLS.

    Fits the model: X(t) = alpha + beta * X(t-1) + epsilon

    The slope is obtained from the closed-form OLS solution
    beta = cov(X(t-1), X(t)) / var(X(t-1)), which is identical to fitting the
    regression with statsmodels but much faster.

    Args:
        x: Time series array (one window)

    Returns:
        AR1 coefficient (beta), or NaN if it cannot be calculated
    """
    # Step 3.3a: Build the lagged pairs (X(t-1), X(t)).
    x = np.asarray(x, dtype=float)
    x_lag = x[:-1]      # X(t-1)
    x_current = x[1:]   # X(t)

    # Step 3.3b: Drop the pairs that contain a missing value.
    valid = np.isfinite(x_lag) & np.isfinite(x_current)
    if valid.sum() < 3:
        return np.nan
    x_lag = x_lag[valid]
    x_current = x_current[valid]

    # Step 3.3c: OLS slope with intercept = covariance / variance.
    x_lag_centered = x_lag - x_lag.mean()
    denominator = np.dot(x_lag_centered, x_lag_centered)
    if denominator == 0:  # constant window: slope is undefined
        return np.nan
    return np.dot(x_lag_centered, x_current - x_current.mean()) / denominator


def sliding_window_analysis(timeseries, series_id, window_size):
    """
    Apply a sliding window to calculate AR1 and variance.

    Args:
        timeseries: Input time series array (STL residuals)
        series_id: Identifier of the time series (used in messages only)
        window_size: Size of the sliding window (number of time steps)

    Returns:
        Tuple of (ar1_array, variance_array) with the same length as the input

    Note:
        Values at the edges (< half_window from start/end) are set to NaN
    """
    timeseries = np.asarray(timeseries, dtype=float)
    ln = timeseries.shape[0]
    half_window = int(window_size / 2)

    # Step 3.3-1: Initialise the output arrays with NaN.
    ar1 = np.full(ln, np.nan)
    variance = np.full(ln, np.nan)

    # Step 3.3-2: Move the window along the series.
    for i in range(half_window, ln - half_window):
        # Extract the window around time step i.
        window_data = timeseries[i - half_window: i + half_window]

        # Calculate the two resilience indicators in this window.
        ar1[i] = calculate_ar1_coefficient(window_data)
        if np.isfinite(window_data).sum() > 1:
            variance[i] = np.nanvar(window_data, ddof=1)  # sample variance

    return ar1, variance


def process_series(input_file, output_dir, vegetation_index, series_id,
                   window_size, date_start, date_end):
    """
    Process a single time series to calculate AR1 and variance.

    Args:
        input_file: Path to the input CSV file with the STL residuals
        output_dir: Base output directory
        vegetation_index: Name of the vegetation index (used in file names)
        series_id: Identifier of the time series
        window_size: Size of the sliding window
        date_start: First date of the analysis period (YYYY-MM-DD)
        date_end: Last date of the analysis period (YYYY-MM-DD)

    Returns:
        True if successful, False otherwise
    """
    try:
        # Step 3.1: Define the output files.
        ar1_name = "_".join([vegetation_index, "AR1", str(window_size)])
        ar1_dir = Path(output_dir) / ar1_name
        ar1_dir.mkdir(parents=True, exist_ok=True)
        ar1_file = ar1_dir / (ar1_name + "_STL_" + str(series_id).zfill(7) + ".csv")

        var_name = "_".join([vegetation_index, "Var", str(window_size)])
        var_dir = Path(output_dir) / var_name
        var_dir.mkdir(parents=True, exist_ok=True)
        var_file = var_dir / (var_name + "_STL_" + str(series_id).zfill(7) + ".csv")

        # Skip the series if both output files already exist.
        if ar1_file.exists() and var_file.exists():
            print(f"Series {series_id} already processed, skipping")
            return True

        # Step 3.2: Load the STL residuals for the analysis period.
        df = pd.read_csv(input_file, index_col=0)
        series_data = df.loc[date_start:date_end, series_id].to_numpy(dtype=float)

        # Step 3.3: Calculate AR1 and variance in a sliding window.
        ar1, var = sliding_window_analysis(series_data, series_id, window_size)

        # Step 3.4: Save the results. Each output file keeps the layout of
        # the input (metadata rows followed by one row per date); only the
        # values inside the analysis period are replaced.
        df_ar1 = df.copy(deep=True)
        df_ar1.loc[date_start:date_end, series_id] = ar1
        df_ar1.to_csv(ar1_file)

        df_var = df.copy(deep=True)
        df_var.loc[date_start:date_end, series_id] = var
        df_var.to_csv(var_file)

        print(f"Successfully processed series {series_id} (window {window_size})")
        return True
    except Exception as e:
        print(f"Failed to process series {series_id}: {e}")
        return False


if __name__ == '__main__':
    # ------------------------------------------------------------------
    # Step 1: Set the analysis parameters and file paths.
    # ------------------------------------------------------------------
    n_processes = max(1, int(multiprocessing.cpu_count() * 0.8))  # parallel workers
    window_sizes = [36, 48, 60, 72, 84]  # sliding-window lengths (months)
    vegetation_index = "kNDVI"           # label used in the file names
    y_start = 1982                       # first year of the analysis period
    y_end = 2022                         # last year of the analysis period
    date_start = "-".join([str(y_start), "01", "01"])
    date_end = "-".join([str(y_end), "12", "01"])

    catalog_file = "./input/kndvi_data.csv"   # table listing the series IDs
    input_dir = "./output/resid/"             # STL residuals (previous script)
    output_dir = "./output/"

    # ------------------------------------------------------------------
    # Step 2: Load the list of series identifiers from the catalog file.
    # ------------------------------------------------------------------
    df = pd.read_csv(catalog_file, index_col=0)
    series_list = df.index.tolist()
    print(f"Loaded {len(series_list)} series from {catalog_file}")

    # ------------------------------------------------------------------
    # Step 3: For every window size, process every series in parallel.
    # ------------------------------------------------------------------
    n_failed = 0
    for window_size in window_sizes:
        MyPool = multiprocessing.Pool(processes=n_processes)
        ResultsList = []
        for series_id in series_list:
            # Locate the residual file written by STL_Decomposition.py.
            input_file = Path(input_dir) / (vegetation_index + "_STL_resid_"
                                            + str(series_id).zfill(7) + ".csv")
            if not input_file.exists():
                print(f"Input file not found: {input_file}")
                n_failed += 1
                continue

            u = MyPool.apply_async(process_series,
                                   (input_file, output_dir, vegetation_index,
                                    str(series_id), window_size, date_start,
                                    date_end,))
            ResultsList.append((series_id, u))
        MyPool.close()
        MyPool.join()

        # --------------------------------------------------------------
        # Step 4: Summarise the run for this window size.
        # --------------------------------------------------------------
        failed = [series_id for series_id, u in ResultsList if not u.get()]
        n_failed += len(failed)
        print(f"Window {window_size}: {len(ResultsList) - len(failed)} "
              f"succeeded, {len(failed)} failed")
        if failed:
            print(f"Failed series: {failed}")

    print(f"Finished: {n_failed} failures in total")
