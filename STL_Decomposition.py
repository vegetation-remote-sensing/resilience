#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
STL Decomposition for Time Series Data

This script performs robust Seasonal-Trend decomposition using LOESS (STL;
Cleveland et al., 1990) on vegetation index time series and saves the trend,
seasonal and residual components of every series (pixel) as separate files.

Workflow:
    Step 1. Set the analysis parameters and file paths.
    Step 2. Load the input table (one row per series, one column per month).
    Step 3. Decompose every series in parallel:
            3.1 extract the series for the analysis period;
            3.2 run the robust STL decomposition;
            3.3 save the trend, seasonal and residual components.
    Step 4. Summarise how many series succeeded or failed.

Usage:
    python STL_Decomposition.py

Requirements:
    - Python >= 3.7
    - pandas >= 1.3.0
    - numpy >= 1.21.0
    - statsmodels >= 0.13.0

Reference:
    Cleveland, R. B., Cleveland, W. S., McRae, J. E. & Terpenning, I. STL: a
    seasonal-trend decomposition procedure based on loess. J. Off. Stat. 6,
    3-73 (1990).
"""

import math
import multiprocessing
from pathlib import Path

import numpy as np
import pandas as pd
from statsmodels.tsa.seasonal import STL


def calc_trend_length(period, seasonal_length):
    """
    Calculate the length of the trend smoother (n_t).

    Based on Cleveland et al., 1990: n_t is the smallest odd integer greater
    than or equal to 1.5 * n_p / (1 - 1.5 / n_s).

    Args:
        period: Period of the seasonal component (n_p)
        seasonal_length: Length of the seasonal smoother (n_s)

    Returns:
        Odd integer for the trend smoother length
    """
    nt = math.ceil((1.5 * period) / (1 - 1.5 / seasonal_length))
    return nt if nt % 2 == 1 else nt + 1


def calc_lowpass_length(period):
    """
    Calculate the length of the low-pass filter (n_l).

    Based on Cleveland et al., 1990: n_l is the smallest odd integer greater
    than or equal to n_p.

    Args:
        period: Period of the seasonal component (n_p)

    Returns:
        Odd integer for the low-pass filter length
    """
    nl = math.ceil(period)
    return nl if nl % 2 == 1 else nl + 1


def perform_stl(series, period, smooth_length, series_id):
    """
    Perform robust STL decomposition on one time series.

    Parameter settings (names follow Cleveland et al., 1990):
        n_p (period)        = period          period of the seasonal component
        n_s (seasonal)      = smooth_length   length of the seasonal smoother
        n_t (trend)         = calc_trend_length(period, smooth_length)
        n_l (low_pass)      = calc_lowpass_length(period)
        degree              = 1               locally linear fitting for the
                                              seasonal, trend and low-pass
                                              smoothers
        jump                = 1               LOESS evaluated at every point
                                              (no interpolation)
        robust              = True            robustness weights down-weight
                                              outliers
        n_i (inner loops)   = 2               statsmodels default when robust
        n_o (outer loops)   = 15              statsmodels default when robust

    Args:
        series: Input time series (no missing values, chronological order)
        period: Period of the seasonal component
        smooth_length: Length of the seasonal smoother
        series_id: Identifier of the time series (used in messages only)

    Returns:
        STL decomposition result (with .trend, .seasonal and .resid),
        or None if the decomposition failed
    """
    try:
        # Step 3.2a: STL cannot handle gaps, so reject series with missing
        # values instead of filling them implicitly.
        values = np.asarray(series, dtype=float).ravel()
        if not np.isfinite(values).all():
            raise ValueError("series contains missing or non-finite values")

        # Step 3.2b: Configure the robust STL model.
        model = STL(values,
                    period=period,
                    seasonal=smooth_length,
                    trend=calc_trend_length(period, smooth_length),
                    low_pass=calc_lowpass_length(period),
                    seasonal_deg=1, trend_deg=1, low_pass_deg=1,
                    seasonal_jump=1, trend_jump=1, low_pass_jump=1,
                    robust=True)

        # Step 3.2c: Fit the model (2 inner and 15 robustness iterations).
        return model.fit()
    except Exception as e:
        print(f"STL decomposition failed for series {series_id}: {e}")
        return None


def decompose_series(vegetation_index, series_id, data, date_start, date_end,
                     period, smooth_length, output_dir):
    """
    Decompose a single time series and save the results.

    Args:
        vegetation_index: Name of the vegetation index (used in file names)
        series_id: Identifier of the series
        data: pd.DataFrame with one column (series_id); the index holds the
              metadata fields (e.g. XLon, YLat) followed by the dates
        date_start: First date of the analysis period (YYYY-MM-DD)
        date_end: Last date of the analysis period (YYYY-MM-DD)
        period: Period of the seasonal component
        smooth_length: Length of the seasonal smoother
        output_dir: Base directory for the output files

    Returns:
        True if successful, False otherwise
    """
    try:
        # Step 3.1: Extract the time series for the analysis period.
        df_series = data.loc[date_start:date_end, series_id].astype(float)

        # Step 3.2: Perform the robust STL decomposition.
        res = perform_stl(df_series, period, smooth_length, str(series_id))
        if res is None:
            return False

        # Step 3.3: Save each component. Every output file keeps the layout
        # of the input (metadata rows followed by one row per date); only
        # the values inside the analysis period are replaced.
        components = {"trend": res.trend,
                      "season": res.seasonal,
                      "resid": res.resid}
        for name, values in components.items():
            df_out = data.copy(deep=True)
            df_out.loc[date_start:date_end, series_id] = values

            outpath = Path(output_dir) / name
            outpath.mkdir(parents=True, exist_ok=True)
            filename = "_".join([vegetation_index, "STL", name,
                                 str(series_id).zfill(7)]) + ".csv"
            df_out.to_csv(outpath / filename, sep=',')

        print(f"Successfully processed series {series_id}")
        return True
    except Exception as e:
        print(f"Failed to decompose series {series_id}: {e}")
        return False


if __name__ == '__main__':
    # ------------------------------------------------------------------
    # Step 1: Set the analysis parameters and file paths.
    # ------------------------------------------------------------------
    n_processes = max(1, int(multiprocessing.cpu_count() * 0.8))  # parallel workers
    period = 12                  # monthly data: 12 observations per year
    smooth_length = 7            # length of the seasonal smoother (odd, >= 7)
    vegetation_index = "kNDVI"   # label used in the output file names
    y_start = 1982               # first year of the analysis period
    y_end = 2022                 # last year of the analysis period
    date_start = "-".join([str(y_start), "01", "01"])
    date_end = "-".join([str(y_end), "12", "01"])

    input_file = "./input/kndvi_data.csv"
    output_dir = "./output/"

    # ------------------------------------------------------------------
    # Step 2: Load the input table.
    # Rows are series (PID); columns are XLon, YLat and one column per month
    # (YYYY-MM-DD). The table is transposed so that each series is a column.
    # ------------------------------------------------------------------
    df = pd.read_csv(input_file, index_col=0)
    df_transposed = df.T
    series_ids = df_transposed.columns.tolist()
    print(f"Loaded {len(series_ids)} series from {input_file}")

    # ------------------------------------------------------------------
    # Step 3: Decompose every series in parallel.
    # ------------------------------------------------------------------
    MyPool = multiprocessing.Pool(processes=n_processes)
    ResultsList = []
    for series_id in series_ids:
        data = df_transposed.loc[:, [series_id]].copy(deep=True)
        u = MyPool.apply_async(decompose_series,
                               (vegetation_index, series_id, data, date_start,
                                date_end, period, smooth_length, output_dir,))
        ResultsList.append((series_id, u))
    MyPool.close()
    MyPool.join()

    # ------------------------------------------------------------------
    # Step 4: Summarise the run and list the series that failed.
    # ------------------------------------------------------------------
    failed = [series_id for series_id, u in ResultsList if not u.get()]
    print(f"Finished: {len(series_ids) - len(failed)} succeeded, "
          f"{len(failed)} failed")
    if failed:
        print(f"Failed series: {failed}")
