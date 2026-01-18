#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
STL Decomposition for Time Series Data

This script performs Seasonal-Trend decomposition using LOESS (STL) on
vegetation index time series data.

Author: Yiling Cai
Institution: School of Geography and Remote Sensing, Guangzhou University, Guangzhou 510006, China;
             Guangdong Guodi Science Technology Co., Ltd, Guangzhou 510075, China;
Date: 2026-01-16

Requirements:
    - Python >= 3.7
    - pandas >= 1.3.0
    - numpy >= 1.21.0
    - statsmodels >= 0.13.0

"""

import gc
import os
import multiprocessing
import shutil

import pandas as pd
from statsmodels.tsa.seasonal import STL
import numpy as np
from pathlib import Path

def calc_trend_length(period, seasonal_length):
    """
    Calcualte the length of the trend smoother.

    based on Cleveland et al., 1990.

    Args:
            period: Period of seasonal component
            seasonal_length: Length of seasonal smoother

    Returns:
        Odd integer for trend smoother length
    """
    nt = (1.5 * period) / (1 - 1.5 * (1 / seasonal_length)) + 1 #Force fractions to be rounded up
    if int(nt) % 2. == 1:
        return int(nt)
    elif int(nt) % 2. == 0:
        return int(nt) + 1

def calc_lowpass_length(period):
    """
       Calculate the length of the low-pass filter.

       Based on Cleveland et al., 1990.

       Args:
           period: Period of seasonal component

       Returns:
           Odd integer for low-pass filter length
    """
    if int(period) % 2. == 1:
        return int(period)
    elif int(period) % 2. == 0:
        return int(period) + 1

def perform_stl(series, period, smooth_length, series_id):
    ### Perform robust STL decomposition on a time series. ###
    # series:              # Input time series data
    # period = 12              # period of seasonal component
    # seasonal_length = 7              # length of seasonal smoother
    # series_id:              # Identifier for the time series
    # nt = calc_trend_length(period, seasonal_length)  # length of trend smoother
    # nl = calc_lowpass_length(period)     # length of low-pass filter
    # isdeg = 1           # Degree of locally-fitted polynomial in seasonal smoothing.
    # itdeg = 1           # Degree of locally-fitted polynomial in trend smoothing.
    # ildeg = 1           # Degree of locally-fitted polynomial in low-pass smoothing.
    # nsjump = None       # Skipping value for seasonal smoothing.
    # ntjump = 1          # Skipping value for trend smoothing. If None, ntjump= 0.1*nt
    # nljump = 1          # Skipping value for low-pass smoothing. If None, nljump= 0.1*nl
    # robust = True       # Flag indicating whether robust fitting should be performed.
    # ni = 1              # Number of loops for updating the seasonal and trend  components.
    # no = 3              # Number of iterations of robust fitting. The value of no should
    #                       be a nonnegative integer. If the data are well behaved without
    #                       outliers, then robustness iterations are not needed. In this case
    #                       set no=0, and set ni=2-5 depending on how much security you want
    #                       that the seasonal-trend looping converges. If outliers are present
    #                       then no=3 is a very secure value unless the outliers are radical,
    #                       in which case no=5 or even 10 might be better. If no>0 then set ni
    #                       to 1 or 2. If None, then no is set to 15 for robust fitting,
    #                       to 0 otherwise.
    ### Returns:             STL decomposition result or None if failed ###
    try:
        res = STL(series, period, seasonal=smooth_length, trend=calc_trend_length(period, smooth_length),
                  low_pass=calc_lowpass_length(period), seasonal_deg=1, trend_deg=1, low_pass_deg=1,
                  seasonal_jump=1,trend_jump=1, low_pass_jump=1, robust=True)
        return res.fit()
    except Exception as e:
        print(f"STL decomposition failed for series {series_id}: {e}")
        return None

def decompose_series(series_id, data, date_start, date_end, period, smooth_length, output_dir):
    """
        Decompose a single time series and save results.

        Args:
            series_id: str, Identifier for the series
            data: pd.DataFrame, Input dataframe containing the time series
            date_start: Starting date of analysis
            date_end: Ending date of analysis
            period: period of seasonal component
            seasonal_length: length of seasonal smoother
            output_dir: str, Base directory for output files

        Returns:
            True if successful, False otherwise
    """
    try:
        # Extract time series for specified date range
        df_series = data.loc[date_start:date_end, [series_id]].copy(deep=True)

        df_trend = data.copy(deep=True)
        df_trend.loc[date_start:date_end, series_id] = np.nan
        df_seasonal = data.copy(deep=True)
        df_seasonal.loc[date_start:date_end, series_id] = np.nan
        df_resid = data.copy(deep=True)
        df_resid.loc[date_start:date_end, series_id] = np.nan

        # Perform STL decomposition
        res = perform_stl(df_series, period, smooth_length, str(series_id))

        df_trend.loc[date_start:date_end, series_id] = res.trend
        df_seasonal.loc[date_start:date_end, series_id] = res.seasonal
        df_resid.loc[date_start:date_end, series_id] = res.resid

        # Save results
        outpath_trend = output_dir + "trend/"
        Path(outpath_trend).mkdir(parents=True, exist_ok=True)
        filename_trend = "_".join(["trend", str(series_id).zfill(7)]) + ".csv"
        df_trend.to_csv(outpath_trend + filename_trend, sep=',')

        outpath_season = output_dir + "\\season\\"
        Path(outpath_season).mkdir(parents=True, exist_ok=True)
        filename_season = "_".join(["season", str(series_id).zfill(7)]) + ".csv"
        df_seasonal.to_csv(outpath_season + filename_season, sep=',')

        outpath_resid = output_dir + "\\resid\\"
        Path(outpath_resid).mkdir(parents=True, exist_ok=True)
        filename_resid = "_".join(["resid", str(series_id).zfill(7)]) + ".csv"
        df_resid.to_csv(outpath_resid + filename_resid, sep=',')

        print(f"Successfully processed series {series_id}")
        del df_series, df_trend, df_seasonal, df_resid, res
        del outpath_trend, filename_trend, outpath_season, filename_season, outpath_resid, outpath_resid
        gc.collect()
        return True
    except Exception as e:
        print(f"Failed to decompose series {series_id}: {e}")
        return False

if __name__ == '__main__':
    # Initialize
    n_processes = 50
    period = 12
    smooth_length = 7
    y_start = 1982
    y_end = 2022
    date_start = "-".join([str(y1), "01", "01"])
    date_end = "-".join([str(y2), "12", "01"])

    input_file = "./path/to/input/data.csv"
    output_dir = "./path/to/output/"

    # Load input data
    df = pd.read_csv(input_file, index_col=0)
    df_transposed = df.T
    series_ids = df_transposed.columns.tolist()

    # Process in parallel
    MyPool = multiprocessing.Pool(processes=n_processes)
    ResultsList = []
    for series_id in series_ids:
            u = MyPool.apply_async(decompose_series, (series_id, data, date_start, date_end,
                                                      period, smooth_length, output_dir,))
            ResultsList.append(u)
    MyPool.close()
    MyPool.join()
