## Overview
This code implements a spatially explicit case-control methodology for detecting abrupt decline (AD) events in forest vegetation and analyzing whether declining resilience serves as an early warning signal. The analysis uses satellite-derived vegetation indices to compare resilience trends between pixels that experienced abrupt decline and neighboring stable control pixels.

## Citation
If you use this code, please cite:
[paper citation - to be added upon publication]

## Requirements
MATLAB R2020a or later
Image Processing Toolbox
Mapping Toolbox
Statistics and Machine Learning Toolbox
Parallel Computing Toolbox (optional, for performance)

## Methods
# 1. Abrupt Decline (AD) Event Detection
AD events are identified using threshold-based anomaly detection:
a. Calculate mean growing-season kNDVI (May-September) for each pixel (2000-2022)
b. Remove linear trends from the time series
c. Define baseline statistics (μ, σ) from baseline period (2000-2009)
d. Classify AD pixels: kNDVI < μ - nσ, where n ∈ [1,6]
e. Classify noAD pixels: kNDVI > μ - nσ/m, where m ∈ [2,6]
# 2. Resilience Trend Analysis
Temporal changes in vegetation resilience are quantified using:
a. AR(1) coefficient (Temporal Autocorrelation - TAC) as resilience indicator
b. Sen's slope estimator (δTAC - non-parametric trend magnitude)
c. Mann-Kendall test (statistical significance, α = 0.05)
d. Rolling window analysis (cumulative from baseline to year t-1)
Output classification:
a. -2: Significant declining resilience
b. -1: Non-significant decline
c. 0: No trend
d. +1: Non-significant increase
e. +2: Significant increasing resilience
# 3. Paired AD-noAD Metric Extraction
For each AD pixel, the code:
a. Identifies year of first AD occurrence (t)
b. Extracts pre-event resilience metrics at t-1 (δTAC and TAC)
c. Defines control group from spatial neighborhood (5×5, 7×7, or 9×9 window)
d. Calculates mean resilience metrics for noAD pixels
e. Creates matched case-control pairs for comparison
Purpose: Test whether declining resilience (AD_δTAC > noAD_δTAC) precedes abrupt decline events.

## Installation
% Add code to MATLAB path
addpath(genpath('./'));
savepath;

## Data Preparation

# Input Requirements
1. Growing Season kNDVI Maps
Format: GeoTIFF (.tif)
Coverage: Annual maps (2000-2022)
Naming: kNDVI_GS_05_09_YYYY.tif
Values: -1 to +1 (normalized)
2. AR(1) Resilience Maps
Format: GeoTIFF (.tif)
Coverage: Annual maps (2000-2022)
Naming: kNDVI_AR1_60_STL_YYYY.tif
Values: -1 to +1 (autocorrelation coefficient)
3. Spatial Reference File
Format: GeoTIFF (.tif)
Content: Land mask (0 = valid pixels)
Must match spatial dimensions of input data

# Directory Structure
project_root/
├── data/
│   ├── kNDVI_growing_season/
│   │   └── kNDVI_GS_05_09_YYYY.tif (2000-2022)
│   ├── AR1/
│   │   └── kNDVI_AR1_60_STL_YYYY.tif (2000-2022)
│   └── reference/
│       └── land_mask.tif
├── detect_abrupt_decline_events.m
├── calculate_resilience_trends.m
├── extract_paired_AD_noAD_resilience.m
└── example_run.m

## Usage
# Quick Start
% Run complete workflow
example_run

# Output Structure
results/
├── AD_noAD_classification/
│   └── AD_noAD_classification_sigma[n]div[m]_YYYY.tif
├── Sen_AR1/
│   └── Sen_kNDVI_AR1_60_STL_2000_YYYY.tif
├── AD_year/
│   └── AD_sigma[n]div[m]_WS[size].tif
├── AD_deltaTAC/
│   └── AD_sigma[n]div[m]_WS[size]_Sen_kNDVI_AR1_60_STL.tif
├── noAD_deltaTAC/
│   └── noAD_sigma[n]div[m]_WS[size]_Sen_kNDVI_AR1_60_STL.tif
└── noAD_count/
    └── noAD_sigma[n]div[m]_WS[size].tif

# File naming:
[n]: Sigma threshold (1-6)
[m]: Buffer parameter (2-6)
[size]: Window size (05, 07, 09)
YYYY: Year

## Contact
Yiling Cai
School of Geography and Remote Sensing, Guangzhou University, Guangzhou 510006, China; Guangdong Guodi Science Technology Co., Ltd, Guangzhou 510075, China;
caiyling6@alumni.sysu.edu.cn

