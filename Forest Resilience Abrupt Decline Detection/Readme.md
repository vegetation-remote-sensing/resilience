# Overview
This code implements a spatially explicit case-control methodology for detecting abrupt decline (AD) events in forest vegetation and analyzing whether declining resilience serves as an early warning signal. The analysis uses satellite-derived vegetation indices to compare resilience trends between pixels that experienced abrupt decline and neighboring stable control pixels.

Every function lists its workflow at the top of the file, and the same step numbers (Step 1, Step 2, ...) mark the corresponding code.

# Requirements
MATLAB R2020a or later  
Image Processing Toolbox  
Mapping Toolbox  
Statistics and Machine Learning Toolbox  
Parallel Computing Toolbox (optional, for performance)  

No non-standard hardware is required.

# Workflow
`example_run.m` runs the three steps in order:

| Step | Function | Input | Output |
|---|---|---|---|
| 1 | `detect_abrupt_decline_events` | growing-season kNDVI maps | AD/noAD classification maps |
| 2 | `calculate_resilience_trends` | annual AR(1) maps | Sen's slope (δTAC) and Mann-Kendall maps |
| 3 | `extract_paired_AD_noAD_resilience` | outputs of Steps 1 and 2, AR(1) maps | paired AD and noAD resilience maps |

# Methods
## 1. Abrupt Decline (AD) Event Detection
AD events are identified using threshold-based anomaly detection:  
a. Read the mean growing-season kNDVI (May-September) of each pixel and year;  
b. Remove the linear trend of each pixel, fitted over the whole period (baseline start to analysis end);  
c. Calculate the baseline statistics (μ, σ) from the baseline period (2000-2009); σ is the population standard deviation (normalised by N);  
d. Classify AD pixels (value 1): kNDVI < μ - nσ, where n ∈ [1,6];  
e. Classify noAD pixels (value 0): kNDVI > μ - nσ/m, where m ∈ [2,6];  
f. Pixels between the two thresholds and masked pixels are undefined (NaN).  
## 2. Resilience Trend Analysis
Temporal changes in vegetation resilience are quantified using:  
a. AR(1) coefficient (Temporal Autocorrelation - TAC) as resilience indicator;  
b. Sen's slope estimator (δTAC - non-parametric trend magnitude, per year);  
c. Mann-Kendall test (two-sided, α = 0.05, with continuity correction and correction for tied values);  
d. Growing window: for each target year t-1, the trend is calculated from all years between the baseline start and t-1.  
At least 3 valid years are required; missing years are skipped and the remaining values keep their original time positions.
### Output classification:  
a. -2: Significant declining resilience;  
b. -1: Non-significant decline;  
c. 0: No trend;  
d. +1: Non-significant increase;  
e. +2: Significant increasing resilience.  
## 3. Paired AD-noAD Metric Extraction
For each AD pixel, the code:  
a. Identifies year of first AD occurrence (t);  
b. Extracts pre-event resilience metrics at t-1 (δTAC and TAC);  
c. Defines control group from spatial neighborhood (5×5, 7×7, or 9×9 window, truncated at the grid edges): the pixels classified as noAD in year t;  
d. Calculates mean resilience metrics at t-1 for the noAD pixels;  
e. Creates matched case-control pairs for comparison;  
Purpose: Test whether declining resilience (AD_δTAC > noAD_δTAC) precedes abrupt decline events.  

# Installation
Add code to MATLAB path:  
addpath(genpath('./'));  
savepath;

# Demo
Demo data are provided in `input/`: a 26 × 26 pixel subset (0.0833° grid, 562 valid pixels) with growing-season kNDVI maps (1982-2022), annual AR(1) maps (1985-2019) and the land mask. The expected outputs are provided in `results/`.

To run the demo, set the MATLAB current folder to this directory and run:

example_run

The script uses the baseline period 2000-2009, the analysis period 2010-2020, n = 1-6, m = 2-6 and neighborhood sizes 5, 7 and 9. It writes 960 GeoTIFF files to `results/` and reproduces the provided values on all valid pixels. The provided result files were exported with ArcGIS: pixels without a value are stored as NoData instead of NaN, classification maps are stored as integers, and the auxiliary files (.tfw, .aux.xml, .xml, .vat.dbf) are ArcGIS sidecar files that the code does not need.

Expected run time on a normal desktop computer: a few minutes or less.

# Data Preparation

## Input Requirements
All maps must share the grid of the spatial reference file. A missing or mismatched input file stops the run with an error.
### 1. Growing Season kNDVI Maps  
Format: GeoTIFF (.tif)  
Coverage: Annual maps (at least baseline start to analysis end, e.g. 2000-2020)  
Naming: kNDVI_GS_05_09_YYYY.tif  
Values: -1 to +1 (normalized)  
### 2. AR(1) Resilience Maps  
Format: GeoTIFF (.tif)  
Coverage: Annual maps (at least baseline start to analysis end - 1, e.g. 2000-2019)  
Naming: kNDVI_AR1_60_STL_YYYY.tif  
Values: -1 to +1 (autocorrelation coefficient)  
### 3. Spatial Reference File  
Format: GeoTIFF (.tif)  
Content: Land mask (0 = valid pixels)  
Must match spatial dimensions of input data

## Directory Structure
Forest Resilience Abrupt Decline Detection/  
├── input/  
│   ├── kNDVI_GS/  
│   │   └── kNDVI_GS_05_09_YYYY.tif  
│   ├── kNDVI_AR1/  
│   │   └── kNDVI_AR1_60_STL_YYYY.tif  
│   └── OL/  
│       └── LC_MASK.tif  
├── results/  (expected outputs of the demo)  
├── detect_abrupt_decline_events.m  
├── calculate_resilience_trends.m  
├── extract_paired_AD_noAD_resilience.m  
└── example_run.m  

# Usage
## Quick Start
% Run complete workflow  
example_run.m

All settings (paths, periods, thresholds, neighborhood sizes) are defined in the configuration blocks of `example_run.m`. Optional settings that are omitted take the default values listed in the header of each function.

## Output Structure
results/  
├── detrended_kNDVI_growing_season/  
│   └── kNDVI_GS_05_09_detrended_YYYY.tif  
├── sigma_thresholds/  
│   ├── threshold_sigma[n]_2000_2009.tif  
│   └── threshold_sigma[n]div[m]_2000_2009.tif  
├── AD_noAD_classification/  
│   └── AD_noAD_classification_sigma[n]div[m]_YYYY.tif  
├── Sen_AR1/  
│   └── Sen_kNDVI_AR1_60_STL_2000_YYYY.tif  
├── MK_AR1/  
│   └── MK_kNDVI_AR1_60_STL_2000_YYYY.tif  
├── SenMK_AR1/  
│   └── SenMK_kNDVI_AR1_60_STL_2000_YYYY.tif  
├── AD_year/  
│   └── AD_sigma[n]div[m]_WS[size].tif  
├── AD_deltaTAC/  
│   └── AD_sigma[n]div[m]_WS[size]_Sen_kNDVI_AR1_60_STL.tif  
├── AD_TAC/  
│   └── AD_sigma[n]div[m]_WS[size]_kNDVI_AR1_60_STL.tif  
├── noAD_count/  
│   └── noAD_sigma[n]div[m]_WS[size].tif  
├── noAD_deltaTAC/  
│   └── noAD_sigma[n]div[m]_WS[size]_Sen_kNDVI_AR1_60_STL.tif  
└── noAD_TAC/  
    └── noAD_sigma[n]div[m]_WS[size]_kNDVI_AR1_60_STL.tif  

## File naming:
[n]: Sigma threshold (1-6)  
[m]: Buffer parameter (2-6)  
[size]: Window size (05, 07, 09)  
YYYY: Year

# References
Sen, P. K. Estimates of the regression coefficient based on Kendall's tau. J. Am. Stat. Assoc. 63, 1379-1389 (1968).  
Mann, H. B. Nonparametric tests against trend. Econometrica 13, 245-259 (1945).  
Kendall, M. G. Rank Correlation Methods (Griffin, 1975).

# Contact
Yiling Cai  
School of Geography and Remote Sensing, Guangzhou University, Guangzhou 510006, China; Guangdong Guodi Science Technology Co., Ltd, Guangzhou 510075, China  
caiyling6@alumni.sysu.edu.cn

# Version History
v1.1 (2026-10-10): Demo data and expected results added; numbered step annotations; default settings now applied; missing input files stop the run; faster trend calculation; progress messages fixed.  
v1.0 (2026-01-16): Initial release
