# Overview

Code for analysing forest resilience and abrupt declines (ADs) from
satellite vegetation-index time series. The repository contains:

1. **STL decomposition** of monthly kNDVI series into trend, seasonal and
   residual components (`STL_Decomposition.py`).
2. **Resilience indicators**: sliding-window lag-1 autocorrelation (AR1, also
   called temporal autocorrelation, TAC) and variance of the STL residuals
   (`AR1_Var_calculating.py`).
3. **Driver analysis**: an XGBoost regression model of TAC preceding abrupt
   declines, interpreted with SHAP values (`XGBoost_SHAP_ADs.py`).
4. **Abrupt-decline detection and paired AD/noAD comparison** in MATLAB
   (`Forest Resilience Abrupt Decline Detection/`, documented in its own
   `Readme.md`).
   
Each Python script is self-contained, and every step is numbered in the code
comments (Step 1, Step 2, ...), matching the workflow listed at the top of the
script. 

# Project Structure

```
├── STL_Decomposition.py        # Step 1: STL decomposition
├── AR1_Var_calculating.py      # Step 2: sliding-window AR1 and variance
├── XGBoost_SHAP_ADs.py         # XGBoost model and SHAP interpretation
├── Forest Resilience Abrupt Decline Detection/   # MATLAB code (see its Readme.md)
├── tests/test_pipeline.py      # automated tests for Steps 1 and 2
├── requirements.txt            # tested package versions
├── input/                      # demo input data
│   ├── kndvi_data.csv                  # five example kNDVI series
│   └── xgb_data_sigma4div2_WS07.csv    # samples for the XGBoost model
├── output/                     # demo outputs
└── LICENSE
```


# 1. System Requirements

| | |
|---|---|
| Operating system | Any OS that runs Python (Windows, Linux, macOS). Tested on Linux (x86-64). |
| Python | >= 3.9; tested with 3.13.16 |
| Python packages | numpy 2.5.3, pandas 3.0.5, scipy 1.18.1, statsmodels 0.15.0, scikit-learn 1.9.1, xgboost 3.4.1, shap 0.53.0 (pytest 9.1.1 for the tests); see `requirements.txt` |
| MATLAB (part 4 only) | R2024a or later with the toolboxes listed in the MATLAB `Readme.md` |
| Hardware | No non-standard hardware is required. All Python scripts run on a normal desktop computer; they use several CPU cores when available. |


# 2. Installation Guide

```bash
git clone <repository-url>
cd resilience-main
python -m venv .venv
# Windows: .venv\Scripts\activate      Linux/macOS: source .venv/bin/activate
pip install -r requirements.txt
```

Typical install time on a normal desktop computer: a few minutes.

# 3. Demo

The demo data are in `input/`. All scripts read from `./input/` and write to
`./output/` by default, so they can be run without editing:

```bash
python STL_Decomposition.py
python AR1_Var_calculating.py
python XGBoost_SHAP_ADs.py
```

## Expected output

```
output/
├── trend/   kNDVI_STL_trend_<id>.csv        STL trend
├── season/  kNDVI_STL_season_<id>.csv       STL seasonal component
├── resid/   kNDVI_STL_resid_<id>.csv        STL residual
├── kNDVI_AR1_<w>/ kNDVI_AR1_<w>_STL_<id>.csv    AR1 for window size w (36, 48, 60, 72, 84 months)
├── kNDVI_Var_<w>/ kNDVI_Var_<w>_STL_<id>.csv    variance for window size w
├── xgb_random_state.txt        random seed used by the XGBoost run
├── xgb_model.json              fitted XGBoost model
├── xgb_model_summary.json      settings, selected hyperparameters, RMSE/R2, software versions
├── xgb_train_test_split.csv    training/test membership of every sample
├── xgb_cv_results.csv          cross-validation score of every hyperparameter candidate
├── shap_values.csv             SHAP value of every feature for every sample
└── shap_importance.csv         mean absolute SHAP value of every feature
```

The STL, AR1 and variance files for the five demo series (60-month window)
are provided in `output/` for comparison; re-running the two scripts
reproduces them to floating-point precision.

## Expected run time

- `STL_Decomposition.py` and `AR1_Var_calculating.py`: a few seconds for the
  demo data.
- `XGBoost_SHAP_ADs.py`: 50 hyperparameter candidates × 5 cross-validation
  folds are fitted, so the run time depends mainly on the number of CPU cores
  (`n_jobs`).

# 4. Instructions for Use

All settings are in the block marked **Step 1** of each script.

## 4.1 STL_Decomposition.py

| Variable | Meaning | Default |
|---|---|---|
| `input_file` | input table | `./input/kndvi_data.csv` |
| `output_dir` | output folder | `./output/` |
| `y_start`, `y_end` | analysis period | 1982, 2022 |
| `period` | observations per seasonal cycle (n_p) | 12 |
| `smooth_length` | length of the seasonal smoother (n_s) | 7 |
| `vegetation_index` | label used in file names | `kNDVI` |
| `n_processes` | parallel workers | 80% of the CPU cores |

A series that cannot be decomposed (for example because it contains missing
values) is listed at the end of the run and does not stop the other series.

## 4.2 AR1_Var_calculating.py

Run after `STL_Decomposition.py`.

| Variable | Meaning | Default |
|---|---|---|
| `catalog_file` | table listing the series identifiers | `./input/kndvi_data.csv` |
| `input_dir` | STL residuals | `./output/resid/` |
| `output_dir` | output folder | `./output/` |
| `window_sizes` | sliding-window lengths (months) | 36, 48, 60, 72, 84 |
| `y_start`, `y_end` | analysis period | 1982, 2022 |
| `n_processes` | parallel workers | 80% of the CPU cores |

Series whose AR1 and variance files already exist are skipped, so an
interrupted run can be resumed.

## 4.3 XGBoost_SHAP_ADs.py

| Variable | Meaning | Default |
|---|---|---|
| `input_file` | sample table | `./input/xgb_data_sigma4div2_WS07.csv` |
| `output_dir` | output folder | `./output/` |
| `id_column`, `y_column` | sample identifier and target | `PID`, `TAC_AD` |
| `x_columns` | predictors | 61 features (listed in the script) |
| `random_state` | seed of every random step | `None` (drawn at random from 0–5000 and recorded) |
| `test_size` | proportion of samples held out for testing | 0.2 |
| `n_iter` | hyperparameter candidates sampled | 50 |
| `n_splits` | cross-validation folds | 5 |
| `n_jobs` | parallel workers of the search | 40 (reduce to the number of available cores) |

**Reproducing a run.** With `random_state = None` a new seed is drawn at each
run, printed, and saved to `output/xgb_random_state.txt` and
`output/xgb_model_summary.json`. To reproduce a run exactly, set
`random_state` to the recorded value and run the script again with the same
package versions, or load the saved model directly:

```python
import xgboost as xgb
model = xgb.XGBRegressor()
model.load_model("output/xgb_model.json")
```

Each run overwrites the XGBoost and SHAP files in `output/`.

## 4.4 MATLAB code

See `Forest Resilience Abrupt Decline Detection/Readme.md`.


# Data Format

## Time series table (`input/kndvi_data.csv`)

One row per series. The first column is the series identifier (`PID`);
columns whose header is a date (`YYYY-MM-DD`) are monthly time steps in
chronological order without gaps; the other columns (`XLon`, `YLat`) are
metadata that are copied to the outputs.

```
PID,XLon,YLat,1982-01-01,1982-02-01,1982-03-01,...
3,-140.875,62.708,0.0034,0.0053,0.0132,...
1,-123.042,62.792,0.0023,0.0106,0.0110,...
```

Each output file of Steps 1 and 2 holds one series as a single column, with
the metadata rows followed by one row per date.

## Sample table (`input/xgb_data_sigma4div2_WS07.csv`)

One row per sample, with the identifier (`PID`), the target (`TAC_AD`) and
the predictor columns.

```
PID,TAC_AD,TAC_mean,FC,FL,FG,...
516142,0.3546,0.1514,57.8624,0.0288,2.0814,...
387407,0.2421,0.0808,87.4868,0.0306,0.4815,...
```

# References

- Cleveland, R. B., Cleveland, W. S., McRae, J. E. & Terpenning, I. STL: a
  seasonal-trend decomposition procedure based on loess. *J. Off. Stat.* **6**,
  3–73 (1990).
- Scheffer, M. et al. Early-warning signals for critical transitions. *Nature*
  **461**, 53–59 (2009).
- Chen, T. & Guestrin, C. XGBoost: a scalable tree boosting system. In *Proc.
  22nd ACM SIGKDD Int. Conf. Knowledge Discovery and Data Mining* 785–794
  (2016).
- Lundberg, S. M. et al. From local explanations to global understanding with
  explainable AI for trees. *Nat. Mach. Intell.* **2**, 56–67 (2020).

# Contact

Author: Yiling Cai

Institution: School of Geography and Remote Sensing, Guangzhou University, Guangzhou 510006, China;
            Guangdong Guodi Science Technology Co., Ltd, Guangzhou 510075, China

Email: caiyling6@alumni.sysu.edu.cn

# License

This project is licensed under the GNU GPLv3 License - see the LICENSE file for details.

# Version History

v1.1.0 (2026-10-10): Numbered step annotations in all Python scripts; recorded random seed and saved model.

v1.0.0 (2026-01-16): Initial release
