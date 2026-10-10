#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
XGBoost Regression with SHAP Interpretation

This script trains an XGBoost regression model that predicts the temporal
autocorrelation (TAC) preceding abrupt declines from environmental drivers,
evaluates it on an independent test set and explains it with SHAP values.

Workflow:
    Step 1. Set the parameters, file paths and random seed.
    Step 2. Load the data and select the features and the target.
    Step 3. Split the data into a training set and a test set.
    Step 4. Search the hyperparameters by cross-validation on the training
            set only (the test set is never used for model selection).
    Step 5. Evaluate the final model on the training and the test set.
    Step 6. Save the model and everything needed to reproduce it.
    Step 7. Calculate SHAP values of the final model for all samples.
    Step 8. Save the SHAP values and the SHAP feature importance.

Reproducibility:
    All random steps (train/test split, cross-validation folds, parameter
    sampling, XGBoost row/column subsampling) use one seed. By default the
    seed is drawn at random from 0-5000 and recorded in
    "xgb_random_state.txt" and "xgb_model_summary.json". To reproduce a run,
    set random_state in Step 1 to the recorded value. The fitted model, its
    hyperparameters, the train/test membership of every sample and the
    software versions are also written to the output directory, so the
    reported results can be regenerated either by re-running this script
    with the recorded seed or by loading the saved model ("xgb_model.json").

Usage:
    python XGBoost_SHAP_ADs.py

Requirements:
    - Python >= 3.8
    - pandas >= 1.3.0
    - numpy >= 1.21.0
    - xgboost >= 1.6.0
    - scikit-learn >= 1.0
    - shap >= 0.41.0

References:
    Chen, T. & Guestrin, C. XGBoost: a scalable tree boosting system. In
    Proc. 22nd ACM SIGKDD 785-794 (2016).
    Lundberg, S. M. et al. From local explanations to global understanding
    with explainable AI for trees. Nat. Mach. Intell. 2, 56-67 (2020).
"""

import json
import platform
from pathlib import Path

import numpy as np
import pandas as pd
import shap
import sklearn
import xgboost as xgb
from sklearn.metrics import mean_squared_error, r2_score
from sklearn.model_selection import KFold, RandomizedSearchCV, train_test_split


def evaluate(model, x, y):
    """
    Calculate the root-mean-square error (RMSE) and R2 of a fitted model.

    Args:
        model: Fitted regression model
        x: Feature matrix
        y: Observed target values

    Returns:
        Tuple of (rmse, r2)
    """
    pred = model.predict(x)
    return float(np.sqrt(mean_squared_error(y, pred))), float(r2_score(y, pred))


if __name__ == '__main__':
    # ------------------------------------------------------------------
    # Step 1: Set the parameters, file paths and random seed.
    # ------------------------------------------------------------------
    random_state = None   # seed for every random step:
                          #   None    -> draw a new seed at random (0-5000)
                          #   integer -> reuse a recorded seed to reproduce
                          #              an earlier run exactly
    test_size = 0.2       # proportion of samples held out for testing
    n_iter = 50           # number of hyperparameter combinations sampled
    n_splits = 5          # number of cross-validation folds
    n_jobs = 40            # parallel workers of the hyperparameter search

    input_file = "./input/xgb_data_sigma4div2_WS07.csv"
    output_dir = Path("./output/")
    output_dir.mkdir(parents=True, exist_ok=True)

    # Draw the seed (if not given) and record it immediately, so that the
    # run can be reproduced even if it is interrupted later.
    if random_state is None:
        random_state = int(np.random.default_rng().integers(0, 5001))
        seed_source = "drawn at random from 0-5000"
    else:
        seed_source = "set by the user"
    (output_dir / "xgb_random_state.txt").write_text(f"{random_state}\n")
    print(f"random_state = {random_state} ({seed_source})")

    id_column = "PID"     # sample (pixel) identifier
    y_column = "TAC_AD"   # target variable
    x_columns = [  # predictors (61 features)
        "TAC_mean",
        "FC", "FL", "FG", "IF", "RF", "PF", "CH", "TD", "FA",
        "LNC", "LPC", "SLA", "LDMC", "SD", "NSR", "RD",
        "PDSI_mean", "sen_PDSI", "PDSI_cv",
        "PRE_mean", "sen_PRE", "PRE_cv",
        "TMP_mean", "sen_TMP", "TMP_cv",
        "Tmax_mean", "sen_Tmax", "Tmax_cv",
        "Tmin_mean", "sen_Tmin", "Tmin_cv",
        "SOL_mean", "sen_SOL", "SOL_cv",
        "LST_mean", "sen_LST", "LST_cv",
        "ET_mean", "sen_ET", "ET_cv",
        "WTD",
        "VPD_mean", "sen_VPD", "VPD_cv",
        "SM_mean", "sen_SM", "SM_cv",
        "CEC", "BDOD", "Clay", "Sand", "PH", "TN", "SOC", "OCD", "AWC",
        "Fire",
        "ELEV", "SLP", "ASP",
    ]

    # Hyperparameter search space.
    param_dist = {
        'n_estimators': np.arange(50, 1000, 10),
        'max_depth': np.arange(3, 8),
        'min_child_weight': np.arange(1, 10),
        'gamma': loguniform(1e-8, 1e-2),
        'subsample': np.linspace(0.5, 1.0, 5),
        'colsample_bytree': np.linspace(0.5, 1.0, 5),
        'reg_alpha': loguniform(1e-3, 1e1),
        'reg_lambda': loguniform(1e-3, 1e1),
        'learning_rate': [0.01, 0.1],
    }
    # random_state, n_estimators, max_depth, min_child_weight, gamma, subsample, colsample_bytree, reg_alpha, reg_lambda, learning_rate
    # 4016, 240, 5, 1, 0.0031,  0.8, 0.8, 0.02, 3, 0.1
    # 4568, 150, 6, 3, 2.7e-07, 0.7, 0.8, 0.4, 8, 0.1
    # 4346, 280, 5, 5, 0.006,   0.8, 0.8, 0.01, 7, 0.1
    # 4398, 110, 5, 8, 0.013,   0.8, 0.8, 0.02, 3, 0.1
    # 3398, 180, 4, 6, 0.02,    0.8, 1, 0.3, 0.6, 0.1
    # 4914, 160, 6, 9, 0.0061,  0.8, 0.8, 0.07, 0.2, 0.1
    # 71, 170, 5, 4, 0.0052,    0.8, 0.8, 0.3, 3, 0.1
    # 1593, 130, 5, 4, 0.0067,  0.8, 0.8, 0.1, 0.2, 0.1
    # 1617, 160, 6, 6, 1e-06,   0.8, 0.8, 1, 3, 0.1
    # 3022, 110, 7, 8, 0.0001,  0.8, 0.8, 0.3, 0.3, 0.1
    # 672, 90, 5, 1, 0.0053,    0.8, 0.8, 0, 1, 0.1

    # ------------------------------------------------------------------
    # Step 2: Load the data and select the features and the target.
    # ------------------------------------------------------------------
    data = pd.read_csv(input_file)
    missing = [c for c in [id_column, y_column] + x_columns if c not in data.columns]
    if missing:
        raise KeyError(f"Columns not found in {input_file}: {missing}")

    x_all = data[x_columns]
    y_all = data[y_column]
    print(f"Loaded {len(data)} samples and {len(x_columns)} features")

    # ------------------------------------------------------------------
    # Step 3: Split the data into a training set and a test set.
    # The test set is set aside here and is used only in Step 5.
    # ------------------------------------------------------------------
    x_train, x_test, y_train, y_test = train_test_split(
        x_all, y_all, test_size=test_size, random_state=random_state)
    print(f"Training samples: {len(x_train)}, test samples: {len(x_test)}")

    # ------------------------------------------------------------------
    # Step 4: Search the hyperparameters on the training set only.
    # Each candidate is scored by k-fold cross-validation within the
    # training set; the best candidate is then refitted on the whole
    # training set (refit=True) and becomes the final model.
    # ------------------------------------------------------------------
    xgb_model = xgb.XGBRegressor(objective="reg:squarederror",
                                 tree_method="hist",
                                 random_state=random_state,
                                 n_jobs=1)
    cv = KFold(n_splits=n_splits, shuffle=True, random_state=random_state)
    random_search = RandomizedSearchCV(xgb_model,
                                       param_distributions=param_dist,
                                       n_iter=n_iter,
                                       scoring='neg_mean_squared_error',
                                       cv=cv,
                                       refit=True,
                                       n_jobs=n_jobs,
                                       verbose=1,
                                       random_state=random_state)
    random_search.fit(x_train, y_train)

    best_model = random_search.best_estimator_
    best_params = {k: v.item() if hasattr(v, "item") else v
                   for k, v in random_search.best_params_.items()}
    cv_rmse = float(np.sqrt(-random_search.best_score_))
    print("Best parameters found: ", best_params)
    print(f"Cross-validated RMSE (training set): {cv_rmse}")

    # ------------------------------------------------------------------
    # Step 5: Evaluate the final model.
    # The test-set scores measure performance on data that played no part
    # in training or hyperparameter selection.
    # ------------------------------------------------------------------
    train_rmse, train_r2 = evaluate(best_model, x_train, y_train)
    test_rmse, test_r2 = evaluate(best_model, x_test, y_test)
    print(f"Train RMSE: {train_rmse}, Test RMSE: {test_rmse}")
    print(f"Train R2: {train_r2}, Test R2: {test_r2}")

    # ------------------------------------------------------------------
    # Step 6: Save the model and everything needed to reproduce it.
    # ------------------------------------------------------------------
    # 6.1 The fitted model (reload with xgb.XGBRegressor().load_model()).
    model_file = output_dir / "xgb_model.json"
    best_model.save_model(model_file)

    # 6.2 Train/test membership of every sample.
    split = pd.DataFrame({id_column: data[id_column], "set": "train"})
    split.loc[x_test.index, "set"] = "test"
    split.to_csv(output_dir / "xgb_train_test_split.csv", index=False)

    # 6.3 Scores of every hyperparameter combination that was tried.
    pd.DataFrame(random_search.cv_results_).to_csv(
        output_dir / "xgb_cv_results.csv", index=False)

    # 6.4 Settings, selected hyperparameters, performance and versions.
    summary = {
        "input_file": str(input_file),
        "target": y_column,
        "features": x_columns,
        "n_samples": int(len(data)),
        "n_train": int(len(x_train)),
        "n_test": int(len(x_test)),
        "random_state": random_state,
        "random_state_source": seed_source,
        "test_size": test_size,
        "n_iter": n_iter,
        "n_splits": n_splits,
        "best_params": best_params,
        "cv_rmse": cv_rmse,
        "train_rmse": train_rmse,
        "train_r2": train_r2,
        "test_rmse": test_rmse,
        "test_r2": test_r2,
        "versions": {"python": platform.python_version(),
                     "numpy": np.__version__,
                     "pandas": pd.__version__,
                     "scikit-learn": sklearn.__version__,
                     "xgboost": xgb.__version__,
                     "shap": shap.__version__},
    }
    with open(output_dir / "xgb_model_summary.json", "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2)

    # 6.5 Check that the saved model reproduces the predictions exactly.
    reloaded = xgb.XGBRegressor()
    reloaded.load_model(model_file)
    if not np.array_equal(reloaded.predict(x_all), best_model.predict(x_all)):
        raise RuntimeError("Saved model does not reproduce the predictions")
    print(f"Model saved to {model_file} and verified")

    # ------------------------------------------------------------------
    # Step 7: Calculate SHAP values of the final model for all samples.
    # The model evaluated in Step 5 is explained on the full data set
    # (training + test samples), so every sample receives an attribution.
    # ------------------------------------------------------------------
    explainer = shap.TreeExplainer(best_model)
    shap_values = explainer.shap_values(x_all)

    # ------------------------------------------------------------------
    # Step 8: Save the SHAP values and the SHAP feature importance.
    # ------------------------------------------------------------------
    # 8.1 SHAP value of every feature for every sample (rows in the same
    #     order as the input file).
    shap_values_df = pd.DataFrame(shap_values, columns=x_columns)
    shap_values_df.to_csv(output_dir / "shap_values.csv", index=False)

    # 8.2 Feature importance = mean absolute SHAP value over all samples.
    shap_importance = np.abs(shap_values).mean(axis=0)
    shap_importance_df = pd.DataFrame({"features": x_columns,
                                       "shap_impt": shap_importance})
    shap_importance_df = shap_importance_df.sort_values(by="shap_impt")
    shap_importance_df.to_csv(output_dir / "shap_importance.csv", index=False)

    # 8.3 Show the ten most important features.
    print(shap_importance_df.tail(10).iloc[::-1])
