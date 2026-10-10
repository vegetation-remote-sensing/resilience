%% Example usage script for abrupt decline detection, resilience trend analysis, and paired AD-noAD resilience extraction.
% This script demonstrates how to use the standardized functions
% for analyzing vegetation resilience and detecting abrupt declines.
%
% Workflow:
%   Step 1. Detect abrupt decline (AD) events
%           (detect_abrupt_decline_events).
%   Step 2. Calculate resilience trends (deltaTAC)
%           (calculate_resilience_trends).
%   Step 3. Extract paired AD-noAD resilience metrics
%           (extract_paired_AD_noAD_resilience).
%
% The three steps must be run in this order: Step 3 reads the AD/noAD maps
% of Step 1 and the Sen's slope maps of Step 2.
%
% Demo data (provided in ./input/, 26 x 26 pixel subset; see Readme.md):
%   ./input/kNDVI_GS/kNDVI_GS_05_09_YYYY.tif     growing-season kNDVI
%   ./input/kNDVI_AR1/kNDVI_AR1_60_STL_YYYY.tif  annual AR(1) maps
%   ./input/OL/LC_MASK.tif                       land mask (0 = valid pixel)
%
% Run this script from the folder that contains it. The outputs are written
% to ./results/; the expected outputs are provided in that folder.


clear; clc;

%% ====================================================================
%% Step 1: Detect abrupt decline (AD) events
%% ====================================================================

%% Step 1.1: Configuration for AD detection
config_ad = struct();

% === Input directories ===
% growing season kNDVI maps
config_ad.input_dir = './input/kNDVI_GS/';

% Spatial reference (mask value = 0)
config_ad.reference_file = './input/OL/LC_MASK.tif';

% === Output directory ===
config_ad.output_dir = './results/';

% === Time periods ===
config_ad.baseline_start = 2000;   % Start of baseline period
config_ad.baseline_end = 2009;     % End of baseline period
config_ad.analysis_start = 2010;   % Start of analysis period
config_ad.analysis_end = 2020;     % End of analysis period

% === Data identifiers ===
config_ad.vegetation_index = 'kNDVI';
config_ad.growing_season = 'GS_05_09';  % May-September

% === Analysis parameters ===
config_ad.sigma_values = 1:6;              % AD threshold severity levels (n)
config_ad.sigma_denominator = 2:6;          % denominator of noAD threshold (m)

%% Step 1.2: Run AD detection
tic;
detect_abrupt_decline_events(config_ad);
elapsed = toc;

fprintf('\n=== Analysis Complete ===\n');
fprintf('AD results: %s\n', config_ad.output_dir);
fprintf('\nTotal processing time: %.1f minutes\n', elapsed/60);

%% ====================================================================
%% Step 2: Calculate resilience trends (deltaTAC)
%% ====================================================================

%% Step 2.1: Configuration for trend analysis
config_trend = struct();

% === Input directories ===
% Resilience maps (TAC = AR1 coefficient)
config_trend.input_dir = './input/kNDVI_AR1/';

% === Output directory ===
config_trend.output_dir = './results/';

% Spatial reference (mask value = 0)
config_trend.reference_file = './input/OL/LC_MASK.tif';

% === Time periods ===
% Resilience before an AD event in year t is measured up to year t-1, so
% the target years of the trend analysis are shifted back by one year
% relative to the AD analysis period (2010-2020 -> 2009-2019).
config_trend.baseline_start = 2000;       % First year used for the trends
config_trend.analysis_start = 2010 - 1;   % First target year
config_trend.analysis_end = 2020 - 1;     % Last target year

% === Data identifiers ===
config_trend.vegetation_index = 'kNDVI';
config_trend.resilience_indicator = 'AR1';
config_trend.decomposition_method = 'STL';
config_trend.temporal_window = 60;

% === Analysis parameters ===
config_trend.significance_level = 0.05;   % alpha of the Mann-Kendall test

%% Step 2.2: Run trend analysis
tic;
calculate_resilience_trends(config_trend);
elapsed = toc;

fprintf('\n=== Analysis Complete ===\n');
fprintf('Trend results: %s\n', config_trend.output_dir);
fprintf('\nTotal processing time: %.1f minutes\n', elapsed/60);

%% ====================================================================
%% Step 3: Extract paired AD-noAD resilience metrics
%% ====================================================================

%% Step 3.1: Configuration for paired extraction
config_paired = struct();

% === Input directories ===
% Resilience maps (TAC = AR1 coefficient)
config_paired.resilience_dir = './input/kNDVI_AR1/';

% Spatial reference (mask value = 0)
config_paired.reference_file = './input/OL/LC_MASK.tif';

% AD/noAD classification maps (output of Step 1)
config_paired.ad_classification_dir = './results/AD_noAD_classification/';

% Resilience trend maps (δTAC = Sen's slope; output of Step 2)
config_paired.trend_dir = './results/Sen_AR1/';


% === Output directory ===
config_paired.output_dir = './results/';

% === Time periods (same as Step 1) ===
config_paired.baseline_start = 2000;   % Start of baseline period
config_paired.analysis_start = 2010;   % Start of analysis period
config_paired.analysis_end = 2020;     % End of analysis period

% === Analysis parameters ===
config_paired.sigma_values = 1:6;              % AD threshold severity levels (n)
config_paired.sigma_denominator = 2:6;          % denominator of noAD threshold (m)
config_paired.window_sizes = [5, 7, 9];        % Neighborhood sizes (5×5, 7×7, 9×9)

% === Data identifiers ===
config_paired.vegetation_index = 'kNDVI';
config_paired.resilience_indicator = 'AR1';
config_paired.decomposition_method = 'STL';
config_paired.temporal_window = 60;

%% Step 3.2: Run extraction
fprintf('Starting paired AD-noAD metric extraction...\n');


tic;
extract_paired_AD_noAD_resilience(config_paired);
elapsed = toc;

fprintf('\n=== Analysis Complete ===\n');
fprintf('AD-noAD resilience results: %s\n', config_paired.output_dir);
fprintf('\nTotal processing time: %.1f minutes\n', elapsed/60);
