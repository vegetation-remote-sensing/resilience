%% Example usage script for abrupt decline detection, resilience trend analysis, and paired AD-noAD resilience extraction.
% This script demonstrates how to use the standardized functions
% for analyzing vegetation resilience and detecting abrupt declines


clear; clc;

%% Configuration for AD detection
config_ad = struct();

% === Input directories ===
% growing season kNDVI maps
config_ad.input_dir = './data/kNDVI_growing_season/';

% Spatial reference (mask value = 0)
config_ad.reference_file = './data/reference/land_mask.tif';

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
config_ad.sigma_values = 1:6;              % AD threshold severity levels
config_ad.sigma_denominator = 2:6;          % denominator of noAD threshold (m)

%% Run AD detection
tic;
detect_abrupt_decline_events(config_ad);
elapsed = toc;

fprintf('\n=== Analysis Complete ===\n');
fprintf('AD results: %s\n', config_ad.output_dir);
fprintf('\nTotal processing time: %.1f minutes\n', elapsed/60);

%% Configuration for trend analysis
config_trend = struct();

% === Input directories ===
% Resilience maps (TAC = AR1 coefficient)
config_trend.input_dir = './data/AR1/';

% === Output directory ===
config_trend.output_dir = './results/';

% Spatial reference (mask value = 0)
config_trend.reference_file = './data/reference/land_mask.tif';

% === Time periods ===
config_trend.baseline_start = 2000;   % Start of baseline period
config_trend.baseline_end = 2009;     % End of baseline period
config_trend.analysis_start = 2010 - 1;   % Start of analysis period
config_trend.analysis_end = 2020 - 1;     % End of analysis period

% === Data identifiers ===
config_trend.vegetation_index = 'kNDVI';
config_trend.resilience_indicator = 'AR1';
config_trend.decomposition_method = 'STL';
config_trend.temporal_window = 60;

% === Analysis parameters ===
config_trend.significance_level = 0.05;

%% Run trend analysis
tic;
calculate_resilience_trends(config_trend);
elapsed = toc;

fprintf('\n=== Analysis Complete ===\n');
fprintf('Trend results: %s\n', config_trend.output_dir);
fprintf('\nTotal processing time: %.1f minutes\n', elapsed/60);

%% Configuration
config_paired = struct();

% === Input directories ===
% Resilience maps (TAC = AR1 coefficient)
config_paired.resilience_dir = './data/AR1/';

% Spatial reference (mask value = 0)
config_paired.reference_file = './data/reference/land_mask.tif';

% AD/noAD classification maps (output from previous step)
config_paired.ad_classification_dir = './results/AD_noAD_classification/';

% Resilience trend maps (δTAC = Sen's slope)
config_paired.trend_dir = './results/Sen_AR1/';


% === Output directory ===
config_paired.output_dir = './results/';

% === Time periods ===
config_paired.baseline_start = 2000;   % Start of baseline period
config_paired.baseline_end = 2009;     % End of baseline period
config_paired.analysis_start = 2010;   % Start of analysis period
config_paired.analysis_end = 2020;     % End of analysis period

% === Analysis parameters ===
config_paired.sigma_values = 1:6;              % AD threshold severity levels
config_paired.sigma_denominator = 2:6;          % denominator of noAD threshold (m)
config_paired.window_sizes = [5, 7, 9];        % Neighborhood sizes (5×5, 7×7, 9×9)

% === Data identifiers ===
config_paired.vegetation_index = 'kNDVI';
config_paired.resilience_indicator = 'AR1';
config_paired.decomposition_method = 'STL';
config_paired.temporal_window = 60;

%% Run extraction
fprintf('Starting paired AD-noAD metric extraction...\n');


tic;
extract_paired_AD_noAD_resilience(config_paired);
elapsed = toc;

fprintf('\n=== Analysis Complete ===\n');
fprintf('AD-noAD resilience results: %s\n', config_paired.output_dir);
fprintf('\nTotal processing time: %.1f minutes\n', elapsed/60);





