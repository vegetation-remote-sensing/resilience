function calculate_resilience_trends(config)
%% CALCULATE_RESILIENCE_TRENDS
% Computes temporal trends in vegetation resilience
%
% Description:
%   This function calculates the Sen's slope and the Mann-Kendall
%   significance of the AR(1) resilience indicator (temporal
%   autocorrelation, TAC). For every target year t of the analysis period,
%   the trend is calculated from all annual values between baseline_start
%   and t (a window that grows by one year at each step).
%
% Workflow:
%   Step 1. Validate the configuration and set defaults.
%   Step 2. Define the time periods.
%   Step 3. Load the reference grid and identify the valid pixels.
%   Step 4. Create the output directories.
%   Step 5. Load the annual AR(1) maps.
%   Step 6. For every target year:
%           6.1 select the data from baseline_start to the target year;
%           6.2 calculate Sen's slope and the Mann-Kendall Z of every pixel;
%           6.3 classify the trend by direction and significance;
%           6.4 save the slope, Z and classification maps.
%
% Input:
%   config - Structure containing configuration parameters:
%       .input_dir                : Directory with TAC (AR1 coefficient) data
%       .output_dir               : Directory for output results
%       .reference_file           : Path to reference GeoTIFF for spatial info
%                                   (land mask; value 0 = valid pixel)
%       .baseline_start           : First year used for the trends (e.g., 2000)
%       .analysis_start           : First target year (e.g., 2009)
%       .analysis_end             : Last target year (e.g., 2019)
%       .vegetation_index         : Index name (default: 'kNDVI')
%       .resilience_indicator     : Indicator name (default: 'AR1')
%       .temporal_window          : Rolling window parameter (default: 60)
%       .decomposition_method     : Decomposition method (default: 'STL')
%       .significance_level       : Alpha for MK test (default: 0.05)
%
%   Input file naming: [index]_[indicator]_[window]_[method]_[year].tif
%
% Output (one file per target year, named ..._[baseline_start]_[year].tif):
%   - Sen_AR1/   : Sen's slope (trend magnitude, per year)
%   - MK_AR1/    : Mann-Kendall Z-statistic
%   - SenMK_AR1/ : Combined classification
%                  (-2 significant decrease, -1 non-significant decrease,
%                    0 no trend, 1 non-significant increase,
%                    2 significant increase)
%
% Requirements: MATLAB R2020a or later, Mapping Toolbox,
%               Image Processing Toolbox,
%               Statistics and Machine Learning Toolbox,
%               Parallel Computing Toolbox (optional)
%
% Author: Yiling Cai
% Institution: School of Geography and Remote Sensing, Guangzhou University, Guangzhou 510006, China;  Guangdong Guodi Science Technology Co., Ltd, Guangzhou 510075, China;
% Date: 2026-01-16
% Version: 1.1

    %% Step 1: Validate the configuration and set defaults
    config = validateTrendConfig(config);

    %% Step 2: Define the time periods
    years = config.baseline_start:config.analysis_end;           % data years
    analysis_years = config.analysis_start:config.analysis_end;  % target years

    fprintf('Calculating resilience trends from %d to %d\n', ...
            config.baseline_start, config.analysis_end);

    %% Step 3: Load the reference grid and identify the valid pixels
    [ref_data, spatial_ref, geo_info] = loadReferenceData(config.reference_file);
    [n_rows, n_cols] = size(ref_data);
    valid_pixels = find(ref_data(:) == 0);  % 0 = valid pixel

    %% Step 4: Create the output directories
    dir_sen = fullfile(config.output_dir, 'Sen_AR1');
    dir_mk = fullfile(config.output_dir, 'MK_AR1');
    dir_senmk = fullfile(config.output_dir, 'SenMK_AR1');
    createOutputDirs({dir_sen, dir_mk, dir_senmk});

    %% Step 5: Load the annual AR(1) maps
    % Result: matrix of size [n_pixels x n_years]
    fprintf('\n=== Loading AR(1) data ===\n');
    ar1_data = loadAR1Data(config, years, n_rows, n_cols);

    %% Step 6: Calculate the trends for every target year
    fprintf('\n=== Computing trends ===\n');
    calculateRollingTrends(ar1_data, years, analysis_years, ...
                          valid_pixels, n_rows, n_cols, ...
                          config, spatial_ref, geo_info, ...
                          dir_sen, dir_mk, dir_senmk);

    fprintf('\nTrend analysis complete!\n');
end

%% Helper Functions

function config = validateTrendConfig(config)
    % Check the required fields and fill in the default values

    required = {'input_dir', 'output_dir', 'reference_file', ...
                'baseline_start', 'analysis_start', 'analysis_end'};

    for i = 1:length(required)
        if ~isfield(config, required{i})
            error('Missing required field: %s', required{i});
        end
    end

    % Set defaults
    if ~isfield(config, 'vegetation_index')
        config.vegetation_index = 'kNDVI';
    end
    if ~isfield(config, 'resilience_indicator')
        config.resilience_indicator = 'AR1';
    end
    if ~isfield(config, 'decomposition_method')
        config.decomposition_method = 'STL';
    end
    if ~isfield(config, 'temporal_window')
        config.temporal_window = 60;
    end
    if ~isfield(config, 'significance_level')
        config.significance_level = 0.05;
    end

    % The periods must be in chronological order
    if config.analysis_start < config.baseline_start || ...
       config.analysis_end < config.analysis_start
        error('Periods must satisfy baseline_start <= analysis_start <= analysis_end.');
    end
end

function ar1_matrix = loadAR1Data(config, years, n_rows, n_cols)
    % Load the annual AR(1) maps into a [n_pixels x n_years] matrix
    n_years = length(years);
    ar1_matrix = NaN(n_rows * n_cols, n_years);

    for i = 1:n_years
        year = years(i);
        filename = sprintf('%s_%s_%d_%s_%04d.tif', ...
                          config.vegetation_index, ...
                          config.resilience_indicator, ...
                          config.temporal_window, ...
                          config.decomposition_method, ...
                          year);
        filepath = fullfile(config.input_dir, filename);

        % Every year is required: a missing year would bias the trends
        if ~exist(filepath, 'file')
            error('File not found: %s', filepath);
        end

        temp_data = imread(filepath);
        if ~isequal(size(temp_data), [n_rows, n_cols])
            error('Grid size of %s does not match the reference file.', filename);
        end
        ar1_matrix(:, i) = temp_data(:);

        if mod(i, 5) == 0
            fprintf('  Loaded %d/%d files\n', i, n_years);
        end
    end
end

function calculateRollingTrends(data, all_years, analysis_years, ...
                                valid_pixels, n_rows, n_cols, ...
                                config, spatial_ref, geo_info, ...
                                dir_sen, dir_mk, dir_senmk)
    % Calculate Sen-MK trends for each target year using data up to that year

    n_valid = length(valid_pixels);

    for target_year = analysis_years
        % Step 6.1: Select the data from baseline_start to the target year
        year_idx = find(all_years == target_year);
        subset_data = data(valid_pixels, 1:year_idx);  % valid pixels only
        n_timepoints = year_idx;

        fprintf('Processing year %d (n=%d)...\n', target_year, n_timepoints);

        % Step 6.2: Sen's slope and Mann-Kendall Z of every valid pixel
        % (runs in parallel if the Parallel Computing Toolbox is available)
        slopes = NaN(n_valid, 1);
        zstats = NaN(n_valid, 1);

        parfor i = 1:n_valid
            [slopes(i), zstats(i)] = senMannKendall(subset_data(i, :));
        end

        % Map the results back to the 2D grid
        slope_map = NaN(n_rows, n_cols);
        zstat_map = NaN(n_rows, n_cols);
        slope_map(valid_pixels) = slopes;
        zstat_map(valid_pixels) = zstats;

        % Step 6.3: Classify the trends by direction and significance
        combined_map = classifyTrends(slope_map, zstat_map, config.significance_level);

        % Step 6.4: Save the results
        saveTrendMaps(slope_map, zstat_map, combined_map, ...
                     target_year, config, spatial_ref, geo_info, ...
                     dir_sen, dir_mk, dir_senmk);
    end
end

function [slope, z_stat] = senMannKendall(time_series)
    % Calculate Sen's slope and the Mann-Kendall Z-statistic of one series
    %
    % Sen's slope  : median of the slopes between all pairs of time points
    %                (Sen, 1968), in units per time step (year)
    % Mann-Kendall : Z-statistic of the test for a monotonic trend
    %                (Mann, 1945; Kendall, 1975), with continuity correction
    %                and variance correction for tied values
    %
    % Missing values (NaN) are skipped; the remaining values keep their
    % original time positions. At least 3 valid values are required.

    % Time positions and values of the valid observations
    time_series = time_series(:);
    t = find(~isnan(time_series));
    x = time_series(t);
    n = length(x);

    if n < 3
        slope = NaN;
        z_stat = NaN;
        return;
    end

    % Differences between all pairs (later minus earlier observation)
    delta_value = x' - x;                 % delta_value(j, i) = x(i) - x(j)
    delta_time = t' - t;                  % delta_time(j, i)  = t(i) - t(j)
    pairs = triu(true(n), 1);             % each pair once (j < i)

    % Sen's slope = median of all pairwise slopes
    slope = median(delta_value(pairs) ./ delta_time(pairs));

    % Mann-Kendall S = number of increases minus number of decreases
    sign_sum = sum(sign(delta_value(pairs)));

    % Variance of S, corrected for groups of tied values
    [~, ~, group] = unique(x);
    ties = accumarray(group, 1);
    var_s = (n * (n - 1) * (2*n + 5) - ...
             sum(ties .* (ties - 1) .* (2*ties + 5))) / 18;

    % Z-statistic with continuity correction
    if sign_sum == 0 || var_s == 0
        z_stat = 0;
    elseif sign_sum > 0
        z_stat = (sign_sum - 1) / sqrt(var_s);
    else
        z_stat = (sign_sum + 1) / sqrt(var_s);
    end
end

function combined = classifyTrends(slope_map, zstat_map, alpha)
    % Classify trends by direction and significance
    %
    % Output codes:
    %   -2: Significant decreasing trend (p < alpha)
    %   -1: Non-significant decreasing trend
    %    0: No trend
    %    1: Non-significant increasing trend
    %    2: Significant increasing trend (p < alpha)

    % Critical Z value for two-tailed test
    z_crit = norminv(1 - alpha/2);  % e.g., 1.96 for alpha=0.05

    % Direction (sign of Sen's slope)
    trend_dir = NaN(size(slope_map));
    trend_dir(slope_map == 0) = 0;
    trend_dir(slope_map > 0) = 1;
    trend_dir(slope_map < 0) = -1;

    % Significance (Mann-Kendall test)
    abs_z = abs(zstat_map);
    trend_sig = ones(size(abs_z));  % Default: non-significant
    trend_sig(abs_z > z_crit) = 2;  % Significant

    % Combine
    combined = trend_dir .* trend_sig;
end

function saveTrendMaps(slope, zstat, combined, year, config, ...
                       spatial_ref, geo_info, dir_sen, dir_mk, dir_senmk)
    % Save the trend maps of one target year

    year_range = sprintf('%04d_%04d', config.baseline_start, year);

    % Sen's slope
    filename = sprintf('Sen_%s_%s_%d_%s_%s.tif', ...
                      config.vegetation_index, ...
                      config.resilience_indicator, ...
                      config.temporal_window, ...
                      config.decomposition_method, ...
                      year_range);
    filepath = fullfile(dir_sen, filename);
    geotiffwrite(filepath, slope, spatial_ref, ...
                'GeoKeyDirectoryTag', geo_info.GeoTIFFTags.GeoKeyDirectoryTag);

    % Mann-Kendall Z
    filename = sprintf('MK_%s_%s_%d_%s_%s.tif', ...
                      config.vegetation_index, ...
                      config.resilience_indicator, ...
                      config.temporal_window, ...
                      config.decomposition_method, ...
                      year_range);
    filepath = fullfile(dir_mk, filename);
    geotiffwrite(filepath, zstat, spatial_ref, ...
                'GeoKeyDirectoryTag', geo_info.GeoTIFFTags.GeoKeyDirectoryTag);

    % Combined classification
    filename = sprintf('SenMK_%s_%s_%d_%s_%s.tif', ...
                      config.vegetation_index, ...
                      config.resilience_indicator, ...
                      config.temporal_window, ...
                      config.decomposition_method, ...
                      year_range);
    filepath = fullfile(dir_senmk, filename);
    geotiffwrite(filepath, combined, spatial_ref, ...
                'GeoKeyDirectoryTag', geo_info.GeoTIFFTags.GeoKeyDirectoryTag);
end

function [data, spatial_ref, geo_info] = loadReferenceData(filepath)
    % Load the reference GeoTIFF (grid definition and land mask)
    if ~exist(filepath, 'file')
        error('Reference file not found: %s', filepath);
    end

    [data, spatial_ref] = geotiffread(filepath);
    geo_info = geotiffinfo(filepath);
end

function createOutputDirs(dir_list)
    % Create output directories if they don't exist
    for i = 1:length(dir_list)
        if ~exist(dir_list{i}, 'dir')
            mkdir(dir_list{i});
            fprintf('Created directory: %s\n', dir_list{i});
        end
    end
end
