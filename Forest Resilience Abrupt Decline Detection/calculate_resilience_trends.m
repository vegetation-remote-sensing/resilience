function calculate_resilience_trends(config)
%% CALCULATE_RESILIENCE_TRENDS 
% Computes temporal trends in vegetation resilience
%
% Description:
%   This function calculates Sen's slope and Mann-Kendall significance for
%   AR(1) resilience indicators over rolling time windows. 
%
% Input:
%   config - Structure containing configuration parameters:
%       .input_dir                : Directory with TAC (AR1 coefficient) data
%       .output_dir               : Directory for output results
%       .reference_file           : Path to reference GeoTIFF for spatial info
%       .baseline_start           : Start of baseline period (e.g., 2000)
%       .baseline_end             : End of baseline period (e.g., 2009)
%       .analysis_end             : End of analysis period (e.g., 2022)
%       .vegetation_index         : Index name (e.g., 'kNDVI')
%       .resilience_indicator     : Indicator name (e.g., 'AR1')
%       .temporal_window          : Rolling window parameter (e.g., 60)
%       .decomposition_method     : Decomposition method (e.g., 'STL')
%       .significance_level       : Alpha for MK test (default: 0.05)
%
% Output:
%   - Sen's slope GeoTIFF files (trend magnitude)
%   - Mann-Kendall Z-statistic GeoTIFF files
%   - Combined Sen-MK classification GeoTIFF files
%
% Author: Yiling Cai
% Institution: School of Geography and Remote Sensing, Guangzhou University, Guangzhou 510006, China;  Guangdong Guodi Science Technology Co., Ltd, Guangzhou 510075, China;
% Date: 2026-01-16
% Version: 1.0

    %% Validate configuration
    validateTrendConfig(config);
    
    %% Initialize parameters
    years = config.baseline_start:config.analysis_end;
    analysis_years = config.analysis_start:config.analysis_end;
    
    fprintf('Calculating resilience trends from %d to %d\n', ...
            config.baseline_start, config.analysis_end);
    
    %% Load reference data
    [ref_data, spatial_ref, geo_info] = loadReferenceData(config.reference_file);
    [n_rows, n_cols] = size(ref_data);
    valid_pixels = find(ref_data(:) == 0);
    
    %% Create output directories
    dir_sen = fullfile(config.output_dir, 'Sen_AR1');
    dir_mk = fullfile(config.output_dir, 'MK_AR1');
    dir_senmk = fullfile(config.output_dir, 'SenMK_AR1');
    createOutputDirs({dir_sen, dir_mk, dir_senmk});
    
    %% Load AR(1) time series
    fprintf('\n=== Loading AR(1) data ===\n');
    ar1_data = loadAR1Data(config, years, n_rows * n_cols);
    
    %% Calculate trends for rolling windows
    fprintf('\n=== Computing trends ===\n');
    calculateRollingTrends(ar1_data, years, analysis_years, ...
                          valid_pixels, n_rows, n_cols, ...
                          config, spatial_ref, geo_info, ...
                          dir_sen, dir_mk, dir_senmk);
    
    fprintf('\nTrend analysis complete!\n');
end

%% Helper Functions

function validateTrendConfig(config)
    % Validate required configuration fields
    required = {'input_dir', 'output_dir', 'reference_file', ...
                'baseline_start', 'analysis_end', 'baseline_end'};
    
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
end

function ar1_matrix = loadAR1Data(config, years, n_pixels)
    % Load AR(1) coefficient time series
    n_years = length(years);
    ar1_matrix = NaN(n_pixels, n_years);
    
    for i = 1:n_years
        year = years(i);
        filename = sprintf('%s_%s_%d_%s_%04d.tif', ...
                          config.vegetation_index, ...
                          config.resilience_indicator, ...
                          config.temporal_window, ...
                          config.decomposition_method, ...
                          year);
        filepath = fullfile(config.input_dir, filename);
        
        if ~exist(filepath, 'file')
            warning('File not found: %s', filepath);
            continue;
        end
        
        temp_data = imread(filepath);
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
    % Calculate Sen-MK trends for each year using data up to that year
    
    n_valid = length(valid_pixels);
    
    for target_year = analysis_years
        year_idx = find(all_years == target_year);
        subset_data = data(:, 1:year_idx);
        n_timepoints = year_idx;
        
        fprintf('Processing year %d (n=%d)...\n', target_year, n_timepoints);
        
        % Initialize output arrays
        slope_map = NaN(n_rows, n_cols);
        zstat_map = NaN(n_rows, n_cols);
        
        % Parallel processing over valid pixels
        slopes = NaN(n_valid, 1);
        zstats = NaN(n_valid, 1);
        
        parfor i = 1:n_valid
            pixel_idx = valid_pixels(i);
            time_series = subset_data(pixel_idx, :);
            
            if sum(~isnan(time_series)) < 3
                continue;  % Need at least 3 points
            end
            
            [slopes(i), zstats(i)] = senMannKendall(time_series);
        end
        
        % Map results back to 2D
        slope_map(valid_pixels) = slopes;
        zstat_map(valid_pixels) = zstats;
        
        % Classify trends
        combined_map = classifyTrends(slope_map, zstat_map, config.significance_level);
        
        % Save results
        saveTrendMaps(slope_map, zstat_map, combined_map, ...
                     target_year, config, spatial_ref, geo_info, ...
                     dir_sen, dir_mk, dir_senmk);
    end
end

function [slope, z_stat] = senMannKendall(time_series)
    % Calculate Sen's slope and Mann-Kendall Z-statistic
    %
    % Sen's slope: median of all pairwise slopes
    % Mann-Kendall: tests for monotonic trend
    
    % Remove NaN values
    valid_data = time_series(~isnan(time_series));
    n = length(valid_data);
    
    if n < 3
        slope = NaN;
        z_stat = NaN;
        return;
    end
    
    % Calculate all pairwise slopes
    slopes = [];
    sign_sum = 0;
    
    for i = 2:n
        for j = 1:(i-1)
            delta_value = valid_data(i) - valid_data(j);
            delta_time = i - j;
            slopes = [slopes; delta_value / delta_time];
            
            % Mann-Kendall sign
            if delta_value > 0
                sign_sum = sign_sum + 1;
            elseif delta_value < 0
                sign_sum = sign_sum - 1;
            end
        end
    end
    
    % Sen's slope = median of all slopes
    slope = median(slopes);
    
    % Mann-Kendall Z-statistic
    var_s = n * (n - 1) * (2*n + 5) / 18;
    
    if sign_sum == 0
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
    
    % Direction
    trend_dir = NaN(size(slope_map));
    trend_dir(slope_map == 0) = 0;
    trend_dir(slope_map > 0) = 1;
    trend_dir(slope_map < 0) = -1;
    
    % Significance
    abs_z = abs(zstat_map);
    trend_sig = ones(size(abs_z));  % Default: non-significant
    trend_sig(abs_z > z_crit) = 2;  % Significant
    
    % Combine
    combined = trend_dir .* trend_sig;
end

function saveTrendMaps(slope, zstat, combined, year, config, ...
                       spatial_ref, geo_info, dir_sen, dir_mk, dir_senmk)
    % Save trend analysis results
    
    year_range = sprintf('%04d_%04d', config.baseline_start, year);
    
    % Sen's slope
    filename = sprintf('Sen_%s_%s_%d_%s_%s.tif', ...
                      config.vegetation_index, ...
                      config.resilience_indicator, ...
                      config.temporal_window, ...
                      config.decomposition_method, ...
                      config.baseline_start, ...
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
                      config.baseline_start, ...
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
                      config.baseline_start, ...
                      year_range);
    filepath = fullfile(dir_senmk, filename);
    geotiffwrite(filepath, combined, spatial_ref, ...
                'GeoKeyDirectoryTag', geo_info.GeoTIFFTags.GeoKeyDirectoryTag);
end

function [data, spatial_ref, geo_info] = loadReferenceData(filepath)
    % Load reference GeoTIFF file
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

