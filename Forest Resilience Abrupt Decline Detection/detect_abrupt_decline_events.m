function detect_abrupt_decline_events(config)
% DETECT_ABRUPT_DECLINE_EVENTS Identifies abrupt decline events in vegetation
%
% Description:
%   This function detects abrupt decline (AD) events by calculating anomalies
%   in growing season kNDVI (kernel Normalized Difference Vegetation Index).
%   AD events are defined as periods when kNDVI falls below the baseline mean
%   by more than n standard deviations (kNDVI < μ - nσ).
%
% Input:
%   config - Structure containing configuration parameters:
%       .input_dir          : Directory containing input kNDVI data
%       .output_dir         : Directory for output results
%       .reference_file     : Path to reference GeoTIFF for spatial info
%       .baseline_start     : Start year of baseline period (e.g., 2000)
%       .baseline_end       : End year of baseline period (e.g., 2009)
%       .analysis_start     : Start year of analysis period (e.g., 2010)
%       .analysis_end       : End year of analysis period (e.g., 2022)
%       .growing_season     : String describing season (e.g., 'GS_05_09')
%       .sigma_values       : AD threshold n (e.g., 1:6)
%       .sigma_denominator  : denominator of noAD threshold m (e.g., 2:6)
%       .vegetation_index   : Index name (e.g., 'kNDVI')
%
% Output:
%   - Detrended kNDVI GeoTIFF files
%   - Sigma threshold GeoTIFF files
%   - Binary AD/no-AD classification GeoTIFF files
%
% Author: Yiling Cai
% Institution: School of Geography and Remote Sensing, Guangzhou University, Guangzhou 510006, China;  Guangdong Guodi Science Technology Co., Ltd, Guangzhou 510075, China;
% Date: 2026-01-16
% Version: 1.0

    %% Validate input configuration
    validateADConfig(config);
    
    %% Initialize parameters
    baseline_years = config.baseline_start:config.baseline_end;
    analysis_years = config.analysis_start:config.analysis_end;
    all_years = config.baseline_start:config.analysis_end;
    
    fprintf('Processing AD detection from %d to %d\n', ...
            config.analysis_start, config.analysis_end);
    fprintf('Baseline period: %d-%d\n', ...
            config.baseline_start, config.baseline_end);
    
    %% Load reference spatial information
    [ref_data, spatial_ref, geo_info] = loadReferenceData(config.reference_file);
    [n_rows, n_cols] = size(ref_data);
    n_pixels = n_rows * n_cols;
    
    % Identify valid pixels (assuming 0 indicates valid data)
    valid_pixels = find(ref_data(:) == 0);
    fprintf('Valid pixels: %d (%.2f%%)\n', ...
            length(valid_pixels), 100*length(valid_pixels)/n_pixels);
    
    %% Create output directories
    dir_detrend = fullfile(config.output_dir, 'detrended_kNDVI_growing_season');
    dir_sigma = fullfile(config.output_dir, 'sigma_thresholds');
    dir_ad = fullfile(config.output_dir, 'AD_noAD_classification');
    createOutputDirs({dir_detrend, dir_sigma, dir_ad});
    
    %% Load and detrend kNDVI time series
    fprintf('\n=== Loading kNDVI data ===\n');
    kndvi_data = loadTimeSeriesData(config, all_years, n_pixels);
    
    fprintf('\n=== Detrending time series ===\n');
    detrended_data = detrendTimeSeries(kndvi_data, valid_pixels);
    
    % Save detrended data
    saveDetrended(detrended_data, all_years, dir_detrend, ...
                  config, n_rows, n_cols, spatial_ref, geo_info);
    
    %% Calculate baseline statistics
    fprintf('\n=== Computing baseline statistics ===\n');
    baseline_idx = 1:length(baseline_years);
    baseline_data = detrended_data(baseline_idx, :);
    
    mu = mean(baseline_data, 1, 'omitnan');
    sigma = std(baseline_data, 1, 1, 'omitnan');
    
    %% Detect AD events using multiple thresholds
    fprintf('\n=== Detecting AD events ===\n');
    detectADEvents(detrended_data, mu, sigma, analysis_years, ...
                   config, n_rows, n_cols, ...
                   spatial_ref, geo_info, dir_sigma, dir_ad);
    
    fprintf('\nProcessing complete!\n');
end

%% Helper Functions

function validateADConfig(config)
    % Validate that all required fields are present
    required_fields = {'input_dir', 'output_dir', 'reference_file', ...
                      'baseline_start', 'baseline_end', ...
                      'analysis_start', 'analysis_end', 'growing_season'};
    
    for i = 1:length(required_fields)
        if ~isfield(config, required_fields{i})
            error('Missing required field: %s', required_fields{i});
        end
    end
    
    % Set defaults
    if ~isfield(config, 'sigma_values')
        config.sigma_values = 1:6;
    end
    if ~isfield(config, 'sigma_denominator')
        config.sigma_denominator = 2:6;
    end
    if ~isfield(config, 'vegetation_index')
        config.vegetation_index = 'kNDVI';
    end
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

function data_matrix = loadTimeSeriesData(config, years, n_pixels)
    % Load time series of kNDVI data
    n_years = length(years);
    data_matrix = NaN(n_years, n_pixels);
    
    for i = 1:n_years
        year = years(i);
        filename = sprintf('%s_%s_%04d.tif', ...
                          config.vegetation_index,config.growing_season, year);
        filepath = fullfile(config.input_dir, filename);
        
        if ~exist(filepath, 'file')
            warning('File not found: %s', filepath);
            continue;
        end
        
        temp_data = imread(filepath);
        data_matrix(i, :) = temp_data(:);
        
        if mod(i, 5) == 0
            fprintf('  Loaded %d/%d files\n', i, n_years);
        end
    end
    
    fprintf('  Loaded %d/%d files total\n', n_years, n_years);
end

function detrended = detrendTimeSeries(data, valid_pixels)
    % Remove linear trends from time series
    [n_years, n_pixels] = size(data);
    detrended = NaN(n_years, n_pixels);
    
    % Detrend only valid pixels
    detrended(:, valid_pixels) = detrend(data(:, valid_pixels));
    
    fprintf('  Detrended %d pixels\n', length(valid_pixels));
end

function saveDetrended(data, years, output_dir, config, ...
                       n_rows, n_cols, spatial_ref, geo_info)
    % Save detrended data as GeoTIFF files
    for i = 1:length(years)
        year = years(i);
        filename = sprintf('%s_%s_detrended_%04d.tif', ...
                          config.vegetation_index,config.growing_season, year);
        filepath = fullfile(output_dir, filename);
        
        data_2d = reshape(data(i, :), n_rows, n_cols);
        geotiffwrite(filepath, data_2d, spatial_ref, ...
                    'GeoKeyDirectoryTag', geo_info.GeoTIFFTags.GeoKeyDirectoryTag);
    end
    
    fprintf('  Saved %d detrended files\n', length(years));
end

function detectADEvents(data, mu, sigma, years, config, ...
                        n_rows, n_cols, ...
                        spatial_ref, geo_info, dir_sigma, dir_ad)
    % Detect AD events using multiple sigma thresholds
    
    year_offset = config.analysis_start - config.baseline_start;
    
    for n = config.sigma_values
        for m = config.sigma_denominator
            % Calculate thresholds
            threshold_ad = mu - sigma * n;
            threshold_noad = mu - sigma * n / m;
            
            % Save threshold maps
            saveThreshold(threshold_ad, n, 0, dir_sigma, config, ...
                         n_rows, n_cols, spatial_ref, geo_info);
            saveThreshold(threshold_noad, n, m, dir_sigma, config, ...
                         n_rows, n_cols, spatial_ref, geo_info);
            
            % Classify each year
            for i = 1:length(years)
                year = years(i);
                year_idx = year_offset + i;
                
                ad_map = classifyAD(data(year_idx, :), ...
                                   threshold_ad, threshold_noad, ...
                                   n_rows, n_cols);
                
                saveADMap(ad_map, year, n, m, dir_ad, ...
                         spatial_ref, geo_info);
            end
        end
        
        fprintf('  Completed sigma=%d threshold\n', n);
    end
end

function ad_map = classifyAD(data_vector, thresh_ad, thresh_noad, n_rows, n_cols)
    % Classify pixels as AD (1) or no-AD (0)
    ad_vector = NaN(size(data_vector));
    
    % AD event: below lower threshold
    ad_vector(data_vector < thresh_ad) = 1;
    
    % No AD: above higher threshold
    ad_vector(data_vector > thresh_noad) = 0;
    
    % Reshape to 2D
    ad_map = reshape(ad_vector, n_rows, n_cols);
end

function saveThreshold(threshold, n, m, output_dir, config, ...
                       n_rows, n_cols, spatial_ref, geo_info)
    % Save threshold map
    if m == 0
        filename = sprintf('threshold_sigma%d_%04d_%04d.tif', ...
                          n, ...
                          config.baseline_start, config.baseline_end);
    else
        filename = sprintf('threshold_sigma%ddiv%d_%04d_%04d.tif', ...
                          n, m, ...
                          config.baseline_start, config.baseline_end);
    end
    
    filepath = fullfile(output_dir, filename);
    data_2d = reshape(threshold, n_rows, n_cols);
    
    geotiffwrite(filepath, data_2d, spatial_ref, ...
                'GeoKeyDirectoryTag', geo_info.GeoTIFFTags.GeoKeyDirectoryTag);
end

function saveADMap(ad_map, year, n, m, output_dir, spatial_ref, geo_info)
    % Save AD classification map
    filename = sprintf('AD_noAD_classification_sigma%ddiv%d_%04d.tif', n, m, year);
    filepath = fullfile(output_dir, filename);
    
    geotiffwrite(filepath, ad_map, spatial_ref, ...
                'GeoKeyDirectoryTag', geo_info.GeoTIFFTags.GeoKeyDirectoryTag);
end
