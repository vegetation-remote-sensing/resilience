function extract_paired_AD_noAD_resilience(config)
%% EXTRACT_PAIRED_AD_NOAD_RESILIENCE Extracts 
% resilience metrics for paired AD and noAD pixels
%
% Description:
%   For each pixel that experiences an abrupt decline (AD) event, this function:
%   1. Identifies the first year of AD occurrence
%   2. Extracts δTAC (resilience trend) up to year t-1
%   3. Extracts TAC (resilience level) at year t-1
%   4. Identifies noAD control pixels in a spatial neighborhood
%   5. Calculates mean δTAC and TAC for noAD pixels
%
%   This creates matched case-control pairs for analyzing abrupt decline events.
%
% Input:
%   config - Structure containing configuration parameters:
%       .ad_classification_dir  : Directory with AD/noAD classification maps
%       .trend_dir              : Directory with δTAC (Sen's slope) data
%       .resilience_dir         : Directory with TAC (AR1 coefficient) data
%       .reference_file         : Path to reference GeoTIFF
%       .output_dir             : Base directory for outputs
%       .baseline_start         : Baseline start year (e.g., 2000)
%       .analysis_start         : Analysis start year (e.g., 2010)
%       .analysis_end           : Analysis end year (e.g., 2020)
%       .sigma_values           : AD threshold n (e.g., 1:6)
%       .sigma_denominator       : denominator of noAD threshold m (e.g., 2:6)
%       .window_sizes           : Neighborhood sizes (e.g., [5 7 9])
%       .vegetation_index       : Index name (e.g., 'kNDVI')
%       .resilience_indicator   : Indicator name (e.g., 'AR1')
%       .decomposition_method   : Method (e.g., 'STL')
%       .temporal_window        : Window parameter (e.g., 60)
%
% Output:
%   For each configuration (sigma, m, window size):
%   - AD_year maps: Year of first AD occurrence
%   - AD_δTAC maps: Resilience trend before AD
%   - AD_TAC maps: Resilience level before AD
%   - noAD_count maps: Number of noAD pixels in neighborhood
%   - noAD_δTAC maps: Mean resilience trend of noAD pixels
%   - noAD_TAC maps: Mean resilience level of noAD pixels
%
% Author: Yiling Cai
% Institution: School of Geography and Remote Sensing, Guangzhou University, Guangzhou 510006, China;  Guangdong Guodi Science Technology Co., Ltd, Guangzhou 510075, China;
% Date: 2026-01-16
% Version: 1.0

    %% Validate configuration
    validatePairingConfig(config);
    
    %% Initialize
    fprintf('=== Paired AD-noAD Resilience Extraction ===\n');
    fprintf('Period: %d-%d (baseline: %d)\n', ...
            config.analysis_start, config.analysis_end, config.baseline_start);
    
    %% Load spatial reference
    [ref_data, spatial_ref, geo_info] = loadReferenceData(config.reference_file);
    [n_rows, n_cols] = size(ref_data);
    n_years = config.analysis_end - config.analysis_start + 1;
    
    fprintf('Grid dimensions: %d × %d\n', n_rows, n_cols);
    
    %% Create output directory structure
    output_dirs = createOutputStructure(config.output_dir);
    
    %% Load time series data
    fprintf('\n=== Loading temporal data ===\n');
    
    % Load δTAC (resilience trends)
    trend_data = loadTrendTimeSeries(config, n_rows, n_cols, n_years);
    
    % Load TAC (resilience levels)
    resilience_data = loadResilienceTimeSeries(config, n_rows, n_cols, n_years);
    
    %% Process each parameter combination
    total_configs = length(config.sigma_values) * ...
                   length(config.sigma_denominator) * ...
                   length(config.window_sizes);
    config_count = 0;
    
    for sigma_n = config.sigma_values
        for sigma_m = config.sigma_denominator
            
            % Load AD/noAD classification time series
            ad_classification = loadADClassificationSeries(config, sigma_n, sigma_m, ...
                                                          n_rows, n_cols, n_years);
            
            for window_idx = 1:length(config.window_sizes)
                config_count = config_count + 1;
                window_size = config.window_sizes(window_idx);
                
                fprintf('\n[%d/%d] Processing: n=%d, m=%d, window=%d×%d\n', ...
                        config_count, total_configs, sigma_n, sigma_m, ...
                        window_size, window_size);
                
                % Extract paired metrics
                [ad_metrics, noad_metrics] = extractPairedMetrics(...
                    ad_classification, trend_data, resilience_data, ...
                    window_size, config, n_rows, n_cols);
                
                % Save results
                savePairedResults(ad_metrics, noad_metrics, ...
                                 sigma_n, sigma_m, window_size, ...
                                 config, spatial_ref, geo_info, output_dirs);
            end
        end
    end
    
    fprintf('\n=== Processing complete ===\n');
    fprintf('Results saved to: %s\n', config.output_dir);
end

%% Configuration Validation

function validatePairingConfig(config)
    % Validate required configuration fields
    
    required_fields = {
        'ad_classification_dir', 'trend_dir', 'resilience_dir', ...
        'reference_file', 'output_dir', 'baseline_start', ...
        'analysis_start', 'analysis_end'
    };
    
    for i = 1:length(required_fields)
        if ~isfield(config, required_fields{i})
            error('Missing required configuration field: %s', required_fields{i});
        end
    end
    
    % Set defaults
    if ~isfield(config, 'sigma_values')
        config.sigma_values = 1:6;
    end
    if ~isfield(config, 'sigma_denominator')
        config.sigma_denominator = 2:6;
    end
    if ~isfield(config, 'window_sizes')
        config.window_sizes = [5, 7, 9];
    end
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
    
    % Validate directories exist
    if ~exist(config.ad_classification_dir, 'dir')
        error('AD classification directory not found: %s', config.ad_classification_dir);
    end
    if ~exist(config.trend_dir, 'dir')
        error('Trend directory not found: %s', config.trend_dir);
    end
    if ~exist(config.resilience_dir, 'dir')
        error('Resilience directory not found: %s', config.resilience_dir);
    end
end

%% Directory Management

function output_dirs = createOutputStructure(base_dir)
    % Create standardized output directory structure
    
    subdirs = {
        'AD_year', 'AD_deltaTAC', 'AD_TAC', ...
        'noAD_count', 'noAD_deltaTAC', 'noAD_TAC'
    };
    
    output_dirs = struct();
    
    for i = 1:length(subdirs)
        dir_path = fullfile(base_dir, subdirs{i});
        if ~exist(dir_path, 'dir')
            mkdir(dir_path);
        end
        % Store path using valid field name
        field_name = strrep(subdirs{i}, '_', '');
        output_dirs.(field_name) = dir_path;
    end
    
    fprintf('Created output directory structure in: %s\n', base_dir);
end

%% Data Loading Functions

function trend_array = loadTrendTimeSeries(config, n_rows, n_cols, n_years)
    % Load δTAC (resilience trend) time series
    %
    % File naming: Sen_[index]_[indicator]_[window]_[method]_[year0]_[year1].tif
    
    fprintf('  Loading δTAC time series...\n');
    
    trend_array = NaN(n_rows, n_cols, n_years);
    
    for i = 1:n_years
        year_end = config.analysis_start + i - 2;  % t-1 for year t
        
        filename = sprintf('Sen_%s_%s_%d_%s_%d_%d.tif', ...
                          config.vegetation_index, ...
                          config.resilience_indicator, ...
                          config.temporal_window, ...
                          config.decomposition_method, ...
                          config.baseline_start, ...
                          year_end);
        
        filepath = fullfile(config.trend_dir, filename);
        
        if ~exist(filepath, 'file')
            warning('Trend file not found: %s', filename);
            continue;
        end
        
        trend_array(:, :, i) = imread(filepath);
    end
    
    fprintf('    Loaded %d trend maps\n', n_years);
end

function resilience_array = loadResilienceTimeSeries(config, n_rows, n_cols, n_years)
    % Load TAC (resilience indicator) time series
    %
    % File naming: [index]_[indicator]_[window]_[method]_[year].tif
    
    fprintf('  Loading TAC time series...\n');
    
    resilience_array = NaN(n_rows, n_cols, n_years);
    
    for i = 1:n_years
        year = config.analysis_start + i - 2;  % t-1 for year t
        
        filename = sprintf('%s_%s_%d_%s_%d.tif', ...
                          config.vegetation_index, ...
                          config.resilience_indicator, ...
                          config.temporal_window, ...
                          config.decomposition_method, ...
                          year);
        
        filepath = fullfile(config.resilience_dir, filename);
        
        if ~exist(filepath, 'file')
            warning('Resilience file not found: %s', filename);
            continue;
        end
        
        resilience_array(:, :, i) = imread(filepath);
    end
    
    fprintf('    Loaded %d resilience maps\n', n_years);
end

function ad_array = loadADClassificationSeries(config, sigma_n, sigma_m, ...
                                               n_rows, n_cols, n_years)
    % Load AD/noAD classification time series
    %
    % File naming: AD_noAD_classification_sigma[n]div[m]_[year].tif
    % Values: 1 = AD, 0 = noAD, NaN = invalid
    
    fprintf('  Loading AD/noAD classifications (σ=%d, m=%d)...\n', sigma_n, sigma_m);
    
    ad_array = NaN(n_rows, n_cols, n_years);
    
    for i = 1:n_years
        year = config.analysis_start + i - 1;
        
        filename = sprintf('AD_noAD_classification_sigma%ddiv%d_%d.tif', ...
                          sigma_n, sigma_m, year);
        
        filepath = fullfile(config.ad_classification_dir, filename);
        
        if ~exist(filepath, 'file')
            warning('AD/noAD classification not found: %s', filename);
            continue;
        end
        
        ad_array(:, :, i) = imread(filepath);
    end
    
    fprintf('    Loaded %d classification maps\n', n_years);
end

%% Core Extraction Algorithm

function [ad_metrics, noad_metrics] = extractPairedMetrics(...
    ad_classification, trend_data, resilience_data, ...
    window_size, config, n_rows, n_cols)
    % Extract paired AD and noAD metrics
    %
    % For each pixel:
    %   1. Find first year of AD (if any)
    %   2. Extract δTAC and TAC at t-1
    %   3. Find noAD pixels in neighborhood at same time
    %   4. Calculate mean noAD metrics
    
    % Initialize output structures
    ad_metrics = struct();
    ad_metrics.year = NaN(n_rows, n_cols);        % Year of first AD
    ad_metrics.deltaTAC = NaN(n_rows, n_cols);     % δTAC before AD
    ad_metrics.TAC = NaN(n_rows, n_cols);         % TAC before AD
    
    noad_metrics = struct();
    noad_metrics.count = NaN(n_rows, n_cols);     % Number of noAD neighbors
    noad_metrics.deltaTAC = NaN(n_rows, n_cols);   % Mean δTAC of noAD
    noad_metrics.TAC = NaN(n_rows, n_cols);       % Mean TAC of noAD
    
    % Window parameters
    radius = floor(window_size / 2);
    
    % Progress tracking
    total_pixels = n_rows * n_cols;
    processed = 0;
    last_percent = 0;
    
    % Process each pixel
    for row = 1:n_rows
        for col = 1:n_cols
            
            % Find first AD year for this pixel
            [ad_occurred, ad_year_idx, ad_year] = findFirstAD(...
                ad_classification(row, col, :), config.analysis_start);
            
            if ~ad_occurred
                continue;  % No AD event, skip
            end
            
            % Store AD year
            ad_metrics.year(row, col) = ad_year;
            
            % Extract AD metrics at t-1
            ad_metrics.deltaTAC(row, col) = trend_data(row, col, ad_year_idx);
            ad_metrics.TAC(row, col) = resilience_data(row, col, ad_year_idx);
            
            % Find noAD pixels in neighborhood
            [noad_count, mean_deltaTAC, mean_TAC] = findNeighborhoodNoAD(...
                row, col, ad_year_idx, radius, ...
                ad_classification, trend_data, resilience_data, ...
                n_rows, n_cols);
            
            % Store noAD metrics
            if noad_count > 0
                noad_metrics.count(row, col) = noad_count;
                noad_metrics.deltaTAC(row, col) = mean_deltaTAC;
                noad_metrics.TAC(row, col) = mean_TAC;
            end
            
            % Progress update
            processed = processed + 1;
            current_percent = floor(100 * processed / total_pixels);
            if current_percent > last_percent && mod(current_percent, 100) == 0
                fprintf('    Progress: %d%%\n', current_percent);
                last_percent = current_percent;
            end
        end
    end
    
    % Summary statistics
    n_ad_pixels = sum(~isnan(ad_metrics.year(:)));
    n_paired = sum(~isnan(noad_metrics.count(:)));
    
    fprintf('    AD pixels: %d, Successfully paired: %d (%.1f%%)\n', ...
            n_ad_pixels, n_paired, 100 * n_paired / max(1, n_ad_pixels));
end

function [occurred, year_idx, year_value] = findFirstAD(pixel_timeseries, start_year)
    % Find the first year when AD occurred for a pixel
    %
    % Input:
    %   pixel_timeseries - 1×n_years array of AD classifications
    %   start_year - First year of analysis period
    %
    % Output:
    %   occurred - Boolean, true if AD happened
    %   year_idx - Index of first AD year (1-based)
    %   year_value - Calendar year of first AD
    
    occurred = false;
    year_idx = NaN;
    year_value = NaN;
    
    % Find first occurrence of AD (value = 1)
    ad_indices = find(pixel_timeseries == 1);
    
    if ~isempty(ad_indices)
        occurred = true;
        year_idx = ad_indices(1);
        year_value = start_year + year_idx - 1;
    end
end

function [noad_count, mean_deltaTAC, mean_TAC] = findNeighborhoodNoAD(...
    center_row, center_col, time_idx, radius, ...
    ad_classification, trend_data, resilience_data, ...
    n_rows, n_cols)
    % Find noAD pixels in spatial neighborhood at specific time
    %
    % Returns mean δTAC and TAC for noAD pixels
    
    % Define window bounds
    row_min = max(1, center_row - radius);
    row_max = min(n_rows, center_row + radius);
    col_min = max(1, center_col - radius);
    col_max = min(n_cols, center_col + radius);
    
    % Extract neighborhood data
    neighborhood_ad = ad_classification(row_min:row_max, col_min:col_max, time_idx);
    neighborhood_trend = trend_data(row_min:row_max, col_min:col_max, time_idx);
    neighborhood_resilience = resilience_data(row_min:row_max, col_min:col_max, time_idx);
    
    % Identify noAD pixels (value = 0)
    noad_mask = (neighborhood_ad == 0);
    
    % Count and calculate means
    noad_count = sum(noad_mask(:));
    
    if noad_count > 0
        noad_trend_values = neighborhood_trend(noad_mask);
        noad_resilience_values = neighborhood_resilience(noad_mask);
        
        mean_deltaTAC = mean(noad_trend_values, 'omitnan');
        mean_TAC = mean(noad_resilience_values, 'omitnan');
    else
        mean_deltaTAC = NaN;
        mean_TAC = NaN;
    end
end

%% Save Results

function savePairedResults(ad_metrics, noad_metrics, ...
                          sigma_n, sigma_m, window_size, ...
                          config, spatial_ref, geo_info, output_dirs)
    % Save all paired metrics as GeoTIFF files
    
    % Construct base filename
    base_name = sprintf('sigma%ddiv%d_WS%02d', sigma_n, sigma_m, window_size);
    
    % Construct full filenames with descriptive components
    trend_suffix = sprintf('_%s_%s_%d_%s', ...
                          config.vegetation_index, ...
                          config.resilience_indicator, ...
                          config.temporal_window, ...
                          config.decomposition_method);
    
    % Save AD metrics
    saveGeoTIFF(ad_metrics.year, ...
               fullfile(output_dirs.ADyear, ['AD_', base_name, '.tif']), ...
               spatial_ref, geo_info);
    
    saveGeoTIFF(ad_metrics.deltaTAC, ...
               fullfile(output_dirs.ADdeltaTAC, ['AD_', base_name, '_Sen', trend_suffix, '.tif']), ...
               spatial_ref, geo_info);
    
    saveGeoTIFF(ad_metrics.TAC, ...
               fullfile(output_dirs.ADTAC, ['AD_', base_name, trend_suffix, '.tif']), ...
               spatial_ref, geo_info);
    
    % Save noAD metrics
    saveGeoTIFF(noad_metrics.count, ...
               fullfile(output_dirs.noADcount, ['noAD_', base_name, '.tif']), ...
               spatial_ref, geo_info);
    
    saveGeoTIFF(noad_metrics.deltaTAC, ...
               fullfile(output_dirs.noADdeltaTAC, ['noAD_', base_name, '_Sen', trend_suffix, '.tif']), ...
               spatial_ref, geo_info);
    
    saveGeoTIFF(noad_metrics.TAC, ...
               fullfile(output_dirs.noADTAC, ['noAD_', base_name, trend_suffix, '.tif']), ...
               spatial_ref, geo_info);
    
    fprintf('    Saved 6 output files\n');
end

function saveGeoTIFF(data, filepath, spatial_ref, geo_info)
    % Save data as GeoTIFF with proper georeferencing
    
    geotiffwrite(filepath, data, spatial_ref, ...
                'GeoKeyDirectoryTag', geo_info.GeoTIFFTags.GeoKeyDirectoryTag);
end

%% Utility Functions

function [data, spatial_ref, geo_info] = loadReferenceData(filepath)
    % Load reference GeoTIFF for spatial information
    
    if ~exist(filepath, 'file')
        error('Reference file not found: %s', filepath);
    end
    
    [data, spatial_ref] = geotiffread(filepath);
    geo_info = geotiffinfo(filepath);
end
