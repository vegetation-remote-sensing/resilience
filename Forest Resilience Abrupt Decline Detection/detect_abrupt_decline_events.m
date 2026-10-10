function detect_abrupt_decline_events(config)
% DETECT_ABRUPT_DECLINE_EVENTS Identifies abrupt decline events in vegetation
%
% Description:
%   This function detects abrupt decline (AD) events from anomalies in
%   growing-season kNDVI (kernel Normalized Difference Vegetation Index).
%   For every pixel, the annual series is linearly detrended and compared
%   with the mean (mu) and standard deviation (sigma) of a baseline period:
%       AD   (value 1): kNDVI < mu - n*sigma
%       noAD (value 0): kNDVI > mu - n*sigma/m
%       undefined (NaN): values between the two thresholds, masked pixels
%
% Workflow:
%   Step 1. Validate the configuration and set defaults.
%   Step 2. Define the baseline and analysis periods.
%   Step 3. Load the reference grid and identify the valid pixels.
%   Step 4. Create the output directories.
%   Step 5. Load the annual growing-season kNDVI maps.
%   Step 6. Remove the linear trend of every pixel and save the result.
%   Step 7. Calculate the baseline mean and standard deviation.
%   Step 8. Calculate the AD/noAD thresholds and classify every year of the
%           analysis period, for every combination of n and m.
%
% Input:
%   config - Structure containing configuration parameters:
%       .input_dir          : Directory containing input kNDVI data
%       .output_dir         : Directory for output results
%       .reference_file     : Path to reference GeoTIFF for spatial info
%                             (land mask; value 0 = valid pixel)
%       .baseline_start     : Start year of baseline period (e.g., 2000)
%       .baseline_end       : End year of baseline period (e.g., 2009)
%       .analysis_start     : Start year of analysis period (e.g., 2010)
%       .analysis_end       : End year of analysis period (e.g., 2022)
%       .growing_season     : String describing season (e.g., 'GS_05_09')
%       .sigma_values       : AD threshold n (default: 1:6)
%       .sigma_denominator  : denominator of noAD threshold m (default: 2:6)
%       .vegetation_index   : Index name (default: 'kNDVI')
%
%   Input file naming: [index]_[season]_[year].tif
%
% Output:
%   - detrended_kNDVI_growing_season/ : Detrended kNDVI GeoTIFF files
%   - sigma_thresholds/               : AD and noAD threshold GeoTIFF files
%   - AD_noAD_classification/         : AD (1) / noAD (0) GeoTIFF files
%
% Notes:
%   - The linear trend is fitted over the whole period
%     (baseline_start to analysis_end).
%   - sigma is the population standard deviation of the baseline period
%     (normalised by N, not N-1).
%
% Requirements: MATLAB R2020a or later, Mapping Toolbox,
%               Image Processing Toolbox
%
% Author: Yiling Cai
% Institution: School of Geography and Remote Sensing, Guangzhou University, Guangzhou 510006, China;  Guangdong Guodi Science Technology Co., Ltd, Guangzhou 510075, China;
% Date: 2026-01-16
% Version: 1.1

    %% Step 1: Validate the configuration and set defaults
    config = validateADConfig(config);

    %% Step 2: Define the baseline and analysis periods
    baseline_years = config.baseline_start:config.baseline_end;
    analysis_years = config.analysis_start:config.analysis_end;
    all_years = config.baseline_start:config.analysis_end;

    fprintf('Processing AD detection from %d to %d\n', ...
            config.analysis_start, config.analysis_end);
    fprintf('Baseline period: %d-%d\n', ...
            config.baseline_start, config.baseline_end);

    %% Step 3: Load the reference grid and identify the valid pixels
    [ref_data, spatial_ref, geo_info] = loadReferenceData(config.reference_file);
    [n_rows, n_cols] = size(ref_data);
    n_pixels = n_rows * n_cols;

    % Valid pixels are marked with 0 in the reference (land mask) file
    valid_pixels = find(ref_data(:) == 0);
    fprintf('Valid pixels: %d (%.2f%%)\n', ...
            length(valid_pixels), 100*length(valid_pixels)/n_pixels);

    %% Step 4: Create the output directories
    dir_detrend = fullfile(config.output_dir, 'detrended_kNDVI_growing_season');
    dir_sigma = fullfile(config.output_dir, 'sigma_thresholds');
    dir_ad = fullfile(config.output_dir, 'AD_noAD_classification');
    createOutputDirs({dir_detrend, dir_sigma, dir_ad});

    %% Step 5: Load the annual growing-season kNDVI maps
    % Result: matrix of size [n_years x n_pixels]
    fprintf('\n=== Loading kNDVI data ===\n');
    kndvi_data = loadTimeSeriesData(config, all_years, n_rows, n_cols);

    %% Step 6: Remove the linear trend of every pixel and save the result
    fprintf('\n=== Detrending time series ===\n');
    detrended_data = detrendTimeSeries(kndvi_data, valid_pixels);

    saveDetrended(detrended_data, all_years, dir_detrend, ...
                  config, n_rows, n_cols, spatial_ref, geo_info);

    %% Step 7: Calculate the baseline mean and standard deviation
    fprintf('\n=== Computing baseline statistics ===\n');
    baseline_idx = 1:length(baseline_years);
    baseline_data = detrended_data(baseline_idx, :);

    mu = mean(baseline_data, 1, 'omitnan');
    sigma = std(baseline_data, 1, 1, 'omitnan');  % population SD (divide by N)

    %% Step 8: Calculate the thresholds and classify AD/noAD
    fprintf('\n=== Detecting AD events ===\n');
    detectADEvents(detrended_data, mu, sigma, analysis_years, ...
                   config, n_rows, n_cols, ...
                   spatial_ref, geo_info, dir_sigma, dir_ad);

    fprintf('\nProcessing complete!\n');
end

%% Helper Functions

function config = validateADConfig(config)
    % Check the required fields and fill in the default values

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

    % The periods must be in chronological order
    if config.baseline_end < config.baseline_start || ...
       config.analysis_start < config.baseline_start || ...
       config.analysis_end < config.analysis_start
        error(['Periods must satisfy baseline_start <= baseline_end, ', ...
               'baseline_start <= analysis_start <= analysis_end.']);
    end
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

function data_matrix = loadTimeSeriesData(config, years, n_rows, n_cols)
    % Load the annual kNDVI maps into a [n_years x n_pixels] matrix
    n_years = length(years);
    data_matrix = NaN(n_years, n_rows * n_cols);

    for i = 1:n_years
        year = years(i);
        filename = sprintf('%s_%s_%04d.tif', ...
                          config.vegetation_index,config.growing_season, year);
        filepath = fullfile(config.input_dir, filename);

        % Every year is required: a missing year would bias the trend
        % and the baseline statistics
        if ~exist(filepath, 'file')
            error('File not found: %s', filepath);
        end

        temp_data = imread(filepath);
        if ~isequal(size(temp_data), [n_rows, n_cols])
            error('Grid size of %s does not match the reference file.', filename);
        end
        data_matrix(i, :) = temp_data(:);

        if mod(i, 5) == 0
            fprintf('  Loaded %d/%d files\n', i, n_years);
        end
    end

    fprintf('  Loaded %d/%d files total\n', n_years, n_years);
end

function detrended = detrendTimeSeries(data, valid_pixels)
    % Remove the linear trend (least-squares line) from every valid pixel
    [n_years, n_pixels] = size(data);
    detrended = NaN(n_years, n_pixels);

    % Report valid pixels with missing years (they cannot be classified)
    n_incomplete = sum(any(isnan(data(:, valid_pixels)), 1));
    if n_incomplete > 0
        warning('%d valid pixels contain missing values.', n_incomplete);
    end

    % Detrend only valid pixels (each column is one pixel)
    detrended(:, valid_pixels) = detrend(data(:, valid_pixels));

    fprintf('  Detrended %d pixels\n', length(valid_pixels));
end

function saveDetrended(data, years, output_dir, config, ...
                       n_rows, n_cols, spatial_ref, geo_info)
    % Save the detrended data as one GeoTIFF file per year
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
    % Detect AD events for every combination of n (severity) and m (buffer)

    % Row index of the first analysis year in the data matrix
    year_offset = config.analysis_start - config.baseline_start;

    for n = config.sigma_values
        % Step 8.1: AD threshold (mu - n*sigma)
        threshold_ad = mu - sigma * n;
        saveThreshold(threshold_ad, n, 0, dir_sigma, config, ...
                     n_rows, n_cols, spatial_ref, geo_info);

        for m = config.sigma_denominator
            % Step 8.2: noAD threshold (mu - n*sigma/m)
            threshold_noad = mu - sigma * n / m;
            saveThreshold(threshold_noad, n, m, dir_sigma, config, ...
                         n_rows, n_cols, spatial_ref, geo_info);

            % Step 8.3: Classify every year of the analysis period
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
    % Classify pixels as AD (1), noAD (0) or undefined (NaN)
    ad_vector = NaN(size(data_vector));

    % AD event: below the lower threshold
    ad_vector(data_vector < thresh_ad) = 1;

    % No AD: above the higher threshold
    ad_vector(data_vector > thresh_noad) = 0;

    % Reshape to 2D
    ad_map = reshape(ad_vector, n_rows, n_cols);
end

function saveThreshold(threshold, n, m, output_dir, config, ...
                       n_rows, n_cols, spatial_ref, geo_info)
    % Save a threshold map (m = 0: AD threshold; m > 0: noAD threshold)
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
    % Save an AD/noAD classification map
    filename = sprintf('AD_noAD_classification_sigma%ddiv%d_%04d.tif', n, m, year);
    filepath = fullfile(output_dir, filename);

    geotiffwrite(filepath, ad_map, spatial_ref, ...
                'GeoKeyDirectoryTag', geo_info.GeoTIFFTags.GeoKeyDirectoryTag);
end
