% =====================================================================
% EEG Data Preprocessing Script
% Author: Juncheng
% Date: Sep 2024
% Description:
%   This script processes EEG data collected under different conditions,
%   including loading, preprocessing (e.g., filtering, EOG correction),
%   and saving the results for further analysis. 
% =====================================================================

clear all; close all; clc;

%% ========================== Initialization ==========================
% Define paths (hidden for sharing purposes)
top_dir_name = '<main_data_directory>'; % Main directory for input EEG data
addpath('<path_to_EEGLab>');            % Path to EEGLab library
addpath('<path_to_custom_scripts>');    % Path to custom scripts

% Set condition and output directory (adjust as needed)
condition = '100';
out_dir_path = '<output_directory>'; 

% Load empty EEG volume structure for initialization
load('EEGVolumeEmpty.mat');

% Define spot names (adjust list as necessary)
spot_names = {  'commercial_aribnb_30s';
                'commercial_att_30s';
                'commercial_carscom_30s';
                'commercial_cookies_30s';
                'commercial_dominos_30s';
                'commercial_doritos_30s';
                'commercial_expedia_30s';
                'commercial_google_pixel_30s';
                'commercial_hr_block_30s';
                'commercial_jersey_mikes_30s';
                'commercial_lego_30s';
                'commercial_meta_quest_2_30s';
                'commercial_milk_30s';
                'commercial_progressive_30s';
                'commercial_publix_30s';
                'commercial_puma_30s';
                'commercial_starbucks_30s';
                'commercial_under_armour_30s';
                'health_alcohol_30s';
                'health_alzheimers_30s';
                'health_covid_vaccine_30s';
                'health_diet_30s';
                'health_drunk_driving_30s';
                'health_fitness_30s';
                'health_kidney_30s';
                'health_mantherapy_30s';
                'health_prediabetes_30s';
                'health_stroke_30s';
                'health_vaping_30s';
                'health_weight_30s'};
n_spots = length(spot_names); % Number of spots

%% ====================== Processing Each Spot =======================
for current_spot = 1:n_spots
    % Search for files matching the current spot name and condition
    search_string = strcat('/*/*', char(spot_names(current_spot)), '*', condition, '*.set');
    file_List = dir([top_dir_name search_string]);
    n_files = length(file_List);   % Total number of files for the current spot
    
    fprintf("Processing Spot: %s --- Total Files: %i\n", char(spot_names(current_spot)), n_files);
    
    % Initialize storage for EEG data
    n_samples = 7250;   % Number of time samples
    n_channels = 24;    % Number of EEG channels
    result_eeg = zeros(n_channels, n_samples, n_files);
    
    % Iterate over all files for the current spot
    for curr_sub = 1:n_files
        current_file_name = file_List(curr_sub).name;
        current_file_path = file_List(curr_sub).folder;
        
        % Load EEG data using EEGLab function
        EEG = pop_loadset('filename', current_file_name, 'filepath', current_file_path, 'verbose', 'off', 'check', 'off');
        fprintf("Subject %i --- Channels: %i, Samples: %i\n", curr_sub, size(EEG.data, 1), size(EEG.data, 2));
        
        % Check for consistent sample size
        if size(EEG.data, 2) ~= n_samples
            warning('Mismatch in sample size for file: %s', current_file_name);
        else
            result_eeg(:,:,curr_sub) = EEG.data;
        end
    end
    
    % Reshape data for preprocessing
    X = permute(result_eeg, [2, 1, 3]);
    
    %% ======================= Preprocessing ==========================
    % Preprocessing Parameters
    opts.fs = 250;            % Sampling rate
    opts.eogchannels = [1 2 11 12]; % Channels for EOG correction (Fp1, Fp2, F7, F8)
    opts.badchannels = -1;    % Use default bad channel detection
    opts.kIQD = 3;            % IQR multiplier for outlier detection
    opts.SaveBadChs = 0;      % Do not save bad channel info
    
    % Run preprocessing function (custom implementation required)
    disp('Performing preprocessing...');
    X = nceeg_preprocess_v9c(X, opts.eogchannels, opts.badchannels, opts, char(spot_names(current_spot)));
    disp('Preprocessing complete.');
    
    %% ==================== Save Processed Data =======================
    out_file_name = strcat('EEGVolume_', char(spot_names(current_spot)), '_', condition, '_', string(n_files), '_viewers_supermatrix.mat');
    out_name = strcat(out_dir_path, out_file_name);
    
    save(out_name, 'X');
    fprintf('Processed data saved to: %s\n', out_name);
end
