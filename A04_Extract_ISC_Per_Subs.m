% =====================================================================
% EEG ISC Per Participant Extraction Script
% Author: Juncheng
% Date: Sep 2024
% Description:
%   This script extracts ISC (Inter-Subject Correlation) values for 
%   Component 1 (C1) across all participants and videos, formats them 
%   into a table, and saves the results as a CSV file.
% =====================================================================

clearvars; close all; clc;

%% ========================== Initialization ==========================
% Set main paths (paths anonymized for sharing purposes)
opts.mainpath = '<main_results_directory>';

% Set input and output paths
opts.EOG        = 'EOG4'; % EOG channel configuration
opts.Preprostr  = ['_noFilt-Epoch_250Hz_BadCh-1_' opts.EOG '_PreproV9c'];
opts.pathin     = [opts.mainpath '\05_isc_results\']; % Input directory for ISC results
opts.pathout    = [opts.mainpath '\07_extract_isc_per_subs\']; % Output directory for extracted data

% Ensure output directory exists
if ~exist(opts.pathout, 'dir')
    mkdir(opts.pathout);
    disp('Output folder created.');
end

% Set EOG location files
if strcmp(opts.EOG, 'EOG2')
    opts.ChannelLocs = [opts.mainpath '\scripts\NCMobi_Smarting_24Ch_EOG2_locationfile.mat'];
    opts.EOGLocs = [1 2];
elseif strcmp(opts.EOG, 'EOG4')
    opts.ChannelLocs = [opts.mainpath '\scripts\NCMobi_Smarting_24Ch_EOG4_locationfile.mat'];
    opts.EOGLocs = [1 2 11 12];
else
    error('Check EOG definition');
end

% Define parameters
opts.nSubs      = 12; % Number of participants
opts.nComps     = 3;  % Number of components to extract
opts.nVids      = 30; % Number of videos
Ncomp           = 3;  % Number of components
conds           = {'all'}; % Condition labels

%% ==================== Load ISC Data and Channel Locations ====================
disp('Loading ISC data and channel locations...');
load([opts.pathin 'ISC_results_50wd.mat'], 'A', 'ISC', 'W', 'ISC_persubject'); % Load ISC results
load(opts.ChannelLocs, 'EEG'); % Load EEG channel locations
disp('Data loaded successfully.');

%% ===================== Extract Data for Component 1 =====================
disp('Extracting ISC values for Component 1 (C1)...');

% Prepare data structure for all videos
for j = 1:opts.nVids
    DataISC_persubject(j, :, :) = ISC_persubject{j}; % Store ISC values
end
disp('.done.');

% Extract ISC values for Component 1 (C1)
C1_values = squeeze(DataISC_persubject(:, 1, :)); % [Videos × Participants]
disp('C1 ISC values extracted:');
disp(C1_values);

%% ==================== Format and Save C1 Data =====================
% Format extracted values as a table
C1_perParticipant = array2table(C1_values, ...
    'VariableNames', arrayfun(@(x) ['Participant_' num2str(x)], 1:size(C1_values, 2), 'UniformOutput', false), ...
    'RowNames', arrayfun(@(x) ['Video_' num2str(x)], 1:size(C1_values, 1), 'UniformOutput', false));

disp('C1 ISC values formatted as a table:');
disp(C1_perParticipant);

% Save table to a CSV file
output_file = fullfile(opts.pathout, 'C1_Participant_50wd.csv');
writetable(C1_perParticipant, output_file, 'WriteRowNames', true);
disp(['C1 data saved as CSV file: ', output_file]);

%% ========================= Final Cleanup =========================
clearvars -except opts;
disp('All tasks completed. Results saved in output directory.');
