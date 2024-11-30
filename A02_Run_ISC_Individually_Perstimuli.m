% =====================================================================
% EEG ISC Analysis Script
% Author: Juncheng
% Date: Sep 2024
% Description:
%   This script processes preprocessed EEG data, computes covariance 
%   matrices for individual videos, averages within- and between-subject 
%   covariance matrices, and calculates ISC (Inter-Subject Correlation).
%   Results include forward projections (A), correlated components (W), 
%   and ISC values resolved by subjects and stimuli.
% =====================================================================

clear all; close all; clc;

%% ========================== Initialization ==========================
% Define paths (hidden for sharing purposes)
main_script_dir = '<path_to_scripts>'; % Directory containing this script
data_dir = '<path_to_combined_struct>'; % Directory with EEG data
output_dir = '<path_to_results>'; % Output directory for results

cd(main_script_dir); % Change to the script directory

% Define EEG data file paths (replace with actual file paths)
spot_path_names = {
    '<path_to_data>/EEGVolume_commercial_aribnb.mat';
    '<path_to_data>/EEGVolume_commercial_att.mat';
    ... % Add more files as needed
};
n_spots = length(spot_path_names); % Number of videos (stimuli)

%% ==================== Process Each Video Stimulus ====================
disp('Loading and processing EEG data for each stimulus...')

for current_spot = 1:n_spots
    % Load EEG data
    EEGVolume_path = spot_path_names{current_spot};
    load(EEGVolume_path, 'X'); % Load EEG data matrix (X)
    
    % Compute covariance matrix for each spot
    [~, D, N] = size(X); % Dimensions: D = channels, N = subjects
    Rij = permute(reshape(cov(X(:,:)), [D N D N]), [1 3 2 4]);
    
    % Save covariance matrix
    save_path = fullfile(output_dir, sprintf('single_cov_mat_spot_%d.mat', current_spot));
    save(save_path, 'Rij', '-v7.3');
    fprintf('Covariance matrix saved: %s\n', save_path);
end

clear EEGVolume_path X Rij current_spot

%% ===================== Load and Combine Cov Matrices ====================
disp('Loading individual covariance matrices...')

cd(output_dir); % Change to the output directory
file_list = dir('single_cov_mat_spot_*.mat'); % Get all covariance matrix files

Rij_Data = cell(length(file_list), 1); % Initialize storage
for i = 1:length(file_list)
    load(file_list(i).name, 'Rij');
    Rij_Data{i} = Rij; % Store loaded covariance matrix
end

% Get common parameters
Dcommon = size(Rij_Data{1}, 2); % Number of sensors
Nsubs = size(Rij_Data{1}, 3); % Number of subjects
Nstim = length(Rij_Data); % Number of stimuli
gamma = 0.1; % Regularization parameter

%% ===================== Compute Within-Subjects Covariance ====================
disp('Computing within-subjects covariance...')

Rw = 0; % Initialize within-subject covariance matrix
for iVid = 1:Nstim
    Rij = Rij_Data{iVid};
    Rw = Rw + 1 / Nsubs * sum(Rij(:, :, 1:Nsubs + 1:Nsubs * Nsubs), 3) / Nstim;
end

disp('Within-subject covariance computed.')

%% ==================== Compute Between-Subjects Covariance ====================
disp('Computing between-subjects covariance...')

Rb = 0; % Initialize between-subject covariance matrix
for iVid = 1:Nstim
    Rij = Rij_Data{iVid};
    Rb = Rb + 1 / (Nsubs - 1) * 1 / Nsubs * (sum(Rij(:,:,:), 3) - Nsubs * Rw);
end
Rb = Rb / Nstim; % Average across stimuli

disp('Between-subject covariance computed.')

%% ===================== Regularize Within-Subjects Covariance ====================
disp('Regularizing within-subjects covariance...')
Rw_reg = (1 - gamma) * Rw + gamma * mean(eig(Rw)) * eye(size(Rw));
disp('Regularization done.')

%% ======================== Save Covariance Results =========================
save(fullfile(output_dir, 'Rw_avg.mat'), 'Rw');
save(fullfile(output_dir, 'Rb_avg.mat'), 'Rb');
save(fullfile(output_dir, 'Rw_reg_avg.mat'), 'Rw_reg');
disp('Covariance results saved.')

%% ===================== Compute Correlated Components =====================
disp('Computing correlated components...')
[W, ISC] = eig(Rb, Rw_reg); % Compute eigenvectors and eigenvalues
[ISC, idx] = sort(diag(ISC), 'descend'); % Sort by ISC
W = W(:, idx); % Sort eigenvectors accordingly
A = Rw * W / (W' * Rw * W); % Compute forward projections

disp('Correlated components computed.')

%% ==================== Compute ISC Resolved by Subjects ====================
disp('Computing ISC per subject...')
ISC_persubject = zeros(size(W, 2), Nsubs, Nstim);
for iVid = 1:Nstim
    Rij = Rij_Data{iVid};
    for i = 1:Nsubs
        Rw_sub = 0;
        Rb_sub = 0;
        for j = 1:Nsubs
            if i ~= j
                Rw_sub = Rw_sub + 1 / (Nsubs - 1) * (Rij(:, :, i, i) + Rij(:, :, j, j));
                Rb_sub = Rb_sub + 1 / (Nsubs - 1) * (Rij(:, :, i, j) + Rij(:, :, j, i));
            end
        end
        ISC_persubject(:, i, iVid) = diag(W' * Rb_sub * W) ./ diag(W' * Rw_sub * W);
    end
end

disp('ISC per subject computed.')

%% ========================== Save ISC Results ==========================
save(fullfile(output_dir, 'ISC_results.mat'), 'ISC', 'ISC_persubject', 'W', 'A');
disp('ISC results saved.')

%% ============================ Cleanup ==============================
clearvars -except output_dir
disp('Processing complete. Results saved in the output directory.')

