% =====================================================================
% EEG Time-Resolved ISC Computation Script
% Author: Juncheng
% Date: Sep 2024
% Description:
%   This script computes time-resolved Inter-Subject Correlation (ISC)
%   for multiple EEG video stimuli using a sliding window approach.
%   Results are saved for further analysis.
% =====================================================================

clear all; close all; clc;

%% ========================== Initialization ==========================
% Define main paths and parameters
main_path = '<main_results_directory>'; % Base directory for results
time_resolved_isc_output_path = [main_path '/08_time_resolved_isc_results/50nd/'];

% Ensure output directory exists
if ~exist(time_resolved_isc_output_path, 'dir')
    mkdir(time_resolved_isc_output_path);
    disp('Output folder created.');
end

% Define paths to EEG data files (paths anonymized for sharing purposes)
spot_path_names = { ...
    '<path_to_combined_struct>/EEGVolume_commercial_aribnb_30s.mat', ...
    '<path_to_combined_struct>/EEGVolume_commercial_att_30s.mat', ...
    ... % Add more file paths as needed
};
n_spots = length(spot_path_names); % Number of video stimuli

% Define parameters for time-resolved ISC
fs = 250;               % Sampling frequency (Hz)
Nsec = 2;               % Window length (seconds)
ISCres = 0.2;           % Time resolution for sliding window (seconds)

%% ==================== Process Each Video Spot ====================
for current_spot = 1:n_spots
    fprintf('Processing video %d of %d...\n', current_spot, n_spots);

    % Load EEG data for the current video spot
    EEGVolume_path = spot_path_names{current_spot};
    load(EEGVolume_path, 'X'); % Assuming 'X' is the EEG data matrix

    [T, D, N] = size(X); % T: time points, D: channels, N: subjects
    fprintf('Loaded EEG data: Time Points = %d, Channels = %d, Subjects = %d\n', T, D, N);

    %% ======== Precompute Global ISC Components ========
    fprintf('Computing global ISC components for video %d...\n', current_spot);

    % Compute global covariance matrices
    Rij_global = permute(reshape(cov(X(:,:)), [D N D N]), [1 3 2 4]);
    Rw_global = 1 / N * sum(Rij_global(:,:,1:N+1:N*N), 3);
    Rb_global = 1 / (N - 1) / N * (sum(Rij_global(:,:,:), 3) - N * Rw_global);

    % Regularize within-subject covariance
    gamma = 0.1; % Regularization parameter
    Rw_reg_global = (1 - gamma) * Rw_global + gamma * mean(eig(Rw_global)) * eye(size(Rw_global));

    % Compute correlated components (W)
    [W, ~] = eig(Rb_global, Rw_reg_global);
    [~, idx] = sort(diag(W), 'descend'); % Sort components by ISC
    W = W(:, idx);

    %% ======== Compute Time-Resolved ISC ========
    fprintf('Computing time-resolved ISC for video %d...\n', current_spot);

    % Initialize time vector and ISC storage
    timevec = 1:ISCres:(T / fs); % Time vector
    ISC_perTime = zeros(size(W, 2), length(timevec)); % ISC for each component and time window

    for t = 1:length(timevec)
        % Define time window
        start_idx = max(1, round((t - 1) * ISCres * fs - (Nsec * fs) / 2));
        end_idx = min(T, start_idx + Nsec * fs - 1);
        Xt = X(start_idx:end_idx, :, :);

        % Zero-pad if window is smaller than required
        if size(Xt, 1) < Nsec * fs
            pad_size = Nsec * fs - size(Xt, 1);
            Xt = [Xt; zeros(pad_size, D, N)];
        end

        % Compute covariance matrices for the current time window
        Rij = permute(reshape(cov(Xt(:,:)), [D N D N]), [1 3 2 4]);
        Rw = 1 / N * sum(Rij(:,:,1:N+1:N*N), 3);
        Rb = 1 / (N - 1) / N * (sum(Rij(:,:,:), 3) - N * Rw);

        % Compute ISC for each component
        ISC_perTime(:, t) = diag(W' * Rb * W) ./ diag(W' * Rw * W);
    end

    %% ======== Save Results ========
    output_file = fullfile(time_resolved_isc_output_path, sprintf('ISC_time_resolved_video_%d.mat', current_spot));
    save(output_file, 'ISC_perTime', 'timevec', 'W');
    fprintf('Time-resolved ISC saved for video %d: %s\n', current_spot, output_file);
end

%% ========================= Final Cleanup =========================
disp('All tasks completed. Time-resolved ISC results saved in output directory.');


    % Save time-resolved ISC results
    time_resolved_filename = fullfile(time_resolved_isc_output_path, ...
        sprintf('time_resolved_isc_spot_%d.mat', current_spot));
    save(time_resolved_filename, 'ISC_perTime', 'timevec');

    fprintf('Time-resolved ISC saved for video %d.\n', current_spot);
end

fprintf('All videos processed and time-resolved ISC results saved.\n');
