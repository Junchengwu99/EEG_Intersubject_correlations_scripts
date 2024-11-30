% =====================================================================
% EEG ISC Component Plotting Script
% Author: Juncheng
% Date: Sep 2024
% Description:
%   This script generates topographical plots for ISC (Inter-Subject 
%   Correlation) components and saves the results as PDF files.
% =====================================================================

clearvars; close all; clc;

%% ========================== Initialization ==========================
% Set main paths (paths anonymized for sharing purposes)
opts.mainpath = '<main_results_directory>';
addpath('<path_to_eeglab_library>');

% Set input and output paths
opts.EOG        = 'EOG4'; % EOG channel definition
opts.Preprostr  = ['_noFilt-Epoch_250Hz_BadCh-1_' opts.EOG '_PreproV9c'];
opts.pathin     = [opts.mainpath '\05_isc_results\']; % Input directory for ISC results
opts.pathout    = [opts.mainpath '\06_plot_extract_components\']; % Output directory for plots

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
opts.nSubs      = 12; % Number of subjects
opts.nComps     = 3;  % Number of components to plot
opts.nVids      = 30; % Number of videos
opts.Scale      = [-0.2 0.2]; % Color scale for topoplot
opts.recode     = 1; % Recode components for visualization
opts.recodeVal  = 24; % Reference sensor
opts.print      = 1; % Enable printing of plots
opts.eps        = 1; % Save as EPS
Ncomp           = 3; % Number of components to plot
conds           = {'all'}; % Condition labels

%% ==================== Load ISC Data and Channel Locations ====================
disp('Loading ISC data and channel locations...');
load([opts.pathin 'ISC_results_50wd.mat'], 'A', 'ISC', 'W', 'ISC_persubject'); % Load ISC results
load(opts.ChannelLocs, 'EEG'); % Load EEG channel locations
disp('Data loaded successfully.');

%% ===================== Compute ISC per Video =====================
disp('Computing ISC per video...');
for j = 1:opts.nVids
    DataISC_persubject(j, :, :) = ISC_persubject{j}; % Store ISC data per video
end
disp('.done.');

% Compute mean ISC for each component across videos
ISCperVid_C1 = mean(squeeze(DataISC_persubject(:, 1, :)), 2);
ISCperVid_C2 = mean(squeeze(DataISC_persubject(:, 2, :)), 2);
ISCperVid_C3 = mean(squeeze(DataISC_persubject(:, 3, :)), 2);

% Save ISC per video data
save(fullfile(opts.pathout, 'ISC_per_video_data_50wd.mat'), 'ISCperVid_C1', 'ISCperVid_C2', 'ISCperVid_C3');
disp('ISC per video data saved.');

%% ==================== Plot ISC Components ====================
disp('Plotting ISC components...');
figure;
set(gcf, 'Position', [100, 100, 1400, 800]); % Set figure size

for iComp = 1:Ncomp
    subplot(1, Ncomp, iComp); % Create subplot for each component
    hold on;

    % Adjust polarity if recoding is enabled
    if opts.recode == 1
        recodeSens = opts.recodeVal - sum(opts.EOGLocs < opts.recodeVal);
        if A(recodeSens, iComp) <= 0
            Ax = -A(:, iComp); % Invert component polarity
        else
            Ax = A(:, iComp);
        end
    else
        Ax = A(:, iComp);
    end

    % Plot the topographical representation of the component
    topoplot(Ax, EEG.chanlocs, 'electrodes', 'on', 'maplimits', opts.Scale, ...
        'colormap', parula, 'whitebk', 'on', 'style', 'both', 'conv', 'off', ...
        'emarker', {'.', 'k', 6, 1}, 'headrad', 0.5, 'numcontour', 5);

    % Add title to each component plot
    title(['C_', num2str(iComp), ' - ISC 50wd']);
end

drawnow; % Ensure plot rendering is complete
disp('..done plotting ISC components.');

%% ===================== Save and Export Plots =====================
if opts.print == 1
    % Configure figure properties for saving
    set(gcf, 'PaperOrientation', 'landscape');
    set(gcf, 'PaperPosition', [1 1 28 19]);
    cd(opts.pathout);

    % Save plot as PDF
    exportgraphics(gcf, ['ISC_Components_' opts.EOG '_' date '.pdf'], 'ContentType', 'vector', 'Resolution', 600);
    disp('..PDF saved.');
end

% Export additional copy
exportgraphics(gcf, 'ISC_Components_topplot_100.pdf', 'ContentType', 'vector', 'Resolution', 600);
disp('..Export to clipboard done.');

%% ========================= Final Cleanup =========================
clearvars -except opts;
disp('All tasks completed. Results saved in output directory.');

