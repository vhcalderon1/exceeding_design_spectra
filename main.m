% SPDX-FileCopyrightText: 2025 Stanford University
%
% SPDX-License-Identifier: MIT
%
% Main script to call the principal analysis codes for this repository.
%
% Many scripts require access to supporting data in the /data and /src folders.
% Make sure you are running this script from the project root directory.
%
% Results, methodology, and motivations are documented in:
%
% Calderon and Baker (2025), "Observed ground motions that exceeded design response 
% spectra in the Western United States", submitted to Earthquake Spectra (in review).
%
% Project repository: https://github.com/vhcalderon1/exceeding_design_spectra

clear; close all; clc;

disp('Starting main analysis for exceeding_design_spectra...');

% -------------------------------------------------------------------------
% Add subfolders for script and helper function access
% -------------------------------------------------------------------------
addpath(fullfile(pwd, 'scripts'));
addpath(fullfile(pwd, 'src'));

% -------------------------------------------------------------------------
% Each section calls a key analysis or figure script.
% Individual scripts read/write from /data and /results as needed.
% -------------------------------------------------------------------------

%% 1. Disaggregation of UHS residuals
disp('Running disagg_residual_UHS.m...');
run(fullfile('scripts', 'disagg_residual_UHS.m'));

%% 2. Disaggregation-based Magnitude & Distance (UHS)
disp('Running disagg_mag_dist_UHS.m...');
run(fullfile('scripts', 'disagg_mag_dist_UHS.m'));

%% 3. Ratio of SA/UHS vs T/Tmax analysis
disp('Running ratiosa_vs_tmax.m...');
run(fullfile('scripts', 'ratiosa_vs_tmax.m'));

%% 4. Magnitude-Distance plotting and analysis
disp('Running magnitude_distance.m...');
run(fullfile('scripts', 'magnitude_distance.m'));

%% 5. Number of UHS exceedances calculation
disp('Running number_exceedances.m...');
run(fullfile('scripts', 'number_exceedances.m'));

%% 6. Set of RSNs and plotting
disp('Running set_rsns_plot.m...');
run(fullfile('scripts', 'set_rsns_plot.m'));

%% 7. Usable GMs for analysis
disp('Running usable_gms.m...');
run(fullfile('scripts', 'usable_gms.m'));

disp('All analyses completed.');
