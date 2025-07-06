%% Magnitude vs. Distance Analysis for Exceedance Thresholds
% Plots magnitude vs. distance for groups exceeding period thresholds.

clear; clc; % close all

% Add src folder to MATLAB path for helper functions
addpath(fullfile('..', 'src'));

% Define local data and results directories
data_dir = fullfile('..', 'data');
results_dir = fullfile('..', 'results', 'figures');
if ~exist(results_dir, 'dir')
    mkdir(results_dir);
end

%% Load Inputs

% Load RSN numbers and locations
rsn_table = readtable(fullfile(data_dir, 'rsn_set.csv'));
rsn_gm = rsn_table.RecordSequenceNumber;
latitude = rsn_table.StationLatitude;
longitude = rsn_table.StationLongitude;
magnitude_unit = rsn_table.EarthquakeMagnitude;

% Load NGA-West2 metadata
matFileName = fullfile(data_dir, 'NGA_W2_corr_meta_data.mat');
variablesToLoad = {'Sa_RotD100', 'Periods', 'soil_Vs30','station_name','EQ_year','EQ_name','closest_D','magnitude'};
loadRotD = load(matFileName, variablesToLoad{:});
Sa_RotD100_Global = loadRotD.Sa_RotD100;
Periods = loadRotD.Periods;
soil_Vs30 = loadRotD.soil_Vs30;
station_name = loadRotD.station_name;
EQ_year = loadRotD.EQ_year;
EQ_name = loadRotD.EQ_name;
closest_D = loadRotD.closest_D;
magnitude = loadRotD.magnitude;

% Helper function for local .mat loading
load_data = @(fname) load(fullfile(data_dir, fname));

% --- UHS ---
d = load_data('rsn_exceed_NGA_West_UHS.mat');
rsn_great_UHS = d.rsn_great_1;
not_found_RSN_UHS = d.not_found_RSN;
d = load_data('T_vs_RSN_UHS.mat'); T_uni_great_UHS = d.T_uni_great_1;
d = load_data('intervals_exceed_UHS.mat'); t_int_UHS = d.t_int_1;

% --- MCE ---
d = load_data('rsn_exceed_NGA_West_MCE.mat'); rsn_great_Two_MCE = d.rsn_great_1; rsn_great_Multi_MCE = d.rsn_great_2;
not_found_RSN_MCE = d.not_found_RSN;
d = load_data('T_vs_RSN_MCE.mat'); T_uni_great_Two_MCE = d.T_uni_great_1; T_uni_great_Multi_MCE = d.T_uni_great_2; % Periods where exceedance happens
d = load_data('intervals_exceed_MCE.mat'); t_int_Two_MCE = d.t_int_1; t_int_Multi_MCE = d.t_int_2; % 1st col: RSN; 2nd cold: T_min; 3rd:  T_max

% --- RT ---
d = load_data('rsn_exceed_NGA_West_RT.mat'); rsn_great_RT = d.rsn_great_1;
not_found_RSN_RT = d.not_found_RSN;
d = load_data('T_vs_RSN_RT.mat'); T_uni_great_RT = d.T_uni_great_1; % Periods where exceedance happens
d = load_data('intervals_exceed_RT.mat'); t_int_RT = d.t_int_1; % 1st col: RSN; 2nd cold: T_min; 3rd:  T_max

%% Group RSNs by period exceedance threshold

T_excee = 1.0; % Separate per period of exceedance

distance_global = zeros(length(rsn_gm)-1,1);
magnitude_global = zeros(length(rsn_gm)-1,1);
for i = 1: length(rsn_gm)-1
    distance_global(i) = closest_D(rsn_gm(i));
    magnitude_global(i) = magnitude(rsn_gm(i));
end

rsn_thres_Multi_MCE_sup = [];
rsn_thres_Multi_MCE_inf = [];
for i = 1:size(t_int_Multi_MCE,1)
    if t_int_Multi_MCE(i,3) > T_excee
        rsn_thres_Multi_MCE_sup = [rsn_thres_Multi_MCE_sup; t_int_Multi_MCE(i,1)];
    else
        rsn_thres_Multi_MCE_inf = [rsn_thres_Multi_MCE_inf; t_int_Multi_MCE(i,1)];
    end
end

rsn_thres_RT_sup = [];
rsn_thres_RT_inf = [];
for i = 1:size(t_int_RT,1)
    if t_int_RT(i,3) > T_excee
        rsn_thres_RT_sup = [rsn_thres_RT_sup; t_int_RT(i,1)];
    else
        rsn_thres_RT_inf = [rsn_thres_RT_inf; t_int_RT(i,1)];
    end
end

% Extracting data
[~, ~, ~, magnitude_comp_7, distance_comp_7, ~,~,~] = extractEQData(rsn_thres_Multi_MCE_sup, rsn_gm, EQ_name, EQ_year, station_name, magnitude, closest_D, soil_Vs30, latitude, longitude);
[~, ~, ~, magnitude_comp_8, distance_comp_8, ~,~,~] = extractEQData(rsn_thres_Multi_MCE_inf, rsn_gm, EQ_name, EQ_year, station_name, magnitude, closest_D, soil_Vs30, latitude, longitude);
[~, ~, ~, magnitude_comp_11, distance_comp_11, ~,~,~] = extractEQData(rsn_thres_RT_sup, rsn_gm, EQ_name, EQ_year, station_name, magnitude, closest_D, soil_Vs30, latitude, longitude);
[~, ~, ~, magnitude_comp_12, distance_comp_12, ~,~,~] = extractEQData(rsn_thres_RT_inf, rsn_gm, EQ_name, EQ_year, station_name, magnitude, closest_D, soil_Vs30, latitude, longitude);

%% Plotting: Magnitude vs. Distance

color_group7 = [0.2, 0.5, 0.9];    % Multi-Period MCE
color_group11 = [1, 0, 0];         % Probabilistic MCE
color_global = [0.75, 0.75, 0.75]; % All selected
size_global = 5;
size_group7 = 90;
size_group11 = 50;
group_glob_label = 'All selected recordings';

% Define sets for scatter plots
global_set = [distance_global, magnitude_global];
group_7 = [distance_comp_7, magnitude_comp_7];
group_8 = [distance_comp_8, magnitude_comp_8];
group_11 = [distance_comp_11, magnitude_comp_11];
group_12 = [distance_comp_12, magnitude_comp_12];

y_lim_sup = 7.5; y_lim_inf = 3.0; y_lim_int = 1.5;

% --------- FIRST TILE: Below threshold ---------

f = figure;
set(f, 'Units', 'pixels', 'Position', [100, 100, 1100, 500]);

t = tiledlayout(1,2, 'TileSpacing','compact', 'Padding','compact');
ax1 = nexttile; hold(ax1, 'on');

% Plot global set
scatter(ax1, global_set(:,1), global_set(:,2), size_global, 'filled', ...
    'MarkerFaceColor', color_global, 'DisplayName', group_glob_label);

% Plot group 8
scatter(ax1, group_8(:,1), group_8(:,2), size_group7, '>', ...
    'MarkerEdgeColor', color_group7, 'LineWidth', 1.5, ...
    'DisplayName', '$\mathrm{MCE}_{\mathrm{R}}$');

% Plot group 12
scatter(ax1, group_12(:,1), group_12(:,2), size_group11, 'o', ...
    'MarkerEdgeColor', color_group11, 'LineWidth', 1.5, ...
    'DisplayName', 'Probabilistic $\mathrm{MCE}_{\mathrm{R}}$');

% Axis formatting
set(ax1, 'XScale', 'log');
xlim(ax1, [1e-2 1e4]);
ylim(ax1, [y_lim_inf, y_lim_sup]);
yticks(ax1, y_lim_inf:y_lim_int:y_lim_sup);
xlabel(ax1, 'Source-to-site distance, $R_{\mathrm{rup}}$ [km]', ...
    'Interpreter', 'latex', 'FontSize', 15);
ylabel(ax1, 'Magnitude', 'Interpreter', 'latex', 'FontSize', 15);
title(ax1, ['Exceedance Below $T \leq $' num2str(T_excee)], ...
    'Interpreter', 'latex', 'FontSize', 15);
legend(ax1, 'Location', 'southwest', 'FontSize', 12, 'Interpreter', 'latex');
legend(ax1, 'boxoff');
set(ax1, 'FontSize', 14, 'LineWidth', 1.5, 'TickLabelInterpreter', 'latex', 'FontName', 'Times New Roman');
box(ax1, 'on');
hold(ax1, 'off');

% --------- SECOND TILE: Above threshold ---------
ax2 = nexttile; hold(ax2, 'on');

% Plot global set
scatter(ax2, global_set(:,1), global_set(:,2), size_global, 'filled', ...
    'MarkerFaceColor', color_global, 'DisplayName', group_glob_label);

% Plot group 7
scatter(ax2, group_7(:,1), group_7(:,2), size_group7, '>', ...
    'MarkerEdgeColor', color_group7, 'LineWidth', 1.5, ...
    'DisplayName', '$\mathrm{MCE}_{\mathrm{R}}$');

% Plot group 11
scatter(ax2, group_11(:,1), group_11(:,2), size_group11, 'o', ...
    'MarkerEdgeColor', color_group11, 'LineWidth', 1.5, ...
    'DisplayName', 'Probabilistic $\mathrm{MCE}_{\mathrm{R}}$');

% Axis formatting
set(ax2, 'XScale', 'log');
xlim(ax2, [1e-2 1e4]);
ylim(ax2, [y_lim_inf, y_lim_sup]);
yticks(ax2, y_lim_inf:y_lim_int:y_lim_sup);
set(ax2, 'YTickLabel', {});  % Hide the y labels for the right plot
xlabel(ax2, 'Source-to-site distance, $R_{\mathrm{rup}}$ [km]', ...
    'Interpreter', 'latex', 'FontSize', 15);
title(ax2, ['Exceedance Above $T > $' num2str(T_excee)], ...
    'Interpreter', 'latex', 'FontSize', 15);
legend(ax2, 'Location', 'southwest', 'FontSize', 12, 'Interpreter', 'latex');
legend(ax2, 'boxoff');
set(ax2, 'FontSize', 14, 'LineWidth', 1.5, 'TickLabelInterpreter', 'latex', 'FontName', 'Times New Roman');
box(ax2, 'on');
hold(ax2, 'off');

% Set white background
set(gcf, 'Color', 'w');

% Export
exportgraphics(f, fullfile(results_dir, 'magnitude_distance_subplot.pdf'), ...
    'ContentType', 'vector', 'Resolution', 1200, 'BackgroundColor', 'white');