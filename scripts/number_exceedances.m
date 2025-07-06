%% Number of Exceedances Analysis
% Computes and plots the number of exceedances for various design spectra.

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

% Helper function to load local .mat data
load_data = @(fname) load(fullfile(data_dir, fname));

% --- UHS ---
d = load_data('rsn_exceed_NGA_West_UHS.mat');
rsn_great_UHS = d.rsn_great_1;
not_found_RSN_UHS = d.not_found_RSN;
d = load_data('T_vs_RSN_UHS.mat'); T_uni_great_UHS = d.T_uni_great_1; % Periods where exceedance happens
d = load_data('intervals_exceed_UHS.mat'); t_int_UHS = d.t_int_1; % 1st col: RSN; 2nd cold: T_min; 3rd: T_max
% --- DBE ---
d = load_data('rsn_exceed_NGA_West_DBE.mat'); rsn_great_Two_DBE = d.rsn_great_1;
not_found_RSN_DBE = d.not_found_RSN;
d = load_data('T_vs_RSN_DBE.mat'); T_uni_great_Two_DBE = d.T_uni_great_1;
d = load_data('intervals_exceed_DBE.mat'); t_int_Two_DBE = d.t_int_1;
d = load_data('rsn_exceed_NGA_West_DBE.mat'); rsn_great_Multi_DBE = d.rsn_great_2;
d = load_data('T_vs_RSN_DBE.mat'); T_uni_great_Multi_DBE = d.T_uni_great_2;
d = load_data('intervals_exceed_DBE.mat'); t_int_Multi_DBE = d.t_int_2;
% --- MCE ---
d = load_data('rsn_exceed_NGA_West_MCE.mat'); rsn_great_Two_MCE = d.rsn_great_1;
not_found_RSN_MCE = d.not_found_RSN;
d = load_data('T_vs_RSN_MCE.mat'); T_uni_great_Two_MCE = d.T_uni_great_1;
d = load_data('intervals_exceed_MCE.mat'); t_int_Two_MCE = d.t_int_1;
d = load_data('rsn_exceed_NGA_West_MCE.mat'); rsn_great_Multi_MCE = d.rsn_great_2;
d = load_data('T_vs_RSN_MCE.mat'); T_uni_great_Multi_MCE = d.T_uni_great_2;
d = load_data('intervals_exceed_MCE.mat'); t_int_Multi_MCE = d.t_int_2;
% --- RT ---
d = load_data('rsn_exceed_NGA_West_RT.mat'); rsn_great_RT = d.rsn_great_1;
not_found_RSN_RT = d.not_found_RSN;
d = load_data('T_vs_RSN_RT.mat'); T_uni_great_RT = d.T_uni_great_1;
d = load_data('intervals_exceed_RT.mat'); t_int_RT = d.t_int_1;

%% Extract Earthquake Information
% Use helper function from src
[EQ_comp_1, ~, ~, ~, ~, ~, ~, ~] = extractEQData(rsn_great_UHS, rsn_gm, EQ_name, EQ_year, station_name, magnitude, closest_D, soil_Vs30, latitude, longitude);
[EQ_comp_2, ~, ~, ~, ~, ~, ~, ~] = extractEQData(rsn_great_Two_DBE, rsn_gm, EQ_name, EQ_year, station_name, magnitude, closest_D, soil_Vs30, latitude, longitude);
[EQ_comp_3, ~, ~, ~, ~, ~, ~, ~] = extractEQData(rsn_great_Multi_DBE, rsn_gm, EQ_name, EQ_year, station_name, magnitude, closest_D, soil_Vs30, latitude, longitude);
[EQ_comp_4, ~, ~, ~, ~, ~, ~, ~] = extractEQData(rsn_great_Two_MCE, rsn_gm, EQ_name, EQ_year, station_name, magnitude, closest_D, soil_Vs30, latitude, longitude);
[EQ_comp_5, ~, ~, ~, ~, ~, ~, ~] = extractEQData(rsn_great_Multi_MCE, rsn_gm, EQ_name, EQ_year, station_name, magnitude, closest_D, soil_Vs30, latitude, longitude);
[EQ_comp_6, ~, ~, ~, ~, ~, ~, ~] = extractEQData(rsn_great_RT, rsn_gm, EQ_name, EQ_year, station_name, magnitude, closest_D, soil_Vs30, latitude, longitude);

%% Count exceedances using helper from src
[ycom_1, ~, ~, limit_y_1] = exceed_Periods(T_uni_great_Two_DBE, Periods);
[ycom_2, ~, ~, limit_y_2] = exceed_Periods(T_uni_great_Multi_DBE, Periods);
[ycom_3, ~, ~, limit_y_3] = exceed_Periods(T_uni_great_Two_MCE, Periods);
[ycom_4, ~, ~, limit_y_4] = exceed_Periods(T_uni_great_Multi_MCE, Periods);
[ycom_5, ~, ~, limit_y_5] = exceed_Periods(T_uni_great_UHS, Periods);
[ycom_6, ~, ~, limit_y_6] = exceed_Periods(T_uni_great_RT, Periods);

lim_y_max = max([limit_y_1, limit_y_2, limit_y_3, limit_y_4, limit_y_5, limit_y_6]);

colors = lines(7);

%% Plot: Number of Exceedances for All Periods

fig_1 = figure;
hold on;
p_1 = plot(Periods, ycom_1, ":", 'Color', colors(1,:), 'LineWidth', 3);
p_2 = plot(Periods, ycom_2, "-", 'Color', colors(1,:), 'LineWidth', 2);
p_3 = plot(Periods, ycom_3, ":", 'Color', colors(2,:), 'LineWidth', 3);
p_4 = plot(Periods, ycom_4, "-",'Color', colors(2,:), 'LineWidth', 2);
p_5 = plot(Periods, ycom_5, 'Color', colors(4,:), 'LineWidth', 2);
p_6 = plot(Periods, ycom_6, 'Color', colors(7,:), 'LineWidth', 2);
xlabel('Period [s]', 'FontSize', 15, 'FontName', 'Times New Roman');
ylabel('Number of Exceedances','FontSize',15,'FontName','Times New Roman');
legend([p_1 p_2 p_3 p_4 p_6 p_5], ...
    'Two-Period $ \mathrm{DE} $ \quad', ...
    'Multi-Period $ \mathrm{DE} $ \quad', ...
    'Two-Period $\mathrm{MCE}_{\mathrm{R}}$ \quad', ...
    'Multi-Period $\mathrm{MCE}_{\mathrm{R}}$ \quad', ...
    'Probabilistic $\mathrm{MCE}_{\mathrm{R}}$ \quad', ...
    'UHS RP=2475 yrs \quad', ...
    'Interpreter', 'latex');
set(gca, 'FontSize', 14, 'LineWidth', 1.5, 'FontName', 'Times New Roman');
xlim([0 10]);
ylim([0 80]);
yticks(0.0:10:80);
set(gca, 'XScale', 'log')
set(fig_1, 'Color', 'w');
box off;
legend boxoff
hold off;
exportgraphics(fig_1, fullfile(results_dir, 'number_exceedances_all.pdf'), ...
    'ContentType', 'vector', 'Resolution', 1200, 'BackgroundColor', 'white');

%% Figure 1 - Multi period only

fig_1_M = figure;
hold on;
p_2 = plot(Periods, ycom_2, "-", 'Color', colors(1,:), 'LineWidth', 2);
p_4 = plot(Periods, ycom_4, "-",'Color', colors(2,:), 'LineWidth', 2);
p_5 = plot(Periods, ycom_5, 'Color', colors(4,:), 'LineWidth', 2);
p_6 = plot(Periods, ycom_6, 'Color', colors(7,:), 'LineWidth', 2);
xlabel('Period [s]', 'FontSize', 15, 'FontName', 'Times New Roman');
ylabel('Number of Exceedances','FontSize',15,'FontName','Times New Roman');
legend([p_2 p_4 p_5 p_6], ...
    '$\mathrm{DE}$ \quad', ...
    '$\mathrm{MCE}_{\mathrm{R}}$', ...
    '$\mathrm{UHS}$ RP=2475 yrs \quad', ...
    'Probabilistic $\mathrm{MCE}_{\mathrm{R}}$ \quad', ...
    'Interpreter', 'latex');
set(gca, 'FontSize', 14, 'LineWidth', 1.5, 'FontName', 'Times New Roman');
xlim([0 10]);
ylim([0 70]);
yticks(0.0:10:70);
set(gca, 'XScale', 'log')
set(fig_1_M, 'Color', 'w');
box off;
legend boxoff
hold off;
exportgraphics(fig_1_M, fullfile(results_dir, 'number_exceedances_multi.pdf'), ...
    'ContentType', 'vector', 'Resolution', 1200, 'BackgroundColor', 'white');