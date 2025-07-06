%% Ratio SA/UHS vs T/Tmax Analysis
% Analyzes and plots the ratio of SA/UHS as a function of normalized period T/Tmax.

clear; clc; % close all

% Add src folder to path for helper functions (if needed)
addpath(fullfile('..', 'src'));

% Set up directories
data_dir = fullfile('..', 'data');
results_dir = fullfile('..', 'results', 'figures');
if ~exist(results_dir, 'dir')
    mkdir(results_dir);
end

%% Inputs

% Load RSN numbers and locations
rsn_table = readtable(fullfile(data_dir, 'rsn_set.csv'));
rsn_gm = rsn_table.RecordSequenceNumber;
latitude = rsn_table.StationLatitude;
longitude = rsn_table.StationLongitude;
magnitude_unit = rsn_table.EarthquakeMagnitude;

% RotD100
matFileName = 'NGA_W2_corr_meta_data.mat';

% NGA-West2 metadata
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
T = Periods;

% Results storage
resultFile = fullfile(data_dir, 'UHS_Results.mat');
if ~isfile(resultFile)
    save(resultFile, 'rsn_gm');
end
storedResults = matfile(resultFile, 'Writable', true);

% Load exceeding records
load(fullfile(data_dir, 'rsn_exceed_NGA_West_UHS.mat'));
load(fullfile(data_dir, 'T_vs_RSN_UHS.mat'));
load(fullfile(data_dir, 'intervals_exceed_UHS.mat'));

%% Analysis

if ~isfile(fullfile(data_dir, "SA_T_max_UHS.mat"))
    % Initialize output arrays
    T_max = zeros(length(rsn_great_2), 1);
    ratio_dif_SA_1 = zeros(length(rsn_great_2), 105);

    % Loop through each RSN ID
    for i = 1:length(rsn_great_2)
        ratio_max = 0;
        T_max_i = 0;
        rsn_val = rsn_great_2(i);
        SA_target_1 = storedResults.(sprintf('SA_UHS_%d', rsn_val));
        T_target_1  = storedResults.(sprintf('T_UHS_%d', rsn_val));
        Sa_RotD100  = Sa_RotD100_Global(rsn_val, :);
        % Choose appropriate T array
        if max(T) <= max(T_target_1)
            T_a1 = T;
            disp("Caramba")
        else
            T_a1 = T_target_1;
        end

        index_T = find(T == max(T_a1));

        % Loop over periods up to index_T
        for j = 1:index_T
            index_T1 = find(T_a1 == T(j), 1);
            if isempty(index_T1)
                SA_1 = linear_interpol(T(j), T_target_1, SA_target_1);
            else
                SA_1 = SA_target_1(index_T1);
            end

            ratio_dif_SA_1(i,j) = Sa_RotD100(j) / SA_1;

            if j >= 2 && ratio_dif_SA_1(i,j) >= ratio_max
                T_max_i = T(j);
                ratio_max = ratio_dif_SA_1(i,j);
            end
        end
        T_max(i) = T_max_i;
    end

    % Save results
    save(fullfile(data_dir, "SA_T_max_UHS.mat"), "T_max", "ratio_dif_SA_1");
    disp("Saved SA_T_max_UHS.mat");
else
    load(fullfile(data_dir, "SA_T_max_UHS.mat"))
end

%% Interpolation of values

% Obtaining T/Tmax
T_total = T(1:size(ratio_dif_SA_1,2));
T_ratio = (1./T_max)*T_total;
Tratio_min = min(T_ratio,[],"all");
Tratio_max = max(T_ratio,[],"all");
delta_ratio = 0.001;
Tratio_int = Tratio_min:delta_ratio:Tratio_max;
Tratio_min_mat = min(T_ratio,[],2);
Tratio_max_mat = max(T_ratio,[],2);
ratio_intp = nan(length(rsn_great_2), length(Tratio_int));

for i = 1:length(rsn_great_2)
    [~, idx_min] = min(abs(Tratio_int - Tratio_min_mat(i)));
    clst_min = Tratio_int(idx_min);
    [~, idx_max] = min(abs(Tratio_int - Tratio_max_mat(i)));
    clst_max = Tratio_int(idx_max);
    Tratio_int_i =  Tratio_int(idx_min:idx_max);
    ratio_intp(i, idx_min:idx_max) = interp1(T_ratio(i,:), ratio_dif_SA_1(i,:), Tratio_int_i, 'linear', NaN);
end
ratio_mean = mean(ratio_intp, 1, 'omitnan');
georatio_mean = geomean(ratio_intp, 1, 'omitnan');

%% Plot

[max_val, linear_idx] = max(ratio_dif_SA_1(:));
[row, col] = ind2sub(size(ratio_dif_SA_1), linear_idx);

fig = figure;
hold on;
yline(1, '--', 'Color', [1, 102/255, 0], 'LineWidth', 2, 'HandleVisibility', 'off');
index_T = 105;
rows_exceed_1_5 = [];
rows_exceed_2_5 = [];

color_gray = [0.7 0.7 0.7];      % Light gray for background lines
color_median = [56 95 150]/255;    % Blue for median

for i = 1:size(ratio_dif_SA_1, 1)
    % Normalize T by T_max(i)
    rsn_val = rsn_great_2(i);
    T_norm = T(1:index_T) ./ T_max(i);
    y_vals = ratio_dif_SA_1(i, 1:index_T);

    % Save threshold exceedance info
    if any(y_vals > 1.5)
        rows_exceed_1_5(end+1) = i;
    end
    if any(y_vals > 2.5)
        rows_exceed_2_5(end+1) = i;
    end
    if i == 1
        h_random = plot(T_norm, y_vals, '-', 'Color', color_gray, 'LineWidth', 1.5, ...
            'DisplayName', 'Random Record');
    else
        % Plot individual line in light gray, hide from legend
        plot(T_norm, y_vals, '-', 'Color', color_gray, 'LineWidth', 1, ...
            'HandleVisibility', 'off');
    end
end

color_max = [0.55 0.55 0.55];
T_norm_max = T(1:index_T) ./ T_max(row);
y_vals_max = ratio_dif_SA_1(row, 1:index_T);
plot(T_norm_max, y_vals_max, '-', ...
    'Color', color_max, 'LineWidth', 1.5, ...
    'HandleVisibility', 'off');

% Plot the median (or mean) as a thick blue line and show in legend
h_median = plot(Tratio_int, georatio_mean, '-', ...
    'Color', color_median, 'LineWidth', 2.5, ...
    'DisplayName', 'Median');
set(gca, 'XScale', 'log', ...
    'FontSize', 14, 'LineWidth', 1.5, ...
    'FontName', 'Times New Roman');

% Step 1: Set xlim to be slightly wider than your first/last ticks
xlim([3e-3, 2.0e2]);  % Widen the limits a bit beyond 1e-3 and 1e3

% Step 2: Set ticks and labels
xticks([3e-3 1e-2 1e-1 1e0 1e1 1e2 2e2]);
xticklabels({'','10^{-2}','10^{-1}','10^{0}','10^{1}','10^{2}',''});
set(gca, 'TickLabelInterpreter', 'tex', 'FontName', 'Times New Roman'); % TeX for Times

% set(gca, 'FontName', 'Times New Roman')
ylim([0, 4]);  % Lower limit set to 1
yticks(0:1:4);

xlabel('$T/T_{\max}$', 'Interpreter', 'latex', 'FontSize', 15);
ylabel('$\mathit{SA} / \mathit{UHS}$', 'Interpreter', 'latex', 'FontSize', 15);

box off;
% set(gca, 'FontSize', 11, 'LineWidth', 1.2, 'TickLabelInterpreter', 'latex');
set(gca, 'Layer', 'top');
legend('show')
legend([h_random, h_median], {'Individual recordings', 'Median'}, ...
    'Location','best','FontSize',13,'Interpreter','latex','Box','off');

% Save figure
set(fig, 'Units', 'centimeters', 'Position', [2 2 15 12]);
exportgraphics(fig, fullfile(results_dir, 'SA_Ratio_Normalized_with_Mean_UHS.pdf'), ...
    'ContentType', 'vector', 'Resolution', 1200, 'BackgroundColor', 'white');