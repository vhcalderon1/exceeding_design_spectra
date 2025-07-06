%% Usable Ground Motions Analysis
% This script computes the number of usable ground motions for each period
% using NGA-West2 data, following the methods described in the associated paper.

clear; clc; % close all

% Set up relative paths
data_dir = fullfile('..', 'data');
results_dir = fullfile('..', 'results', 'figures');
if ~exist(results_dir, 'dir')
    mkdir(results_dir);
end

% Load required data files
mat_file = fullfile(data_dir, 'NGA_W2_corr_meta_data.mat');
csv_file = fullfile(data_dir, 'rsn_set.csv');

if isfile(mat_file) && isfile(csv_file)
    data = load(mat_file); % Loads struct
    rsn_selected = readtable(csv_file);
    disp('Data loaded successfully.');
else
    error('One or more data files are missing. Please check the data/ folder.');
end

% Extract required variables from loaded data
lowestUsableFreq = data.lowest_usable_freq;
magnitude = data.magnitude;
periods = data.Periods;
rsn_global = 1:numel(lowestUsableFreq);
rsn_list = rsn_selected{:,1};

%% Filter and Count Usable Ground Motions

% Filter records with valid lowest usable frequency
zeroIndices = find(lowestUsableFreq <= 0);
validIndices = ~ismember(rsn_global, zeroIndices);
lowestUsableFreq = lowestUsableFreq(validIndices);
magnitude = magnitude(validIndices);
rsn_global(zeroIndices) = [];
maxUsablePeriod = 1 ./ lowestUsableFreq;

% Number of usable ground motions for each period (global)
n_usable_gm = arrayfun(@(T) sum(maxUsablePeriod >= T), periods);

% Among selected records only
[is_selected, idx_selected] = ismember(rsn_global, rsn_list);
index_rsn_global = find(is_selected);
rsn_w_lf = rsn_global(is_selected); 
rsn_w_lf = unique(rsn_w_lf);
lowestUsableFreq_list = lowestUsableFreq(index_rsn_global);
maxUsablePeriod_list = 1 ./ lowestUsableFreq_list;
magnitude_list = magnitude(index_rsn_global);

n_usable_analyzed = arrayfun(@(T) sum(maxUsablePeriod_list >= T), periods);

% Save usable period analysis for further use
periods_usablef = periods;
save(fullfile(data_dir, "usable_gmm.mat"), "periods_usablef", "n_usable_analyzed");

% For magnitude thresholds
threshold_1 = 5.0;
threshold_2 = 7.0;
idx_th_1 = find(magnitude_list >= threshold_1);
idx_th_2 = find(magnitude_list >= threshold_2);
maxUsablePeriod_th_1 = 1 ./ lowestUsableFreq(idx_th_1);
maxUsablePeriod_th_2 = 1 ./ lowestUsableFreq(idx_th_2);

n_usable_th_1 = arrayfun(@(T) sum(maxUsablePeriod_th_1 >= T), periods);
n_usable_th_2 = arrayfun(@(T) sum(maxUsablePeriod_th_2 >= T), periods);

%% Plot Results

colors = lines(3); % For clarity with 3 curves
fig_1 = figure;
hold on;
p1 = plot(periods, n_usable_analyzed, '-', 'Color', colors(1,:), 'LineWidth', 3);
p2 = plot(periods, n_usable_th_1, '--', 'Color', colors(2,:), 'LineWidth', 2);
p3 = plot(periods, n_usable_th_2, '-.', 'Color', colors(3,:), 'LineWidth', 2);

set(gca, 'XScale', 'log', 'YScale', 'log');
xlabel('Period [s]', 'FontSize', 15, 'FontName', 'Times New Roman');
ylabel('Number of Usable Ground Motions', 'FontSize', 15, 'FontName', 'Times New Roman');
lgd = legend([p1, p2, p3], ...
    'All selected recordings', ...
    '$M \geq 5.0$ only', ...
    '$M \geq 7.0$ only', ...
    'Interpreter', 'latex', 'Location', 'northeast');
defaultPos = lgd.Position;
lgd.Position = [0.2, 0.8 - defaultPos(4), defaultPos(3), defaultPos(4)];
set(gca, 'FontSize', 14, 'LineWidth', 1.5, 'FontName', 'Times New Roman');
xlim([0 10]);
ylim([1e2 2e4]);
ax = gca; ax.YRuler.Exponent = 4;
set(fig_1, 'Color', 'w'); box off; legend boxoff; grid off; hold off;
fig_1.Position = [100, 100, 650, 500]; % adjust as needed

% Save the figure as a PDF in results/figures/
exportgraphics(fig_1, fullfile(results_dir, 'usable_gms.pdf'), ...
    'ContentType', 'vector', 'Resolution', 1200, ...
    'BackgroundColor', 'white');
disp("")
