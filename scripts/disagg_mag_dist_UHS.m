%% Disaggregation vs Observed Comparison (Magnitude and Distance)
% Compares observed and disaggregation-based event parameters at exceedance.

clear; clc; % close all

% Set up directories
data_dir = fullfile('..', 'data');
results_dir = fullfile('..', 'results', 'figures');
if ~exist(results_dir, 'dir'), mkdir(results_dir); end

%% Inputs

% Load CSV with disaggregation results
dissag_selected_Path = fullfile(data_dir, 'RSN_Exceed_MultiPeriod_MCE_with_Dissag.csv');
if isfile(dissag_selected_Path)
    dissag_table = readtable(dissag_selected_Path, 'NumHeaderLines', 1);
    disp('CSV file loaded successfully.');
else
    error('The specified file does not exist. Please check the data/ folder.');
end

% RSN values
rsn_vals = dissag_table.Var1;
% 1) Real Magnitudes
realMag = dissag_table.Var5;
% 1.1) Real Distances
realDist = dissag_table.Var9;
% 2) Disaggregation Magnitudes for T=0.1, T=1.0, T=5.0
mag_1 = dissag_table.Var13;
mag_2 = dissag_table.Var16;
mag_3 = dissag_table.Var19;
disMag = [mag_1 mag_2 mag_3];
% 2.1) Disaggregation Distances for T=0.1, T=1.0, T=5.0
dist_1 = dissag_table.Var14;
dist_2 = dissag_table.Var17;
dist_3 = dissag_table.Var20;
disDist = [dist_1 dist_2 dist_3];
% 3) Exceedance Period Ranges [exceedMin, exceedMax]
exceedMin = dissag_table.Var11;
exceedMax = dissag_table.Var12;
exceedMid = 0.5 * (exceedMin + exceedMax);  % Average exceedance period
% 4) Disaggregation Periods used: T=0.1, T=1.0, T=5.0
Tvals = [0.1, 1.0, 5.0];
numT = length(Tvals);
fadeScale = [1.5, 1.0, 0.4];  % Smaller = fades quicker, larger = fades slower

% Compute axis limits rounded to nearest 0.5
xMin = min(realMag); 
xMax = max(realMag);
yMin = min(disMag(:));
yMax = max(disMag(:));

axisMax = ceil(max(xMax, yMax) * 2) / 2;  % Round up to nearest 0.5
axisMin = floor(min(xMin, yMin) * 2) / 2;  % Round down to nearest 0.5
axisLimits = [axisMin, axisMax];  % Keep aspect ratio

% Define markers, colors, and linewidth
lineWidth = 1.5;
% Define custom colormap (Red → White)
numColors = 100;
colorMatrix = [linspace(1,1,numColors)', linspace(1,0,numColors)', linspace(1,0,numColors)'];

%=== Filtered Plot: Only Distance > 0.6 ===
distance_thrs = 0.6;

%% Magnitude comparison

jitterAmount = 0.03;
colors = [
    0.00, 0.45, 0.70;  % blue
    0.00, 0.60, 0.50;  % teal
    0.85, 0.33, 0.10   % burnt orange
];
markers = {'o', 's', '^'};  % Circle, square, triangle
markerSize = 70;

% Create figure
fig_combined = figure('Name', 'Combined Filtered Points');

set(fig_combined, 'InvertHardcopy', 'off');
set(fig_combined, 'PaperUnits', 'centimeters');
set(fig_combined, 'PaperSize', [15 14]);
set(fig_combined, 'PaperPositionMode', 'manual');
set(fig_combined, 'PaperPosition', [0 0 15 14]);

hold on; grid on;

% For storing all valid filtered points
x_all = [];
y_all = [];

for j = 1:numT
    distance = abs(exceedMid - Tvals(j)) ./ Tvals(j);
    colorVals = exp(-distance / fadeScale(j));
    idx = colorVals > distance_thrs;

    if sum(idx) == 0
        continue;
    end
    rsn_i = rsn_vals(idx);
    x = realMag(idx) + jitterAmount * randn(sum(idx),1);
    y = disMag(idx,j) + jitterAmount * randn(sum(idx),1);

    scatter(x, y, markerSize, markers{j}, ...
            'MarkerFaceColor', colors(j,:), ...
            'MarkerEdgeColor', 'k', ...
            'LineWidth', lineWidth, ...
            'DisplayName', sprintf('$SA(T \\approx %.1f~\\mathrm{s})$', Tvals(j)));

    x_all = [x_all; x];
    y_all = [y_all; y];
end

% Plot y = x line
plot(axisLimits, axisLimits, 'k--', 'LineWidth', 1.2, 'HandleVisibility', 'off');

% Axes formatting
xlim([axisMin, axisMax]);
ylim([axisMin, axisMax]);
xticks(linspace(axisMin, axisMax, 5));
yticks(linspace(axisMin, axisMax, 5));
axis square;

set(gca, 'FontSize', 13, 'LineWidth', 1.0, 'FontName', 'Times New Roman');
xlabel('Observed $M$', 'Interpreter', 'latex', 'FontSize', 15);
ylabel('$\bar{M}$ from disaggregation', 'Interpreter', 'latex', 'FontSize', 15);
lgd = legend('Interpreter', 'latex', 'FontSize', 11, 'Location', 'southeast');


% Export high-quality figure
exportgraphics(fig_combined, fullfile(results_dir, 'combined_dissag_vs_real_magn.pdf'), ...
    'ContentType', 'vector', 'Resolution', 1200);


%% Distance Plot

axisMin = 1;
axisMax = 100;  % Adjust as needed
lineWidth = 1.2;

jitterX = jitterAmount * randn(size(realDist));
jitterY = jitterAmount * randn(size(disDist));

% === Create figure ===
fig_combined_dist = figure('Name', 'Combined Distance Comparison');

set(fig_combined_dist, 'InvertHardcopy', 'off');
set(fig_combined_dist, 'PaperUnits', 'centimeters');
set(fig_combined_dist, 'PaperSize', [15 14]);
set(fig_combined_dist, 'PaperPositionMode', 'manual');
set(fig_combined_dist, 'PaperPosition', [0 0 15 14]);

hold on; grid on;

for j = 1:numT
    distance = abs(exceedMid - Tvals(j)) ./ Tvals(j);
    colorVals = exp(-distance / fadeScale(j));
    idx = colorVals > distance_thrs;

    if sum(idx) == 0
        continue;
    end

    x = realDist(idx) + jitterAmount * randn(sum(idx),1);
    y = disDist(idx,j) + jitterAmount * randn(sum(idx),1);

    scatter(x, y, markerSize, markers{j}, ...
        'MarkerFaceColor', colors(j,:), ...
        'MarkerEdgeColor', 'k', ...
        'LineWidth', lineWidth, ...
        'DisplayName', sprintf('$SA(T \\approx %.1f~\\mathrm{s})$', Tvals(j)));

end

% Diagonal reference line
plot([axisMin, axisMax], [axisMin, axisMax], 'k--', 'LineWidth', 1.2, 'HandleVisibility', 'off');
set(gca, 'XScale', 'log', 'YScale', 'log');  % Set both axes to log scale
% === Format axes ===
xlim([axisMin, axisMax]);
ylim([axisMin, axisMax]);
xticks([1 10 100]);  % Log ticks for base 10
yticks([1 10 100]);
axis square;
set(gca, 'MinorGridLineStyle', ':', 'MinorGridAlpha', 0.05);  % Optional: adjust visibility

xlabel('Observed $R_{\mathrm{rup}}$ [km]', 'Interpreter', 'latex', 'FontSize', 15);
ylabel('$\bar{R}_{\mathrm{rup}}$ from disaggregation [km]', 'Interpreter', 'latex', 'FontSize', 15);
set(gca, 'FontSize', 13, 'LineWidth', 1.0, 'FontName', 'Times New Roman');
legend('Interpreter', 'latex', 'FontSize', 11, 'Location', 'southeast');

% === Export ===
exportgraphics(fig_combined_dist, fullfile(results_dir, 'combined_disagg_vs_real_distances.pdf'), ...
    'ContentType', 'vector', 'Resolution', 1200);