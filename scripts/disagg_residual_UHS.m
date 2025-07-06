%% Disaggregation-based Magnitude & Distance (UHS cases)
% Computes residuals and compares observed vs. disaggregation epsilon values.

clear; clc; % close all

% Set up directories
addpath(fullfile('..','src'));              % helper functions
data_dir    = fullfile('..','data');
results_dir = fullfile('..','results','figures');
if ~exist(results_dir,'dir'), mkdir(results_dir); end

%% load raw data
load(fullfile(data_dir,'NGA_W2_corr_meta_data.mat'));   % creates Periods, Sa_RotD50, Sa_RotD100, …
load(fullfile(data_dir,'Z1_Z25_updated.mat'));          % Z1 variables
Z1_CVMH(Z1_CVMH>0) = Z1_CVMH(Z1_CVMH>0)/1000;           % convert m → km

% CSV MCE Exceedances and Dissagregation
csvPath = fullfile(data_dir,'RSN_Exceed_UHS_with_Dissag.csv');
if isfile(csvPath)
    dataTable_1   = readtable(csvPath);
    disp('UHS disaggregation CSV loaded successfully.');
else
    error('RSN_Exceed_UHS_with_Dissag.csv not found in data/.');
end
rsn_vals = dataTable_1.RSN;
epsilon_T_0P1 = dataTable_1.mean_epsiolon_T_0P1;
epsilon_T_1P0 = dataTable_1.mean_epsiolon_T_1P0;
epsilon_T_5P0 = dataTable_1.mean_epsiolon_T_5P0;
epsilon_dis = [epsilon_T_0P1 epsilon_T_1P0 epsilon_T_5P0];
load(fullfile(data_dir,'plot_dissagre_UHS_val.mat'));   % supplies idx_vals, etc.
exceedMin  = dataTable_1.T__min_;      exceedMax = dataTable_1.T__max_;
exceedMid  = 0.5*(exceedMin+exceedMax);

%% initialize matrices

% trim long periods so periods match the GMPE range
tIndex = find(Periods<=10); 
Periods = Periods(tIndex);
Sa_RotD50 = Sa_RotD50(:,tIndex);
Sa_RotD100 = Sa_RotD100(:,tIndex);

% get sizes of matrices
% numRecs = length(eqid);
numRecs = length(rsn_vals);
numT = length(Periods);

[eventEQID,eventIdx] = unique(eqid); % get indices of unique events
eventEQID(1) = []; % throw out the first case, which is eqid = -999
eventIdx(1) = []; % throw out the first case, which is eqid = -999
for i = 1:length(eventEQID)
   idx = find(eqid == eventEQID(i),1);
   eventMag(i,1) = magnitude(idx); % get event magnitudes 
end

% initialize residual matrices
% === Residual matrices for CY14 ===
resid_CY14_RotD100Within        = nan(numRecs, numT);
resid_CY14_RotD100Total         = nan(numRecs, numT);
resid_CY14_RotD100BetweenLong   = nan(max(eqid), numT);

% === Residual matrices for ASK14 ===
resid_ASK14_RotD100Within       = nan(numRecs, numT);
resid_ASK14_RotD100Total        = nan(numRecs, numT);
resid_ASK14_RotD100BetweenLong  = nan(max(eqid), numT);

resid_BSSA_RotD100Within = nan(numRecs, numT);
resid_BSSA_RotD100BetweenLong = nan(max(eqid), numT);
resid_BSSA_RotD100Total = nan(numRecs, numT);

resid_CB14_RotD100Within       = nan(numRecs, numT);
resid_CB14_RotD100BetweenLong  = nan(max(eqid), numT);
resid_CB14_RotD100Total        = nan(numRecs, numT);
medianPred_CB14                = nan(numRecs, numT);

% === Additional book-keeping ===
resid_BetweenNumRecs            = nan(max(eqid), 1);
medianPred_CY14                 = nan(numRecs, numT);
medianPred_ASK14                = nan(numRecs, numT);
numUsableRecs                   = nan(numT, 1);


%% get predictions

% get fw/hw terms for CY 2014 GMPE
FwHw = zeros(length(Sa_RotD100),1);
for i = 1:length(Sa_RotD100)
    if strcmp(Fhw{i},'hw')
        FwHw(i) = 1;
    end
end

% Predicted RotD100/RotD50 ratios
[ rotD100Factor, sigmaRotD100Factor , phiRotD100Factor, tauRotD100Factor] = SB_2014_ratios( Periods );

% get max usable period for each ground motion
maxUsableT = 1./lowest_usable_freq;

% save CY2014 median and rotd100 spectral values
median_gmm = zeros(length(rsn_vals), numT);
rotD100_gmm = zeros(length(rsn_vals), numT);

%% Sa residuals
for i = 1:numT % period index
    % i
    disp(Periods(i))

    % allowableFilter = (maxUsableT > Periods(i) & Sa_RotD100(:,i)>0);
    % idxTemp = find(allowableFilter & chiouYoungsUsed); % filter records
    % [~, eqIdx] = sort(eqid(idxTemp)); % sort by event
    % idx = idxTemp(eqIdx);
    idx = rsn_vals;

    numUsableRecs(i) = length(idx);

    %% ===== CY14 MODEL: Median & Sigmas =====
    for k = 1:length(idx)
        % idx(k)
        [cy_median(k,1), cy_sigma(k,1), cy_phi(k,1), cy_tau(k,1)] = CY_2014_nga_mod( ...
            magnitude(idx(k)), Periods(i), closest_D(idx(k)), Rjb(idx(k)), Rx(idx(k)), ...
            Z_tor(idx(k)), dip(idx(k)), rakeAngle(idx(k)), Z1_CVMH(idx(k)), ...
            soil_Vs30(idx(k)), FwHw(idx(k)), 1, Region_BSSA(idx(k)));
        disp(['Inputs for rsn = ' num2str(idx(k)) ':']);
        fprintf('Inputs for k = %d: mag=%.3f, Width=%.3f, T=%.3f, Rrup=%.3f, Rjb=%.3f, Rx=%.3f, Ztor=%.3f, dip=%.3f, rake=%.3f, Z1=%.3f, Vs30=%.3f, FwHw=%d, Er=%d\n', ...
            k, ...
            magnitude(idx(k)), ...
            ruptureWidth(idx(k)),...
            Periods(i), ...
            closest_D(idx(k)), ...
            Rjb(idx(k)), ...
            Rx(idx(k)), ...
            Z_tor(idx(k)), ...
            dip(idx(k)), ...
            rakeAngle(idx(k)), ...
            Z1_CVMH(idx(k)), ...
            soil_Vs30(idx(k)), ...
            FwHw(idx(k)), ...
            Region_BSSA(idx(k)) ...
        );

    end

    % ===== Apply RotD100 correction to CY14 =====
    inputCY.imObservations = Sa_RotD100(idx,i);
    inputCY.imMedian = cy_median .* rotD100Factor(i);
    median_gmm(:,i) = cy_median;
    rotD100_gmm(:,i) = cy_median .* rotD100Factor(i);
    inputCY.totalSigma = sqrt(cy_sigma.^2 + sigmaRotD100Factor(i).^2);
    inputCY.interSigma = sqrt(cy_tau.^2 + tauRotD100Factor(i).^2);
    inputCY.intraSigma = sqrt(cy_phi.^2 + phiRotD100Factor(i).^2);
    inputCY.eventIds = eqid(idx);
    outCY = GetInterIntraEventResiduals(inputCY);

    resid_CY14_RotD100Within(:,i) = outCY.intraEventResidualsNormalized;
    resid_CY14_RotD100Total(:,i) = (log(Sa_RotD100(idx,i)) - log(inputCY.imMedian)) ./ inputCY.totalSigma;

    idxNgms = find(outCY.eventData.eventNumGms <= 1);
    outCY.eventData.interEventResidualNormalized(idxNgms) = nan;
    resid_CY14_RotD100BetweenLong(outCY.eventData.eventId ,i) = outCY.eventData.interEventResidualNormalized;


    % Clean up
    clear cy_median cy_sigma cy_phi cy_tau

end
rsn_gm_UHS = rsn_vals;
Periods_gmm = Periods;
save(fullfile(data_dir,'CY_forecast_UHS.mat'), ...
     "median_gmm","rotD100_gmm","Periods_gmm","rsn_gm_UHS");
% Find indices of desired periods in Periods array
periods_analysis = [0.1, 1.0, 5.0];
[~, idx_analysis] = ismember(periods_analysis, Periods);
% Check for missing periods
if any(idx_analysis == 0)
    warning('Some periods were not found in Periods array.');
end

% Preallocate index array
idx_closest = zeros(size(exceedMid));

% Loop over each value in exceedMid
for i = 1:length(exceedMid)
    [~, idx_closest(i)] = min(abs(Periods - exceedMid(i)));
end


%% Plot
% === Setup ===
colors = [
    0.00, 0.45, 0.70;  % blue
    0.00, 0.60, 0.50;  % teal
    0.85, 0.33, 0.10   % burnt orange
];
markers = {'o', 's', '^'};  % Circle, square, triangle
jitterAmount = 0.00;

axisMin = 0.0;
axisMax = 3.5;
axisLimits = [axisMin, axisMax];
xticks_vals = linspace(axisMin, axisMax, 6);
yticks_vals = xticks_vals;

% === Create Figure ===
fig_eps_combined_v2 = figure('Name', 'Combined Residual Comparison');

set(fig_eps_combined_v2, 'InvertHardcopy', 'off');
set(fig_eps_combined_v2, 'PaperUnits', 'centimeters');
set(fig_eps_combined_v2, 'PaperSize', [15 14]);
set(fig_eps_combined_v2, 'PaperPositionMode', 'manual');
set(fig_eps_combined_v2, 'PaperPosition', [0 0 15 14]);

hold on; grid on;

% === Plot all 3 groups together ===
for i = 1:length(idx_analysis)
    idx_col = idx_analysis(i);
    x_vals = [];
    if idx_col ~= 0
        for j = 1:length(idx_closest)
            x_vals = [x_vals resid_CY14_RotD100Total(j, idx_closest(j))];
        end
        x_vals = x_vals';
        y_vals = epsilon_dis(:, i);
        idx_val = logical(idx_vals(:, i));
        x_vals = x_vals(idx_val) + jitterAmount * randn(sum(idx_val), 1);
        y_vals = y_vals(idx_val) + jitterAmount * randn(sum(idx_val), 1);

        scatter(x_vals, y_vals, 80, markers{i}, ...
            'MarkerFaceColor', colors(i,:), ...
            'MarkerEdgeColor', 'k', ...
            'LineWidth', 1.0, ...
            'DisplayName', sprintf('$SA(T \\approx %.1f~\\mathrm{s})$', periods_analysis(i)));
    end
end

% Diagonal reference line y = x
plot(axisLimits, axisLimits, 'k--', 'LineWidth', 1.2, 'HandleVisibility', 'off');

% === Axes formatting ===
xlim(axisLimits);
ylim(axisLimits);
xticks(xticks_vals);
yticks(yticks_vals);
axis square;

set(gca, 'FontSize', 13, 'LineWidth', 1.2, 'FontName', 'Times New Roman');

xlabel('Total residual $\varepsilon$', 'Interpreter', 'latex', 'FontSize', 15);
ylabel('$\bar{\varepsilon}$ from disaggregation', 'Interpreter', 'latex', 'FontSize', 15);

% Draw custom black axes crossing at origin
plot([0 0], axisLimits, 'k-', 'LineWidth', 1.2, 'HandleVisibility', 'off');  % Vertical axis
plot(axisLimits, [0 0], 'k-', 'LineWidth', 1.2, 'HandleVisibility', 'off');  % Horizontal axis
lgd = legend('Interpreter', 'latex', 'FontSize', 11, 'Location', 'southeast');

box off;
set(gca, 'Layer', 'top');  % Ensures ticks and markers remain visible
set(gca, 'LooseInset', get(gca, 'TightInset'));  % Reduces padding inside figure

% === Export figure ===
exportgraphics(fig_eps_combined_v2, ...
    fullfile(results_dir,'combined_disagg_vs_real_epsilon.pdf'), ...
    'ContentType','vector','Resolution',1200);