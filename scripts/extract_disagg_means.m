%% Extract Disaggregation Data for Multi-Period MCE Exceedance
% Extract disaggregation from USGS Hazard Tools for each ground motion and merges results with metadata.

clear; clc; % close all

% Set up directories
addpath(fullfile('..','src'));              % Add src folder for helper functions
data_dir    = fullfile('..','data');
results_dir = fullfile('..','results','figures');
if ~exist(results_dir,'dir'), mkdir(results_dir); end

%% Inputs

% Define the full path to the .mat file
rsn_selected_Path = fullfile(data_dir, 'RSN_Exceed_MultiPeriod_MCE.xlsx'); % Modify by RSN_Exceed_UHS.xlsx to get disagg of UHS RSNs
rsn_exceed_multi_MCE = readtable(rsn_selected_Path,'NumHeaderLines',1);

rsn_gm = rsn_exceed_multi_MCE.Var1;
latitude = rsn_exceed_multi_MCE.Var7;
longitude = rsn_exceed_multi_MCE.Var8;
vs30 = rsn_exceed_multi_MCE.Var10;

RT        = 2475;                         % Return period
saPeriods = [0.1, 1.0, 5.0];              % Periods of analysis

%% Dissagregation output

dissagregation_meanvals = zeros(length(rsn_gm),length(saPeriods)*3);
for i = 1:length(rsn_gm)
    allValues = SADisaggPoints(latitude(i), longitude(i), vs30(i), RT, saPeriods);
    dissagregation_meanvals(i,:) = allValues;
end

%% Save matrix for reproducibility

original_table = readtable(rsn_selected_Path);
dissagregation_rsn = rsn_gm;
save(fullfile(data_dir, "dissagregation_exceed_MCE.mat"), ...
     "dissagregation_rsn", "dissagregation_meanvals");

%% Create New Column Names for Dissagregation Data

numPeriods = numel(saPeriods);  % Number of periods (expected to be 3)

% Preallocate a cell array for the new column names (9 columns in total)
newColNames = cell(1, numPeriods * 3); 

for i = 1:numPeriods
    % Convert the period to a string with one decimal place, then replace '.' with 'P'
    numStr = sprintf('%.1f', saPeriods(i));
    numStr = strrep(numStr, '.', 'P');
    
    % Create names for the i-th group of 3 columns using the modified number string
    newColNames{3*(i-1)+1} = sprintf('mean_m_T_%s', numStr);
    newColNames{3*(i-1)+2} = sprintf('mean_r_T_%s', numStr);
    newColNames{3*(i-1)+3} = sprintf('mean_epsiolon_T_%s', numStr);
end


%% Merge with original table

dissag_table = array2table(dissagregation_meanvals, 'VariableNames', newColNames);
combined_table = [original_table, dissag_table];

%% Save combined table (CSV in data/)

writetable(combined_table, fullfile(data_dir, "RSN_Exceed_MultiPeriod_MCE_with_Dissag.csv"));
% For UHS save as: RSN_Exceed_UHS_with_Dissag

disp('Disaggregation values computed and combined table saved to data/.');
disp(head(combined_table));
