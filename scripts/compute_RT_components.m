%% Risk-Targeted Spectrum Calculation for Multi-Period MCE Exceedance
% Computes RT, 84th percentile, and lower limit spectra for each ground motion in the MCER exceedance set.

clear; close all; clc;

% Set up directories
addpath(fullfile('..','src'));              % Add src folder for helper functions
data_dir    = fullfile('..','data');
results_dir = fullfile('..','results','figures');
raw_dir  = fullfile(data_dir, 'raw');         % unprocessed source files
if ~exist(results_dir,'dir'), mkdir(results_dir); end

%% Inputs

% Input Location
% targetLat = 37.4277; % Example target latitude
% targetLon = -122.1701; % Example target longitude

% Periods
T_i = [0.0 0.01 0.02 0.03 0.05 0.075 0.1 0.15 0.2 0.25 0.3 0.4 0.5 0.75 1.0...
        1.5 2.0 3.0 4.0 5.0 7.5 10.0];

file = readtable(fullfile(data_dir,'RSN_Exceed_MultiPeriod_MCE.xlsx'),'NumHeaderLines',1);
rsn_gm = file.Var1;
latitude = file.Var7;
longitude = file.Var8;
vs30 = file.Var10;
% vs30_array = num2cell(vs30);

%% Assign site class and look for origin data

[siteClasses,fileNames_RT,fileName_84th_s,rowIndex_siteClassLL] = assignSiteClass_RT(vs30);

% Risk targeted spectrum
folderPath_RT = fullfile(raw_dir, 'Risk Target MCE Conterminous');

% Lower limit SA
folderPath_LL = fullfile(raw_dir, 'Deterministic-Lower-Limit-SAs_NEHRP-2020.csv');
csvData_LL = readtable(folderPath_LL);



%% Lower Limit Spectrum for each ground motion

SA_LL = zeros(length(siteClasses),length(T_i));

for i = 1:numel(siteClasses)
    row_siteClassLL = rowIndex_siteClassLL(i);
    SA_LL(i,:) = csvData_LL{row_siteClassLL, 2:end};
end


%% Compute Risk-Targeted and 84th-Percentile Spectra (spatial interpolation)

SA_RT = zeros(length(rsn_gm),length(T_i));
SA_84th = zeros(length(rsn_gm),length(T_i));
for i = 1:length(rsn_gm)
    fileName_RT = fileNames_RT{i};
    filePath_RT = fullfile(folderPath_RT, fileName_RT);
    fileName_84th = fileName_84th_s{i};
    filePath_84th = fullfile(folderPath_RT, fileName_84th);
    SA_RT(i,:) = spatial_interp_SA_v2(filePath_RT, latitude(i), longitude(i));
    SA_84th(i,:) = spatial_interp_SA_v2(filePath_84th, latitude(i), longitude(i));
end

%% Save results to data directory

rsn_gm_RT = rsn_gm;
T_RT = T_i;

save(fullfile(data_dir,"rsn_exceed_RiskTarget.mat"), ...
    "rsn_gm_RT", "SA_RT", "SA_84th", "SA_LL", "T_RT");

disp('Risk-targeted, 84th, and lower-limit spectra saved in data/ directory.');











