%% MCER Exceedance Results Table Creation
% Compiles summary tables for ground motions exceeding code spectra (MCER).

clear; close all; clc;

% Set up directories
addpath(fullfile('..','src'));              % Add src folder for helper functions
data_dir    = fullfile('..','data');
results_dir = fullfile('..','results','figures');
if ~exist(results_dir,'dir'), mkdir(results_dir); end

%% Inputs

% Load metadata for all RSNs
file = readtable(fullfile(data_dir, 'rsn_set.csv'),'NumHeaderLines',1);
rsn_gm = file.Var1;
latitude = file.Var2;
longitude = file.Var3;
magnitude_unit = file.Var5;

% Load main NGA-West2 database
matFileName = fullfile(data_dir, 'NGA_W2_corr_meta_data.mat');
variablesToLoad = {'Sa_RotD100', 'Periods', 'soil_Vs30','station_name','EQ_year','EQ_name','closest_D','magnitude'}; % Replace with your variable names
loadRotD = load(matFileName, variablesToLoad{:});
Sa_RotD100 = loadRotD.Sa_RotD100;
Sa_RotD100_Global = Sa_RotD100;
Periods = loadRotD.Periods;
soil_Vs30 = loadRotD.soil_Vs30;
station_name = loadRotD.station_name;
EQ_year = loadRotD.EQ_year;
EQ_name = loadRotD.EQ_name;
closest_D = loadRotD.closest_D;
magnitude = loadRotD.magnitude;

% Load exceeding records
% 1st component: two-period MCE
% 2nd component: multi-period MCE
load(fullfile(data_dir,'rsn_exceed_NGA_West_MCE.mat'))
load(fullfile(data_dir,'T_vs_RSN_MCE.mat'))
load(fullfile(data_dir,'intervals_exceed_MCE.mat'))

%% Table creation

% Two period MCE
% Data selection

% Location
% First component
latitude_comp_1 = zeros(length(rsn_great_1),1);
longitude_comp_1 = zeros(length(rsn_great_1),1);
magnitude_unit_comp_1 = cell(length(rsn_great_1), 1);
for i = 1: length(rsn_great_1)
    index_csv = find(rsn_gm == rsn_great_1(i));
    latitude_comp_1(i) = latitude(index_csv);
    longitude_comp_1(i) = longitude(index_csv);
    magnitude_unit_comp_1{i} = magnitude_unit{index_csv};
end

% EQ data
EQ_comp_1 = cell(length(rsn_great_1), 1);
EQyear_comp_1 = zeros(length(rsn_great_1),1);
stationname_comp_1 = cell(length(rsn_great_1), 1);
magnitude_comp_1 = zeros(length(rsn_great_1),1);
distance_comp_1 = zeros(length(rsn_great_1),1);
vs30_comp_1 = zeros(length(rsn_great_1),1);
for i = 1: length(rsn_great_1)
    EQ_comp_1{i} = EQ_name{rsn_great_1(i)};
    EQyear_comp_1(i) = EQ_year{rsn_great_1(i)};
    stationname_comp_1{i} = station_name{rsn_great_1(i)};
    magnitude_comp_1(i) = magnitude(rsn_great_1(i));
    distance_comp_1(i) = closest_D(rsn_great_1(i));
    vs30_comp_1(i) = soil_Vs30(rsn_great_1(i));
end

% Second component
latitude_comp_2 = zeros(length(rsn_great_2),1);
longitude_comp_2 = zeros(length(rsn_great_2),1);
magnitude_unit_comp_2 = cell(length(rsn_great_2), 1);
for i = 1: length(rsn_great_2)
    index_csv = find(rsn_gm == rsn_great_2(i));
    latitude_comp_2(i) = latitude(index_csv);
    longitude_comp_2(i) = longitude(index_csv);
    magnitude_unit_comp_2{i} = magnitude_unit{index_csv};
end

% EQ data
EQ_comp_2 = cell(length(rsn_great_2), 1);
EQyear_comp_2 = zeros(length(rsn_great_2),1);
stationname_comp_2 = cell(length(rsn_great_2), 1);
magnitude_comp_2 = zeros(length(rsn_great_2),1);
distance_comp_2 = zeros(length(rsn_great_2),1);
vs30_comp_2 = zeros(length(rsn_great_2),1);
for i = 1: length(rsn_great_2)
    EQ_comp_2{i} = EQ_name{rsn_great_2(i)};
    EQyear_comp_2(i) = EQ_year{rsn_great_2(i)};
    stationname_comp_2{i} = station_name{rsn_great_2(i)};
    magnitude_comp_2(i) = magnitude(rsn_great_2(i));
    distance_comp_2(i) = closest_D(rsn_great_2(i));
    vs30_comp_2(i) = soil_Vs30(rsn_great_2(i));
end


%% Create and save summary tables (in data folder)

title_table = {'RSN', 'Earthquake Name', 'Year','Station Name', ...
    'Earthquake Magnitude','Magnitude Type','Station Latitude','Station Longitude', ...
    'ClstD (km)','Vs30 (m/s)','T_{min}','T_{max}'};

table_twoperiod = table(rsn_great_1, EQ_comp_1, EQyear_comp_1, stationname_comp_1, ...
    magnitude_comp_1, magnitude_unit_comp_1, latitude_comp_1, longitude_comp_1, ...
    distance_comp_1, vs30_comp_1, t_int_1(:,2), t_int_1(:,3), ...
    'VariableNames', title_table);
filename1 = fullfile(data_dir, 'RSN_Exceed_TwoPeriod_MCE.xlsx');
writetable(table_twoperiod, filename1);

table_multip = table(rsn_great_2, EQ_comp_2, EQyear_comp_2, stationname_comp_2, ...
    magnitude_comp_2, magnitude_unit_comp_2, latitude_comp_2, longitude_comp_2, ...
    distance_comp_2, vs30_comp_2, t_int_2(:,2), t_int_2(:,3), ...
    'VariableNames', title_table);
filename2 = fullfile(data_dir, 'RSN_Exceed_MultiPeriod_MCE.xlsx');
writetable(table_multip, filename2);

disp(" ")
