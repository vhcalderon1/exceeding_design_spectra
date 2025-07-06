%% NGA-West2 RotD100 Exceedance Analysis (DE Level)
% Identifies and records ground motions that exceed the ASCE 7-22 DE spectra.

clear; close all; clc;

% Set up directories
addpath(fullfile('..','src'));              % Add src folder for helper functions
data_dir    = fullfile('..','data');
results_dir = fullfile('..','results','figures');
if ~exist(results_dir,'dir'), mkdir(results_dir); end

%% Inputs

tic;

% Load GM RSN
file = readtable(fullfile(data_dir, 'rsn_set.csv'),'NumHeaderLines',1);
rsn_gm = file.Var1;
latitude = file.Var2;
longitude = file.Var3;

riskCategory = 'III';
gm_title = 'Example';

% Load NGA-West2 spectral & site data
matFileName = fullfile(data_dir, 'NGA_W2_corr_meta_data.mat');
variablesToLoad = {'Sa_RotD100', 'Periods', 'soil_Vs30','station_name','EQ_year','EQ_name','closest_D'}; % Replace with your variable names
loadRotD = load(matFileName, variablesToLoad{:});
Sa_RotD100 = loadRotD.Sa_RotD100;
Sa_RotD100_Global = Sa_RotD100;
Periods = loadRotD.Periods;
soil_Vs30 = loadRotD.soil_Vs30;
station_name = loadRotD.station_name;
EQ_year = loadRotD.EQ_year;
EQ_name = loadRotD.EQ_name;
closest_D = loadRotD.closest_D;
T = Periods;

% Set up storage for DBE results in data folder
resultFile = fullfile(data_dir, 'DBE_Results.mat');
if ~isfile(resultFile)
    save(resultFile, 'rsn_gm'); % Create file with one variable
end
storedResults = matfile(resultFile, 'Writable', true);

%% Loop through each RSN and check for code spectrum exceedances

not_found_RSN = [];
rsn_great_2 = [0];
rsn_great_1 = [0];
ratio_dif_SA_2 = [];
mult_dif_SA_2 = [];
T_great_2 = [];
ratio_dif_SA_1 = [];
mult_dif_SA_1 = [];
T_great_1 = [];
T_uni_great_2 = {};
T_uni_great_1 = {};
t_int_1 = [];
t_int_2 = [];

for j = 1:length(rsn_gm)
    disp(j)
    Tmax_2 = 0;
    Tmax_1 = 0;
    Tmin_1 = 10;
    Tmin_2 = 10;

    rsn_val = rsn_gm(j);

    % Check if results already exist
    if isprop(storedResults, sprintf('twoT_DE_T_%d', rsn_val))
        two_T_DE_T = storedResults.(sprintf('twoT_DE_T_%d', rsn_val));
        two_T_DE_ord = storedResults.(sprintf('twoT_DE_SA_%d', rsn_val));
        multi_T_DE_T = storedResults.(sprintf('multiT_DE_T_%d', rsn_val));
        multi_T_DE_ord = storedResults.(sprintf('multiT_DE_SA_%d', rsn_val));
    else
        % ASCE 7-22 Spectrum
        gm_latitude = latitude(j);
        gm_longitude = longitude(j);
        gm_site_class = assignSiteClass(soil_Vs30(rsn_val));

        % ASCE7 Tool threshold
        [two_T_DE_T,two_T_DE_ord,multi_T_DE_T,multi_T_DE_ord,two_T_MCE_T,two_T_MCE_ord,multi_T_MCE_T,multi_T_MCE_ord] = ASCE722_spectrums_v2(gm_latitude,gm_longitude,riskCategory,gm_site_class,gm_title);

        % Store the results
        storedResults.(sprintf('twoT_DE_T_%d', rsn_val)) = two_T_DE_T;
        storedResults.(sprintf('twoT_DE_SA_%d', rsn_val)) = two_T_DE_ord;
        storedResults.(sprintf('multiT_DE_T_%d', rsn_val)) = multi_T_DE_T;
        storedResults.(sprintf('multiT_DE_SA_%d', rsn_val)) = multi_T_DE_ord;
    end

    % Condition to analyze existing values
    Sa_RotD100 = Sa_RotD100_Global(rsn_val,:);
    if any(Sa_RotD100 > 0)

        % Comparing objective

        T_comp_1 = two_T_DE_T;
        Sa_comp_1 = two_T_DE_ord;
        T_comp_2 = multi_T_DE_T;
        Sa_comp_2 = multi_T_DE_ord;

        % Condition to evaluate

        if max(Sa_RotD100) < min(Sa_comp_1) && max(Sa_RotD100) < min(Sa_comp_2)
            continue
        end

        % Period of analysis
    
        if max(T) <= max(T_comp_1) || max(T) <= max(T_comp_2)
            T_a1 = T;
        else
            if length(T_comp_1) < length(T_comp_2)
                T_a1 = T_comp_1;
                T_a2 = T_comp_2;
            else
                T_a1 = T_comp_2;
                T_a2 = T_comp_1;
            end
        end
    
        index_T = find(T == max(T_a1));
    
        for i = 1:index_T
    
            index_T1 = find(T_a1 == T(i));
            index_T2 = find(T_a2 == T(i));
    
            % Checking if GM spectrum is greater than Design spectrum
    
            if isempty(index_T1) || isempty(index_T2)
    
                % Design sprectrum
                SA_2 = linear_interpol(T(i), T_comp_2, Sa_comp_2);
                SA_1 = linear_interpol(T(i), T_comp_1, Sa_comp_1);

                % Difference condition

                if SA_2 < Sa_RotD100(i)

                    ratio_dif_SA_2 = [ratio_dif_SA_2; (Sa_RotD100(i)-SA_2)/SA_2];
                    mult_dif_SA_2 = [mult_dif_SA_2; Sa_RotD100(i)/SA_2];
                    T_great_2 = [T_great_2; T(i)];

                    T_uni_great_2 = updateCellArray(T_uni_great_2, T(i), rsn_val);

                    if T(i) < Tmin_2
                        Tmin_2 = T(i);
                        Tmax_2 = Tmin_2;
                    end
                    if T(i) > Tmax_2 
                        Tmax_2 = T(i);
                    end

                end

                if SA_1 < Sa_RotD100(i)

                    ratio_dif_SA_1 = [ratio_dif_SA_1; (Sa_RotD100(i)-SA_1)/SA_1];
                    mult_dif_SA_1 = [mult_dif_SA_1; Sa_RotD100(i)/SA_1];
                    T_great_1 = [T_great_1; T(i)];

                    T_uni_great_1 = updateCellArray(T_uni_great_1, T(i), rsn_val);

                    if T(i) < Tmin_1
                        Tmin_1 = T(i);
                        Tmax_1 = Tmin_1;
                    end
                    if T(i) > Tmax_1 
                        Tmax_1 = T(i);
                    end

                end

                % Count of unique contribution
                % Selecting RSN that are greater than code spectrums
    
                if SA_2 < Sa_RotD100(i) && rsn_great_2(end) ~= rsn_val
    
                    rsn_great_2 = [rsn_great_2; rsn_val];

                end

                if SA_1 < Sa_RotD100(i) && rsn_great_1(end) ~= rsn_val

                    rsn_great_1 = [rsn_great_1; rsn_val];
    
                end
    
            else
    
                % Design sprectrum
                if length(Sa_comp_2) == length(T_a2)
                    SA_2 = Sa_comp_2(index_T2);
                    SA_1 = Sa_comp_1(index_T1);
                else
                    SA_2 = Sa_comp_2(index_T1);
                    SA_1 = Sa_comp_1(index_T2);
                end

                % Difference condition

                if SA_2 < Sa_RotD100(i)

                    ratio_dif_SA_2 = [ratio_dif_SA_2; (Sa_RotD100(i)-SA_2)/SA_2];
                    mult_dif_SA_2 = [mult_dif_SA_2; Sa_RotD100(i)/SA_2];
                    T_great_2 = [T_great_2; T(i)];
                    T_uni_great_2 = updateCellArray(T_uni_great_2, T(i), rsn_val);

                    if T(i) < Tmin_2
                        Tmin_2 = T(i);
                        Tmax_2 = Tmin_2;
                    end
                    if T(i) > Tmax_2 
                        Tmax_2 = T(i);
                    end

                end

                if SA_1 < Sa_RotD100(i)

                    ratio_dif_SA_1 = [ratio_dif_SA_1; (Sa_RotD100(i)-SA_1)/SA_1];
                    mult_dif_SA_1 = [mult_dif_SA_1; Sa_RotD100(i)/SA_1];
                    T_great_1 = [T_great_1; T(i)];
                    T_uni_great_1 = updateCellArray(T_uni_great_1, T(i), rsn_val);

                    if T(i) < Tmin_1
                        Tmin_1 = T(i);
                        Tmax_1 = Tmin_1;
                    end
                    if T(i) > Tmax_1 
                        Tmax_1 = T(i);
                    end

                end

                % Selecting RSN that are greater than code spectrums
    
                if SA_2 < Sa_RotD100(i) && rsn_great_2(end) ~= rsn_val
    
                    rsn_great_2 = [rsn_great_2; rsn_val];

                end
    
                if SA_1 < Sa_RotD100(i) && rsn_great_1(end) ~= rsn_val
    
                    rsn_great_1 = [rsn_great_1; rsn_val];
    
                end
    
            end
    
        end

    if Tmin_2<10 && Tmax_2>0
        t_int_2 = [t_int_2; rsn_val Tmin_2 Tmax_2];
    end
    if Tmin_1<10 && Tmax_1>0
        t_int_1 = [t_int_1; rsn_val Tmin_1 Tmax_1];
    end

    else
    
        not_found_RSN = [not_found_RSN; rsn_val];
        
    end

end
toc;

%% Saving data (all outputs go to data folder)
rsn_great_1 = rsn_great_1(2:end);   % Remove initial placeholder zero
rsn_great_2 = rsn_great_2(2:end);

save(fullfile(data_dir,"rsn_exceed_NGA_West_DBE.mat"), ...
    "rsn_great_1","not_found_RSN","rsn_great_2");
save(fullfile(data_dir,"val_exceed_NGA_West_DBE.mat"), ...
    "T_great_2", "ratio_dif_SA_2","mult_dif_SA_2", ...
    "T_great_1", "ratio_dif_SA_1","mult_dif_SA_1");
save(fullfile(data_dir,"T_vs_RSN_DBE.mat"), ...
    "T_uni_great_2","T_uni_great_1");
save(fullfile(data_dir,"intervals_exceed_DBE.mat"), ...
    "t_int_2","t_int_1");

disp('Finished DE exceedance analysis.');