%% Set RSNs and Plot Spectra Comparison
% This script plots the response spectra of recordings that exceed target spectra.

clear; clc; % close all
addpath(fullfile('..', 'src'));
% ==== USER PARAMETERS: Edit these to control script behavior ====
rsn_val = 810;              % Example Record Sequence Number (set as needed)
prompt_type = 'MULTP';         % Example: 'RMC', 'RPLL', 'MULTP', 'GMM', 'ALL'
x_max_plot = 2;             % X axis limit
% ===============================================================

% Define relative paths
data_dir = fullfile('..', 'data');
results_dir = fullfile('..', 'results', 'figures');
if ~exist(results_dir, 'dir')
    mkdir(results_dir);
end

% Load GM RSN info
rsn_table = readtable(fullfile(data_dir, 'rsn_set.csv'));
rsn_gm = rsn_table.RecordSequenceNumber;
latitude = rsn_table.StationLatitude;
longitude = rsn_table.StationLongitude;

% Load NGA-West2 meta-data
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

% Load results and supporting data
storedResults_UHS = matfile(fullfile(data_dir, 'UHS_Results.mat'), 'Writable', false);
storedResults_MCE = matfile(fullfile(data_dir, 'MCE_Results.mat'), 'Writable', false);
storedResults_DBE = matfile(fullfile(data_dir, 'DBE_Results.mat'), 'Writable', false);
load(fullfile(data_dir, 'rsn_exceed_RiskTarget.mat'));
load(fullfile(data_dir, 'rsn_exceed_NGA_West_RT.mat'));
load(fullfile(data_dir, 'CY_forecast_UHS.mat'));

% [Optional] Load colormap if needed, otherwise set colors manually
color_1 = [0,0,0];
color_2 = [0.9290 0.6940 0.1250];
color_3 = [0.8500 0.3250 0.0980];
color_4 = [0.3010 0.7450 0.9330];
color_5 = [0 0.4470 0.7410];
color_5p = [0.7, 0.7, 0.7];
color_6 = [0.6350 0.0780 0.1840];

%% Extract EQ information

[EQ_comp_1, EQyear_comp_1, stationname_comp_1, magnitude_comp_1, ...
 distance_comp_1, vs30_comp_1,latitude_comp_1,longitude_comp_1] = ...
    extractEQData(rsn_val, rsn_gm, EQ_name, EQ_year, station_name, magnitude, closest_D, soil_Vs30, latitude, longitude);

%% Set plotting data based on the prompt_type

switch upper(prompt_type)
    case 'RMC'
        % Plots MCER, RT
        T_spect_1 = T_RT;
        index = find(rsn_gm_RT == rsn_val);
        SA_spect_1 = SA_RT(index,:); 
        T_spect_2 = storedResults_MCE.(sprintf('multiT_MCE_T_%d', rsn_val));
        SA_spect_2 = storedResults_MCE.(sprintf('multiT_MCE_SA_%d', rsn_val));
        line_title_1 = 'Probabilistic $\mathrm{MCE}_{\mathrm{R}}$';
        line_title_2 = '$\mathrm{MCE}_{\mathrm{R}}$';
        plot_case = 'RMC';
    case 'RPLL'
        % MCER, RT, Deterministic 84th, Lower Limit
        T_spect_2 = T_RT;
        index = find(rsn_gm_RT == rsn_val);
        SA_spect_2 = SA_RT(index,:); 
        T_spect_1 = storedResults_MCE.(sprintf('multiT_MCE_T_%d', rsn_val));
        SA_spect_1 = storedResults_MCE.(sprintf('multiT_MCE_SA_%d', rsn_val));
        T_spect_3 = T_RT;
        SA_spect_3 = SA_84th(index,:);
        T_spect_4 = T_RT;
        SA_spect_4 = SA_LL(index,:);
        line_title_1 = '$\mathrm{MCE}_{\mathrm{R}}$';
        line_title_2 = 'Probabilistic $\mathrm{MCE}_{\mathrm{R}}$';
        line_title_3 = 'Deterministic $\mathrm{MCE}_{\mathrm{R}}$';
        line_title_4 = 'Deterministic lower limit';
        plot_case = 'RPLL';
    case 'MULTP'
        % Multiperiod ROTD100, DE, MCER, UHS
        T_spect_1 = T_RT';
        index = find(rsn_gm_RT == rsn_val);
        SA_spect_1 = SA_RT(index,:)';
        T_spect_2 = storedResults_DBE.(sprintf('multiT_DE_T_%d', rsn_val));
        SA_spect_2 = storedResults_DBE.(sprintf('multiT_DE_SA_%d', rsn_val));
        line_title_1 = 'Probabilistic $\mathrm{MCE}_{\mathrm{R}}$';
        line_title_2 = '$\mathrm{DE}$';
        T_spect_4 = storedResults_MCE.(sprintf('multiT_MCE_T_%d', rsn_val));
        SA_spect_4 = storedResults_MCE.(sprintf('multiT_MCE_SA_%d', rsn_val));
        line_title_4 = '$\mathrm{MCE}_{\mathrm{R}}$';
        T_spect_5 = storedResults_UHS.(sprintf('T_UHS_%d', rsn_val));
        SA_spect_5 = storedResults_UHS.(sprintf('SA_UHS_%d', rsn_val));
        line_title_5 = 'UHS RP = 2475 yrs';
        plot_case = 'MULTP';
    case 'GMM'
        % ROTD100, DE, RT, MCER, GMM, UHS
        T_spect_1 = T_RT;
        index_uhs = find(rsn_gm_UHS == rsn_val);
        index = find(rsn_gm_RT == rsn_val);
        SA_spect_1 = SA_RT(index,:); 
        T_spect_2 = storedResults_DBE.(sprintf('multiT_DE_T_%d', rsn_val));
        SA_spect_2 = storedResults_DBE.(sprintf('multiT_DE_SA_%d', rsn_val));
        line_title_1 = 'Probabilistic $\mathrm{MCE}_{\mathrm{R}}$';
        line_title_2 = '$\mathrm{DE}$';
        T_spect_3 = Periods_gmm;
        SA_spect_3 = rotD100_gmm(index_uhs,:);
        line_title_3 = '$\mathrm{CY14}$';
        T_spect_4 = storedResults_MCE.(sprintf('multiT_MCE_T_%d', rsn_val));
        SA_spect_4 = storedResults_MCE.(sprintf('multiT_MCE_SA_%d', rsn_val));
        line_title_4 = '$\mathrm{MCE}_{\mathrm{R}}$';
        T_spect_5 = storedResults_UHS.(sprintf('T_UHS_%d', rsn_val));
        SA_spect_5 = storedResults_UHS.(sprintf('SA_UHS_%d', rsn_val));
        line_title_5 = 'UHS RP = 2475 yrs';
        plot_case = 'GMM';
    case 'ALL'
        % ROTD100, TWO AND MULTI DE, TWO AND MULTI MCER, UHS
        T_spect_1 = storedResults_DBE.(sprintf('twoT_DE_T_%d', rsn_val));
        SA_spect_1 = storedResults_DBE.(sprintf('twoT_DE_SA_%d', rsn_val));
        T_spect_2 = storedResults_DBE.(sprintf('multiT_DE_T_%d', rsn_val));
        SA_spect_2 = storedResults_DBE.(sprintf('multiT_DE_SA_%d', rsn_val));
        line_title_1 = 'Two-period $ \mathrm{DE} $';
        line_title_2 = 'Multi-period $ \mathrm{DE} $';
        T_spect_3 = storedResults_MCE.(sprintf('twoT_MCE_T_%d', rsn_val));
        SA_spect_3 = storedResults_MCE.(sprintf('twoT_MCE_SA_%d', rsn_val));
        T_spect_4 = storedResults_MCE.(sprintf('multiT_MCE_T_%d', rsn_val));
        SA_spect_4 = storedResults_MCE.(sprintf('multiT_MCE_SA_%d', rsn_val));
        line_title_3 = 'Two-period $\mathrm{MCE}_{\mathrm{R}}$ \quad';
        line_title_4 = 'Multi-period $\mathrm{MCE}_{\mathrm{R}}$ \quad ';
        T_spect_5 = storedResults_UHS.(sprintf('T_UHS_%d', rsn_val));
        SA_spect_5 = storedResults_UHS.(sprintf('SA_UHS_%d', rsn_val));
        line_title_5 = 'UHS RP = 2475 yrs';
        plot_case = 'ALL';
    otherwise
        error('Specify prompt');
end

% RotD100
Sa_RotD100 = Sa_RotD100_Global(rsn_val,:);

%% Plot Section

if  strcmp(plot_case, 'RMC')
    color_2 = [56 95 150]/255;
    color_3 = [56 95 150]/255;
    fig_1 = figure;
    p_1 = plot(T,Sa_RotD100,'-','Color', color_1,'LineWidth',2);
    hold on
    p_2 = plot(T_spect_1,SA_spect_1,'--','Color', color_2,'LineWidth',2);
    hold on
    p_3 = plot(T_spect_2,SA_spect_2, '-','Color', color_3,'LineWidth',2);
    hold off
    set(gca, 'FontSize', 14, 'LineWidth', 1.5, 'FontName', 'Times New Roman');
    xlabel('Period [s]','FontSize',15,'FontName','Times New Roman')
    ylabel('Spectral Acceleration [g]','FontSize',15,'FontName','Times New Roman');
    legend('boxoff')
    legend_handle = legend([p_1 p_2 p_3],'Seismic record spectrum',line_title_1,line_title_2);
    set(legend_handle,'FontSize',14, 'FontName', 'Times New Roman', 'Interpreter', 'latex');
    axis tight;
    set(gca, 'linewidth', 1.5);
    box off;
    y_max_plot = ceil(max([max(Sa_RotD100), max(SA_spect_1), max(SA_spect_2)]) * 2) / 2;
    axis([0 x_max_plot 0 y_max_plot]);
    xticks(0:1:x_max_plot);
    yticks(0:0.5:y_max_plot);
    xtickformat('%.0f');
    ytickformat('%.1f');
    set(fig_1, 'Color', 'w');
    exportgraphics(fig_1, fullfile(results_dir, ...
        sprintf('Fig_RMC_%d.pdf', rsn_val)), ...
        'ContentType', 'vector', 'Resolution', 1200, ...
        'BackgroundColor', 'white');
    disp(['Exported Fig_RMC_' num2str(rsn_val) '.pdf']);
elseif strcmp(plot_case, 'RPLL')
    color_5 = [0.00 0.45 0.74];
    color_6 = [0.85 0.33 0.10];
    color_3 = [0.93 0.69 0.13];
    color_5p = [0.00 0.60 0.50];
    fig_1 = figure;
    p_2 = plot(T_spect_1,SA_spect_1,'-o','Color', color_5,'LineWidth',2);
    hold on
    p_3 = plot(T_spect_2,SA_spect_2, '-.','Color', color_6,'LineWidth',2);
    hold on
    p_4 = plot(T_spect_3,SA_spect_3,':','Color', color_5p,'LineWidth',2);
    hold on
    p_5 = plot(T_spect_4,SA_spect_4,'--','Color', color_3,'LineWidth',2);
    hold off
    set(gca, 'FontSize', 14, 'LineWidth', 1.5, 'FontName', 'Times New Roman');
    xlabel('Period [s]','FontSize',15,'FontName','Times New Roman')
    ylabel('Spectral Acceleration [g]','FontSize',15,'FontName','Times New Roman');
    legend('boxoff')
    legend_handle = legend([p_2 p_3 p_4 p_5],line_title_1,line_title_2,line_title_3,line_title_4);
    set(legend_handle, 'FontSize', 13,'FontName', 'Times New Roman', 'Interpreter', 'latex');
    axis tight;
    set(gca, 'linewidth', 1.5);
    box off;
    y_max_plot = ceil(max([max(Sa_RotD100), max(SA_spect_1), max(SA_spect_2),...
        max(SA_spect_3),max(SA_spect_4)]) * 2) / 2;
    axis([0 x_max_plot 0 y_max_plot]);
    xticks(0:1:x_max_plot);
    yticks(0:0.5:y_max_plot);
    xtickformat('%.0f');
    ytickformat('%.1f');
    set(fig_1, 'Color', 'w');
    exportgraphics(fig_1, fullfile(results_dir, ...
        sprintf('Fig_RPLL_%d.pdf', rsn_val)), ...
        'ContentType', 'vector', 'Resolution', 1200, ...
        'BackgroundColor', 'white');
    disp(['Exported Fig_RPLL_' num2str(rsn_val) '.pdf']);
elseif strcmp(plot_case, 'MULTP')
    fig_1 = figure;
    p_1 = plot(T,Sa_RotD100,'-','Color', color_1,'LineWidth',2);
    hold on
    p_3 = plot(T_spect_2,SA_spect_2, '-.','Color', color_2,'LineWidth',2);
    hold on
    p_5 = plot(T_spect_4,SA_spect_4,'-.','Color', color_3,'LineWidth',2);
    hold on
    p_6 = plot(T_spect_5,SA_spect_5,'-','Color', color_6,'LineWidth',2);
    hold on
    p_2 = plot(T_spect_1,SA_spect_1,'-.','Color', color_5,'LineWidth',2);
    hold off
    xlabel('Period [s]','FontSize',14,'FontName','Times New Roman')
    ylabel('Spectral Acceleration [g]','FontSize',14,'FontName','Times New Roman');
    legend('boxoff')
    legend_handle = legend([p_1 p_6 p_2 p_5 p_3],'Seismic record spectrum', line_title_5, line_title_1, line_title_4,...
        line_title_2);
    set(legend_handle,'FontSize',13, 'FontName', 'Times New Roman', 'Interpreter', 'latex');
    axis tight; % Tighten the axes to fit the plot
    set(gca, 'linewidth', 1.5);
    set(gca, 'FontSize', 14, 'LineWidth', 1.5, 'FontName', 'Times New Roman');
    box off;
    y_max_plot = ceil(max([max(Sa_RotD100), max(SA_spect_1), max(SA_spect_2),...
        max(SA_spect_4),max(SA_spect_5)]) * 2) / 2;
    axis([0 x_max_plot 0 y_max_plot]);
    xticks(0:1:x_max_plot);
    yticks(0:0.5:y_max_plot);
    xtickformat('%.0f');
    ytickformat('%.1f');
    set(fig_1, 'Color', 'w');
    exportgraphics(fig_1, fullfile(results_dir, ...
        sprintf('Fig_MULTP_%d.pdf', rsn_val)), ...
        'ContentType', 'vector', 'Resolution', 1200, ...
        'BackgroundColor', 'white');
    disp(['Exported Fig_MULTP_' num2str(rsn_val) '.pdf']);
elseif strcmp(plot_case, 'GMM')
    fig_1 = figure;
    p_1 = plot(T,Sa_RotD100,'-','Color', color_1,'LineWidth',2);
    hold on
    p_3 = plot(T_spect_2,SA_spect_2, '-.','Color', color_2,'LineWidth',2);
    hold on
    p_4 = plot(T_spect_3,SA_spect_3,'--','Color', color_4,'LineWidth',2);
    hold on
    p_5 = plot(T_spect_4,SA_spect_4,'-','Color', color_3,'LineWidth',2);
    hold on
    p_6 = plot(T_spect_5,SA_spect_5,'-','Color', color_6,'LineWidth',2);
    hold on
    p_2 = plot(T_spect_1,SA_spect_1,'-.','Color', color_5,'LineWidth',2);
    hold off
    xlabel('Period [s]','FontSize',15,'FontName','Times New Roman')
    ylabel('Spectral Acceleration [g]','FontSize',15,'FontName','Times New Roman');
    legend('boxoff')
    legend_handle = legend([p_1 p_6 p_2 p_5 p_3 p_4],'Seismic record spectrum', line_title_5, line_title_1, line_title_4,...
        line_title_2, line_title_3);
    set(legend_handle,'FontSize',13, 'FontName', 'Times New Roman', 'Interpreter', 'latex');
    set(legend_handle, 'Location', 'southwest');
    str_cell_array = { ...
    "RSN: " + num2str(rsn_val), ...
    "EQ: " + EQ_comp_1, ...
    "M = " + num2str(magnitude_comp_1), ...
    "Rrup = " + num2str(distance_comp_1) + " km" ...
    };
    
    str_title = strjoin(string(str_cell_array), "; ");
    
    str_cell_array_2 = { ...
        "Vs30 = " + num2str(vs30_comp_1) + " m/s", ...
        "Station: " + stationname_comp_1 ...
    };
    
    str_title_2 = strjoin(string(str_cell_array_2), "; ");
    
    title({ ...
        str_title, ...
        str_title_2 ...
    }, 'FontName', 'Times New Roman', 'FontWeight', 'normal')
    axis tight;
    set(gca, 'linewidth', 1.5);
    set(gca, 'FontSize', 14, 'LineWidth', 1.5, 'FontName', 'Times New Roman');
    set(gca, 'XScale', 'log', 'YScale', 'log');
    box off;
    grid on;
    set(gca, 'MinorGridLineStyle', ':', 'MinorGridAlpha', 0.05);
    y_max_plot = 10;
    axis([1e-2 x_max_plot 1e-3 y_max_plot]);
    xtickformat('%.0f');
    ytickformat('%.1f');
    axis square;
    set(fig_1, 'Color', 'w');
    exportgraphics(fig_1, fullfile(results_dir, ...
        sprintf('Fig_GMM_%d.pdf', rsn_val)), ...
        'ContentType', 'vector', 'Resolution', 1200, ...
        'BackgroundColor', 'white');
    disp(['Exported Fig_GMM_' num2str(rsn_val) '.pdf']);
else
    fig_1 = figure;
    p_1 = plot(T,Sa_RotD100,'-','Color', color_1,'LineWidth',2);
    hold on
    p_2 = plot(T_spect_1,SA_spect_1,'--','Color', color_2,'LineWidth',2);
    hold on
    p_3 = plot(T_spect_2,SA_spect_2, '-.','Color', color_3,'LineWidth',2);
    hold on
    p_4 = plot(T_spect_3,SA_spect_3,'--','Color', color_4,'LineWidth',2);
    hold on
    p_5 = plot(T_spect_4,SA_spect_4,'-.','Color', color_5,'LineWidth',2);
    hold on
    p_6 = plot(T_spect_5,SA_spect_5,'-','Color', color_6,'LineWidth',2);
    hold off
    xlabel('Period [s]','FontSize',15,'FontName','Times New Roman')
    ylabel('Spectral Acceleration [g]','FontSize',15,'FontName','Times New Roman');
    legend('boxoff')
    legend_handle = legend([p_1 p_2 p_3 p_4 p_5 p_6],'Seismic record spectrum',line_title_1,line_title_2,...
        line_title_3,line_title_4,line_title_5);
    set(legend_handle, 'FontSize', 14,'FontName', 'Times New Roman', 'Interpreter', 'latex');
    axis tight;
    set(gca, 'linewidth', 1.5);
    box off;
    set(gca, 'FontSize', 14, 'LineWidth', 1.5, 'FontName', 'Times New Roman');
    y_max_plot = ceil(max([max(Sa_RotD100), max(SA_spect_1), max(SA_spect_2),...
        max(SA_spect_3),max(SA_spect_4),max(SA_spect_5)]) * 2) / 2;
    axis([0 x_max_plot 0 y_max_plot]);
    xticks(0:1:x_max_plot);
    yticks(0:0.5:y_max_plot);
    xtickformat('%.0f');
    ytickformat('%.1f');
    set(fig_1, 'Color', 'w');
        exportgraphics(fig_1, fullfile(results_dir, ...
        sprintf('Fig_ALL_%d.pdf', rsn_val)), ...
        'ContentType', 'vector', 'Resolution', 1200, ...
        'BackgroundColor', 'white');
    disp(['Exported Fig_ALL_' num2str(rsn_val) '.pdf']);
end


disp(" ")