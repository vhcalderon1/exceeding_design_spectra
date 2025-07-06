function [siteclass,fileNames_RT,fileName_84th_s,row_LL] = assignSiteClass(data)
    % Initialize siteclass cell array
    siteclass = cell(length(data),1);
    fileNames_RT = cell(length(data),1);
    fileName_84th_s = cell(length(data),1);
    row_LL = zeros(length(data),1);

    % Loop through the data to classify site
    for i = 1:length(data)
        vs30 = data(i);

        if vs30 > 152.4 && vs30 <= 213.36
            sc_i = 'DE';
            fileName = 'ConUS-2018_MaxDirection-RTSAs_vs30=185-siteClass=DE_NEHRP-2020.csv';
            fileName_84th = 'ConUS-2018_MaxDirection-84thSAs_vs30=185-siteClass=DE_NEHRP-2020.csv';
            rowIndex_siteClassLL = 7;
        elseif vs30 > 213.36 && vs30 <= 304.80
            sc_i = 'D';
            fileName = 'ConUS-2018_MaxDirection-RTSAs_vs30=260-siteClass=D_NEHRP-2020.csv';
            fileName_84th = 'ConUS-2018_MaxDirection-84thSAs_vs30=260-siteClass=D_NEHRP-2020.csv';
            rowIndex_siteClassLL = 6;
        elseif vs30 > 304.80 && vs30 <= 441.96
            sc_i = 'CD';
            fileName = 'ConUS-2018_MaxDirection-RTSAs_vs30=365-siteClass=CD_NEHRP-2020.csv';
            fileName_84th = 'ConUS-2018_MaxDirection-84thSAs_vs30=365-siteClass=CD_NEHRP-2020.csv';
            rowIndex_siteClassLL = 5;
        elseif vs30 > 441.96 && vs30 <= 640.08
            sc_i = 'C';
            fileName = 'ConUS-2018_MaxDirection-RTSAs_vs30=530-siteClass=C_NEHRP-2020.csv';
            fileName_84th = 'ConUS-2018_MaxDirection-84thSAs_vs30=530-siteClass=C_NEHRP-2020.csv';
            rowIndex_siteClassLL = 4;
        elseif vs30 > 640.08 && vs30 <= 914.40
            sc_i = 'BC';
            fileName = 'ConUS-2018_MaxDirection-RTSAs_vs30=760-siteClass=BC_NEHRP-2020.csv';
            fileName_84th = 'ConUS-2018_MaxDirection-84thSAs_vs30=760-siteClass=BC_NEHRP-2020.csv';
            rowIndex_siteClassLL = 3;
        elseif vs30 > 914.40 && vs30 <= 1524.00
            sc_i = 'B';
            fileName = 'ConUS-2018_MaxDirection-RTSAs_vs30=1080-siteClass=B_NEHRP-2020.csv';
            fileName_84th = 'ConUS-2018_MaxDirection-84thSAs_vs30=1080-siteClass=B_NEHRP-2020.csv';
            rowIndex_siteClassLL = 2;
        elseif vs30 > 1524.00
            sc_i = 'A';
            fileName = 'ConUS-2018_MaxDirection-RTSAs_vs30=1500-siteClass=A_NEHRP-2020.csv';
            fileName_84th = 'ConUS-2018_MaxDirection-84thSAs_vs30=1500-siteClass=A_NEHRP-2020.csv';
            rowIndex_siteClassLL = 1;
        else
            sc_i = 'E'; % For values outside the specified brackets
            fileName = 'ConUS-2018_MaxDirection-RTSAs_vs30=150-siteClass=E_NEHRP-2020.csv';
            fileName_84th = 'ConUS-2018_MaxDirection-84thSAs_vs30=150-siteClass=E_NEHRP-2020.csv';
            rowIndex_siteClassLL = 8;
        end

        siteclass{i} = sc_i;
        fileNames_RT{i} = fileName;
        fileName_84th_s{i} = fileName_84th;
        row_LL(i) = rowIndex_siteClassLL;
    end
end
