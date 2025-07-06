function siteclass = assignSiteClass(data)
    % Initialize siteclass cell array
    siteclass = cell(length(data),1);

    % Loop through the data to classify site
    for i = 1:length(data)
        vs30 = data(i);

        if vs30 > 152.4 && vs30 <= 213.36
            sc_i = 'DE';
        elseif vs30 > 213.36 && vs30 <= 304.80
            sc_i = 'D';
        elseif vs30 > 304.80 && vs30 <= 441.96
            sc_i = 'CD';
        elseif vs30 > 441.96 && vs30 <= 640.08
            sc_i = 'C';
        elseif vs30 > 640.08 && vs30 <= 914.40
            sc_i = 'BC';
        elseif vs30 > 914.40 && vs30 <= 1524.00
            sc_i = 'B';
        elseif vs30 > 1524.00
            sc_i = 'A';
        else
            sc_i = 'E'; % For values outside the specified brackets
        end

        siteclass{i} = sc_i;
    end
end
