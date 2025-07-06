function [T_vals, SA_vals] = UHS_Alaska2023_v1(gm_latitude, gm_longitude, gm_vs30, RT)
    % Inputs
    % gm_latitude: latitude of the ground motion site
    % gm_longitude: longitude of the ground motion site
    % gm_vs30: Vs30 value of the ground motion site
    % RT: Return period

    % Outputs
    % T_vals: Period values
    % SA_vals: Spectral acceleration values

    AFE_0 = 1 / RT;
    gm_vs30 = round(gm_vs30);

    % Cap of Vs30
    if gm_vs30 > 1500
        gm_vs30 = 1500;
    elseif gm_vs30 < 150
        gm_vs30 = 150;
    end

    %% Call

    max_retries = 50; % Maximum number of retries
    initial_timeout = 50; % Initial timeout in seconds
    timeout_increment = 50; % Increment timeout by this amount on each retry

    % Control of fails
    T_i = [];
    num_retries = 0;

    while num_retries < max_retries
        try
            baseUrl = 'https://earthquake.usgs.gov/ws/nshmp/alaska-2023/dynamic/hazard/';
            params = append(string(gm_longitude), '/', string(gm_latitude), '/', string(gm_vs30));
            params = convertStringsToChars(params);
            
            % Combine the base URL with the query parameters
            fullUrl = [baseUrl params];
            
            % Use webread to make the request and retrieve data
            options = weboptions('Timeout', initial_timeout + timeout_increment * num_retries);
            data = webread(fullUrl, options);
            
            n_entry_try = 3;
            T_i = data.response.hazardCurves(n_entry_try).imt.display;
            if isempty(T_i)
                num_retries = num_retries + 1;
            else
                break;
            end

        catch exception
            % In case of a connection timeout or any other error, retry
            disp(['Retrying... Attempt ' num2str(num_retries + 1) ' of ' num2str(max_retries)]);
            num_retries = num_retries + 1;
        end
    end

    if isempty(T_i)
        error('Failed to retrieve data from USGS after %d attempts.', max_retries);
    end

    %% Interpolation

    n_data = length(data.response.hazardCurves) - 2;
    T_vals = zeros(n_data, 1);
    SA_vals = zeros(n_data, 1);

    % For one value
    for i = 1:n_data
        n_entry = i + 2; % Corresponds to Period - starts from 3
        T_i = data.response.hazardCurves(n_entry).imt.display;
        T_i = regexp(T_i, '[\d.]+', 'match');
        T_i = str2double(T_i{1});
        
        SA_T = data.response.hazardCurves(n_entry).data(1).values.xs; % In g's
        AFE_T = data.response.hazardCurves(n_entry).data(1).values.ys; % Annual freq of exceedance
        
        % Remove trailing zeros
        non_zero_idx = find(AFE_T ~= 0);
        SA_T = SA_T(non_zero_idx);
        AFE_T = AFE_T(non_zero_idx);

        % Semi-log interpolation
        SA_0 = interp1(log(AFE_T), SA_T, log(AFE_0), 'linear');

        T_vals(i) = T_i;
        SA_vals(i) = SA_0;
    end
end
