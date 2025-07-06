function [two_T_DE_T, two_T_DE_ord, multi_T_DE_T, multi_T_DE_ord, two_T_MCE_T, two_T_MCE_ord, multi_T_MCE_T, multi_T_MCE_ord] = ASCE722_spectrums(gm_latitude, gm_longitude, riskCategory, site_class, gm_title)
    % Inputs
    % gm_latitude: latitude of the ground motion site
    % gm_longitude: longitude of the ground motion site
    % riskCategory: Risk category of the site
    % site_class: Site class of the site
    % gm_title: Title for the ground motion data

    % Outputs
    % two_T_DE_T: Two-period design spectrum periods
    % two_T_DE_ord: Two-period design spectrum ordinates
    % multi_T_DE_T: Multi-period design spectrum periods
    % multi_T_DE_ord: Multi-period design spectrum ordinates
    % two_T_MCE_T: Two-period MCE spectrum periods
    % two_T_MCE_ord: Two-period MCE spectrum ordinates
    % multi_T_MCE_T: Multi-period MCE spectrum periods
    % multi_T_MCE_ord: Multi-period MCE spectrum ordinates

    baseUrl = 'https://earthquake.usgs.gov/ws/designmaps/asce7-22.json?';
    params = append('latitude=', string(gm_latitude), '&longitude=', string(gm_longitude), '&riskCategory=', riskCategory, '&siteClass=', char(site_class), '&title=', gm_title);
    params = convertStringsToChars(params);

    % Combine the base URL with the query parameters
    fullUrl = [baseUrl params];

    % Retry logic
    max_retries = 50; % Maximum number of retries
    initial_timeout = 50; % Initial timeout in seconds
    timeout_increment = 50; % Increment timeout by this amount on each retry

    num_retries = 0;
    data = [];

    while num_retries < max_retries
        try
            % Use webread to make the request and retrieve data
            options = weboptions('Timeout', initial_timeout + timeout_increment * num_retries);
            data = webread(fullUrl, options);
            if isempty(data)
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

    if isempty(data)
        error('Failed to retrieve data from USGS after %d attempts.', max_retries);
    end

    % Saving relevant parameters
    % Two-Period Spectrum
    sd1 = data.response.data.sd1;
    sds = data.response.data.sds;
    t0 = data.response.data.t0;
    tl = data.response.data.tl;
    ts = data.response.data.ts;

    T = [0:0.025:1 1.05:0.05:10];
    two_T_DE_T = T;
    two_T_MCE_T = T;
    two_T_DE_ord = calculateDesignSpectrum(sd1, sds, t0, tl, ts, two_T_DE_T);
    two_T_DE_ord = two_T_DE_ord';

    SA = calculateDesignSpectrum(sd1, sds, t0, tl, ts, two_T_MCE_T);
    two_T_MCE_ord = SA * 1.5;
    two_T_MCE_ord = two_T_MCE_ord';

    % Multi-period
    multi_T_DE_T = data.response.data.multiPeriodDesignSpectrum.periods;
    multi_T_DE_ord = data.response.data.multiPeriodDesignSpectrum.ordinates;

    multi_T_MCE_T = data.response.data.multiPeriodMCErSpectrum.periods;
    multi_T_MCE_ord = data.response.data.multiPeriodMCErSpectrum.ordinates;
end