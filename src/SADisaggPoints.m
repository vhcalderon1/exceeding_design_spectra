function allValues = SADisaggPoints(gm_latitude, gm_longitude, gm_vs30, RT, saPeriods)
%SADISAGGPOINTS  Query USGS disagg at multiple SA periods and extract data with retries.
%
% INPUTS:
%   gm_latitude  - Site latitude
%   gm_longitude - Site longitude
%   gm_vs30      - Vs30 value at site
%   RT           - Return period (e.g., 2475)
%   saPeriods    - Array of spectral periods (e.g., [0.01, 0.2, 5.0])
%
% OUTPUT:
%   allValues    - 1 x (3 * numel(saPeriods)) numeric array of "value" fields
%
% EXAMPLE:
%   vals = SADisaggPoints(37.726, -122.424, 760, 2475, [0.01, 0.2, 5.0])

    % Cap Vs30 within the valid range [150, 1500]
    gm_vs30 = round(gm_vs30);
    gm_vs30 = max(min(gm_vs30, 1500), 150);

    % Base URL for the USGS disaggregation service
    baseUrl = 'https://earthquake.usgs.gov/ws/nshmp/conus-2023/dynamic/disagg/';

    % Number of data points per call (assumed 3 based on the output format)
    nPerCall = 3;
    nPeriods = numel(saPeriods);
    allValues = zeros(1, nPerCall * nPeriods);

    % Retry parameters
    max_retries = 10000;         % Maximum retry attempts
    initial_timeout = 50;     % Initial timeout in seconds
    timeout_increment = 50;   % Increment timeout by this amount per retry

    % Loop over the specified SA periods
    for i = 1:nPeriods
        % Convert SA period to "SAxPx" format
        saStr = strrep(sprintf('%.1f', saPeriods(i)), '.', 'P');  % Guarantees "1.0" → "1P0"
        imtStr = ['SA' saStr];

        % Construct the API call URL
        fullUrl = sprintf('%s%.6f/%.6f/%d/%d?out=DISAGG_DATA,GMM,SOURCE&imt=%s', ...
            baseUrl, gm_longitude, gm_latitude, round(gm_vs30), RT, imtStr);

        % Retry loop
        num_retries = 0;
        dataRetrieved = false;
        while num_retries < max_retries
            try
                % Set the timeout dynamically based on retries
                options = weboptions('Timeout', initial_timeout + timeout_increment * num_retries);
                
                % Query the web service (JSON response)
                dataStruct = webread(fullUrl, options);

                % Extract the required data
                extractedData = dataStruct.response.disaggs.data(1).summary(4).data;
                currentValues = [extractedData.value];  % 1×3 numeric

                % Store values in the correct slice
                idxStart = (i - 1) * nPerCall + 1;
                idxEnd   = i * nPerCall;
                allValues(idxStart:idxEnd) = currentValues;

                % If successful, break out of retry loop
                dataRetrieved = true;
                break;
            catch exception
                % If an error occurs, retry
                disp(['Retrying... Attempt ' num2str(num_retries + 1) ' of ' num2str(max_retries)]);
                num_retries = num_retries + 1;
            end
        end

        % If we fail after max retries, throw an error
        if ~dataRetrieved
            error('Failed to retrieve data from USGS after %d attempts for SA = %s.', max_retries, imtStr);
        end
    end
end
