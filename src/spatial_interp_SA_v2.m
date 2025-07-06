function interpolated_SA = spatial_interp_SA_v2(filePath, targetLat, targetLon)
    % INTERPOLATESA Performs spectral acceleration interpolation at a target location
    % 
    % INPUTS:
    %   filePath   - Full path to the CSV file
    %   targetLat  - Target latitude
    %   targetLon  - Target longitude
    %
    % OUTPUT:
    %   interpolated_SA - Interpolated spectral acceleration values for the target location

    % Periods
    T_i = [0.0 0.01 0.02 0.03 0.05 0.075 0.1 0.15 0.2 0.25 0.3 0.4 0.5 0.75 1.0...
        1.5 2.0 3.0 4.0 5.0 7.5 10.0];

    n_data = 9;

    % Load the CSV file
    data = readtable(filePath); % Read the CSV file
    latitudes = data.Latitude;
    longitudes = data.Longitude;

    %% Distance calculation
    

    % 1st alternative
    % Calculate the distance using the Haversine formula
    R = 6371; % Earth's mean radius in kilometers

    % Convert degrees to radians
    lat1 = deg2rad(targetLat);
    lon1 = deg2rad(targetLon);
    lat2 = deg2rad(latitudes);
    lon2 = deg2rad(longitudes);

    % Compute the Haversine distance
    dlat = lat2 - lat1;
    dlon = lon2 - lon1;
    a = sin(dlat / 2).^2 + cos(lat1) .* cos(lat2) .* sin(dlon / 2).^2;
    c = 2 * atan2(sqrt(a), sqrt(1 - a));
    distances = R * c;

    % Find the four closest locations
    [sortedDistances, indices] = sort(distances);
    closestIndices = indices(1:n_data);

    SA_values = data{closestIndices, 3:end}; % Spectral acceleration values
    closestDistances = sortedDistances(1:n_data);
    

    % 2nd alternative
    % R = 6371; % Earth's radius in kilometers
    % x = R * cosd(latitudes) .* cosd(longitudes);
    % y = R * cosd(latitudes) .* sind(longitudes);
    % z = R * sind(latitudes);
    % 
    % targetX = R * cosd(targetLat) * cosd(targetLon);
    % targetY = R * cosd(targetLat) * sind(targetLon);
    % targetZ = R * sind(targetLat);
    % 
    % % Compute Euclidean distances
    % distances = sqrt((x - targetX).^2 + (y - targetY).^2 + (z - targetZ).^2);
    % 
    % % Find the closest locations
    % [sortedDistances, indices] = sort(distances);
    % closestIndices = indices(1:n_data);
    % closestDistances = sortedDistances(1:n_data);
    % SA_values = data{closestIndices, 3:end};

    %3rd alternative
    % a = 6378.137; % Semi-major axis
    % b = 6356.752; % Semi-minor axis
    % f = 1 / 298.257223563; % Flattening
    % 
    % % Convert degrees to radians
    % phi1 = deg2rad(targetLat);
    % phi2 = deg2rad(latitudes);
    % lambda1 = deg2rad(targetLon);
    % lambda2 = deg2rad(longitudes);
    % 
    % % Compute Vincenty distances
    % U1 = atan((1 - f) * tan(phi1));
    % U2 = atan((1 - f) * tan(phi2));
    % L = lambda2 - lambda1;
    % 
    % sinU1 = sin(U1); cosU1 = cos(U1);
    % sinU2 = sin(U2); cosU2 = cos(U2);
    % 
    % lambda = L;
    % for i = 1:100 % Iterative solution
    %     sinLambda = sin(lambda);
    %     cosLambda = cos(lambda);
    %     sinSigma = sqrt((cosU2 .* sinLambda).^2 + ...
    %                     (cosU1 .* sinU2 - sinU1 .* cosU2 .* cosLambda).^2);
    %     cosSigma = sinU1 .* sinU2 + cosU1 .* cosU2 .* cosLambda;
    %     sigma = atan2(sinSigma, cosSigma);
    %     sinAlpha = cosU1 .* cosU2 .* sinLambda ./ sinSigma;
    %     cos2Alpha = 1 - sinAlpha.^2;
    %     cos2SigmaM = cosSigma - 2 * sinU1 .* sinU2 ./ cos2Alpha;
    %     C = f / 16 * cos2Alpha .* (4 + f * (4 - 3 * cos2Alpha));
    %     lambdaPrev = lambda;
    %     lambda = L + (1 - C) .* f .* sinAlpha .* ...
    %              (sigma + C .* sinSigma .* (cos2SigmaM + C .* cosSigma .* (-1 + 2 * cos2SigmaM.^2)));
    %     if max(abs(lambda - lambdaPrev)) < 1e-12
    %         break;
    %     end
    % end
    % 
    % u2 = cos2Alpha .* ((a^2 - b^2) / b^2);
    % A = 1 + u2 / 16384 .* (4096 + u2 .* (-768 + u2 .* (320 - 175 * u2)));
    % B = u2 / 1024 .* (256 + u2 .* (-128 + u2 .* (74 - 47 * u2)));
    % deltaSigma = B .* sinSigma .* (cos2SigmaM + B / 4 .* (cosSigma .* ...
    %     (-1 + 2 * cos2SigmaM.^2) - B / 6 .* cos2SigmaM .* (-3 + 4 * sinSigma.^2) .* (-3 + 4 * cos2SigmaM.^2)));
    % 
    % distances = b .* A .* (sigma - deltaSigma);
    
    % % Find the closest locations
    % [sortedDistances, indices] = sort(distances);
    % closestIndices = indices(1:n_data);
    % closestDistances = sortedDistances(1:n_data);
    % SA_values = data{closestIndices, 3:end};

    % 4th alternative
    %     % Compute Manhattan distances
    % distances = abs(latitudes - targetLat) + abs(longitudes - targetLon);
    % 
    % % Find the closest locations
    % [sortedDistances, indices] = sort(distances);
    % closestIndices = indices(1:n_data);
    % closestDistances = sortedDistances(1:n_data);
    % 
    % SA_values = data{closestIndices, 3:end};
    %% Interpolaton rules
    % weights = 1 ./ closestDistances; %Inverse distance weighting (IDW)

    weights = (1 ./ closestDistances).^2; % Barycentric Interpolation
    % best

    % sigma = mean(closestDistances); % Adjust sigma as needed for sensitivity
    % weights = exp(-(closestDistances.^2) / (2 * sigma^2)); %Gaussian Weighting 2nd best with 9

    % weights = 1 - (closestDistances / max(closestDistances)); % Linear Distance Normalization

    % weights = 1 ./ (1 + closestDistances); % Harmonic Mean Weighting 3rd best 9
    
    % weights = max(0, 1 - (closestDistances / max(closestDistances))); % Triangular Interpolation

    % lambda = mean(closestDistances); % Adjust lambda as needed
    % weights = exp(-closestDistances / lambda); % Exponential Weighting Formula

    % epsilon = 1e-6;
    % weights = -log(closestDistances + epsilon); %logarithmic weights

    weights = weights / sum(weights); % Normalize weights to sum to 1
    interpolated_SA = sum(SA_values .* weights, 1);
end
