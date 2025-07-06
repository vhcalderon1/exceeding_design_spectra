function [EQ_comp_2, EQyear_comp_2, stationname_comp_2, magnitude_comp_2, distance_comp_2, vs30_comp_2,latitude_comp_2,longitude_comp_2] = extractEQData(rsn_great_2, rsn_gm, EQ_name, EQ_year, station_name, magnitude, closest_D, soil_Vs30, latitude, longitude)
    % Pre-allocate output variables
    EQ_comp_2 = cell(length(rsn_great_2), 1);
    EQyear_comp_2 = zeros(length(rsn_great_2), 1);
    stationname_comp_2 = cell(length(rsn_great_2), 1);
    magnitude_comp_2 = zeros(length(rsn_great_2), 1);
    distance_comp_2 = zeros(length(rsn_great_2), 1);
    vs30_comp_2 = zeros(length(rsn_great_2), 1);
    
    % Extract information from input vectors
    for i = 1:length(rsn_great_2)
        EQ_comp_2{i} = EQ_name{rsn_great_2(i)};
        EQyear_comp_2(i) = EQ_year{rsn_great_2(i)};
        stationname_comp_2{i} = station_name{rsn_great_2(i)};
        magnitude_comp_2(i) = magnitude(rsn_great_2(i));
        distance_comp_2(i) = closest_D(rsn_great_2(i));
        vs30_comp_2(i) = soil_Vs30(rsn_great_2(i));
    end

    latitude_comp_2 = zeros(length(rsn_great_2),1);
    longitude_comp_2 = zeros(length(rsn_great_2),1);
    % magnitude_unit_comp_2 = cell(length(rsn_great_2), 1);

    for i = 1: length(rsn_great_2)
        index_csv = find(rsn_gm == rsn_great_2(i));
        if isempty(index_csv)
            continue
        end
        latitude_comp_2(i) = latitude(index_csv);
        longitude_comp_2(i) = longitude(index_csv);
        % magnitude_unit_comp_2{i} = magnitude_unit{index_csv};
    end

end
