function [ycom_1, xValues_1, numElements_1, limit_y_1] = exceed_Periods(T_uni_great_Two_DBE, Periods)
    % Extract x values and the number of elements in each cell
    xValues_1 = cell2mat(T_uni_great_Two_DBE(:, 1));
    numElements_1 = cellfun(@(v) numel(v), T_uni_great_Two_DBE(:, 2));

    % Sort x values and reorder numElements accordingly
    [xValues_1, sortIdx] = sort(xValues_1);
    numElements_1 = numElements_1(sortIdx);

    % Initialize ycom_1 with zeros, same size as Periods
    ycom_1 = zeros(size(Periods,2),1);

    % Find indices in Periods that match xValues_1
    [~, idx] = ismember(xValues_1, Periods);

    % Assign numElements_1 to ycom_1 at the matching indices
    ycom_1(idx(idx > 0)) = numElements_1(idx > 0);

    ycom_1 = ycom_1';
    
    % Calculate limit_y_1
    limit_y_1 = ceil(max(numElements_1) / 5) * 5;
end
