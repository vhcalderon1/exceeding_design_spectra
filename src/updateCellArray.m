function T = updateCellArray(T, x, y)
    % Check if T is empty and initialize if necessary
    if isempty(T)
        T = {x, y}; % Initialize with the first entry
        return;
    end

    % Extract all the x values from T
    xValues = T(:, 1);
    
    % Check if x is already in T
    idx = find(cell2mat(xValues) == x, 1);
    
    if isempty(idx)
        % x is not in T, add a new entry
        T = [T; {x, y}];
    else
        % x exists in T, append y to the corresponding vector
        T{idx, 2} = [T{idx, 2}; y];
    end
end