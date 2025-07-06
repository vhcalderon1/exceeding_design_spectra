% Inputs

% x_o : The point at which f(x) is to be interpolated

% x : An array of x values where f(x) is known

% fx : An array of corresponding f(x) values

% Outputs

% fx_o : Result of the interpolation

function fx_o = linear_interpol(x_o, x, fx)

% Check if the input is within the bounds of the x values
    if x_o < x(1) || x_o > x(end)
        error('x_o is out of the bounds of the x values.');
    end

    % Find the bracket interval
    for i = 1:length(x)-1
        if x_o >= x(i) && x_o <= x(i+1)
            % Linear interpolation formula
            fx_o = fx(i) + (fx(i+1) - fx(i)) * (x_o - x(i)) / (x(i+1) - x(i));
            return;
        end
    end

    % If x_o exactly matches the last x element, handle it separately
    if x_o == x(end)
        fx_o = fx(end);
    else
        error('No bracket found. This should not happen.');
    end

end