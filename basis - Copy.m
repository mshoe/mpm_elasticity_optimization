function [ normalizer ] = basis( ind, h, x )
%1D piecewise linear basis function

    % index of grid: ind
    % grid-spacing: h
    % position: x

    % Piecewise linear basis function for x
    if (ind - 1)*h <= x && x <= ind * h
        normalizer = 1 + (x - ind * h) / h;
    elseif ind * h <= x && x <= (ind + 1)*h
        normalizer = 1 - (x - ind * h) / h;
    else
        normalizer = 0;
    end
    

end

