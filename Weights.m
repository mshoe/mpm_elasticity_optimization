function [ Wpg ] = Weights( Xp, GridCells, N, m, h )
%WEIGHTS Summary of this function goes here
%   Detailed explanation goes here

    Wpg = zeros(N, m^2);
    for p = 1:N
        for ind = 1:m^2
            Wpg(p, ind) = BilinearInterpolation(Xp(p,1), Xp(p,2), GridCells(ind, 1), ...
                GridCells(ind, 2), h);
        end 
    end
    
end

