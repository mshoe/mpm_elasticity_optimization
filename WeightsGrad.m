function [ WpgGrad ] = WeightsGrad( Xp, GridCells, N, m, h )
%WEIGHTSGRAD Summary of this function goes here
%   Detailed explanation goes here

    WpgGrad = zeros(N, 2, m^2);
    for p = 1:N
        for ind = 1:m^2
            WpgGrad(p, :, ind) = BiLinGrad(Xp(p,1), Xp(p,2), GridCells(ind, 1), ...
                GridCells(ind, 2), h);
        end 
    end

end

