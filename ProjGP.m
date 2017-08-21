function [ Fp ] = ProjGP( Fg , Wpg, N, m)
%PROJGP Summary of this function goes here
%   Detailed explanation goes here
    Fp = zeros(N, 2);
    
    for p = 1:N
        for ind = 1:m^2
            Fp(p, :) = Fp(p, :) +  Wpg(p, ind) * Fg(ind, :);
        end
    end

end

