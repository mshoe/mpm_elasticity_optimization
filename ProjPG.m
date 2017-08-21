function [ Fg ] = ProjPG( Fp, Wpg, N, m )
%PROJPG Summary of this function goes here
%   Detailed explanation goes here

    Fg = zeros(m^2, 1);
    funFg = @(Fg) 0.5 * sum((sum(repmat(Fg, [1 N]) .* Wpg.', 1) - Fp.').^2, 2);
    
    Fg = fminunc(funFg, Fg);

end

