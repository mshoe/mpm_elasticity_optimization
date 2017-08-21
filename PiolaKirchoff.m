function [ PiolaKirchoff ] = PiolaKirchoff( Def, MoE, PoR, N )
%PIOLAKIRCHOFF Summary of this function goes here
%   Detailed explanation goes here
% Calculate Lame coefficients
mew = MoE / (2 * (1 + PoR));
lambda = MoE * PoR / ((1 + PoR) * (1 - 2 * PoR));
J = zeros(N, 1);
PiolaKirchoff = zeros(2, 2, N); % d Energy / d DefGrad
for p = 1:N
    J(p) = det(Def(:, :, p));
    PiolaKirchoff(:, :, p) = mew * Def(:, :, p) + (lambda * log(J(p)) - mew) * inv(Def(:, :, p).'); 
end

end

