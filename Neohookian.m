function [ energy ] = Neohookian( Def, MoE, PoR, N )
%NEOHOOKIAN Summary of this function goes here
%   Detailed explanation goes here
% Def is the deformation gradient list
% MoE is modulus of elasticity
% PoR is Poisson's ratio
% n is # of particles

% Calculate Lame coefficients
mew = MoE / (2 * (1 + PoR));
lambda = MoE * PoR / ((1 + PoR) * (1 - 2 * PoR));
% From wiki article on Neo-Hookian
C1 = mew/2;
D1 = lambda/2;

% Calculate invariants
I1 = zeros(N, 1);
I3 = zeros(N, 1);
J = zeros(N, 1);
I1_bar = zeros(N, 1);
for p = 1:N
    I1(p) = trace(Def(:, :, p)' * Def(:, :, p));
    I3(p) = det(Def(:, :, p))^2;
   J(p) = det(Def(:, :, p));
   I1_bar(p) = J(p)^(-1) * I1(p);
end

% Calculate neohookian energy
energy = zeros(N, 1);
for p = 1:N
%    energy(p) = C1 * (I1_bar(p) - 2) + D1 * (J(p)  - 1)^2;
    energy(p) = mew/2 * (I1(p) - log(I3(p)) - 2) ...
        + lambda/8 * (log(I3(p)))^2;
end
end

