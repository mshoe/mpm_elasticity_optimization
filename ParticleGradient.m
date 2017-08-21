function [FGradp ] = ParticleGradient( Fg, GridCelli, GridCellj, ...
    BasisGrad1, BasisGrad2, BasisGrad3, BasisGrad4, m, N)
%PARTICLEGRADIENT Summary of this function goes here
%   Detailed explanation goes here
    FGradp = zeros(2, 2, N);
    Fg1 = Fg((GridCellj)*m + GridCelli+1, :);
    Fg2 = Fg((GridCellj+1)*m + GridCelli+1, :);
    Fg3 = Fg((GridCellj+1)*m + GridCelli+1+1, :);
    Fg4 = Fg((GridCellj)*m + GridCelli+1+1, :);
    for p = 1:N
        %numPtcls = 1;%GridCellParticleCount(gridCell(p));
        FGradp(:, :, p) = FGradp(:, :, p) + BasisGrad1(p,:).' * Fg1(p,:);%/numPtcls;
        FGradp(:, :, p) = FGradp(:, :, p) + BasisGrad2(p,:).' * Fg2(p,:);%/numPtcls;
        FGradp(:, :, p) = FGradp(:, :, p) + BasisGrad3(p,:).' * Fg3(p,:);%/numPtcls;
        FGradp(:, :, p) = FGradp(:, :, p) + BasisGrad4(p,:).' * Fg4(p,:);%/numPtcls;
    end
end

