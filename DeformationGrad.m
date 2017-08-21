function [ Defp ] = DeformationGrad( Ug, WpgGrad, N, m )
%DEFORMATION Summary of this function goes here
%   Detailed explanation goes here

    UGradp = zeros(2, 2, N);
    
    Defp = zeros(2, 2, N);
  	
    
    for p = 1:N
        for ind = 1:m^2
            UGradp(:, :, p) = UGradp(:, :, p) + WpgGrad(p, :, ind).' * Ug(ind, :);
        end
    end
    
    % Compute deformation gradients
    for p = 1:N
        Defp(:, :, p) = inv(eye(2) - UGradp(:, :, p));
    end
    
end

