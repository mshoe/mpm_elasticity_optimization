function [ Defg ] = DeformationGradGrid( Up, WpgGrad, N, m )
%DEFORMATIONGRADGRID Summary of this function goes here
%   Detailed explanation goes here

    UGradg = zeros(2, 2, m^2);
    
    Defg = zeros(2, 2, m^2);
  	
    
    for p = 1:N
        for ind = 1:m^2
            UGradg(:, :, ind) = UGradg(:, :, ind) + WpgGrad(p, :, ind).' * Up(p, :);
        end
    end
    
    for ind = 1:m^2
        Defg(:, :, ind) = inv(eye(2) - UGradg(:, :, ind));
    end
end

