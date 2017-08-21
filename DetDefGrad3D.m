function [ DetDefp ] = DetDefGrad3D( Defp, N)
%DETDEFGRAD3D Summary of this function goes here
%   Detailed explanation goes here

    DetDefp = zeros(N, 1);
    for p = 1:N
        DetDefp(p) = det(Defp(:,:,p));
    end
end

