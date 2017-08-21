function [ L ] = Lin2Dto1D( F )
%LIN2DTO1D Summary of this function goes here
%   Detailed explanation goes here
    length = size(F, 1) * size(F, 2);
    L = zeros(length, 1);
    A = 1:2:length;
    B = 2:2:length;
    L(A) = F(:,1);
    L(B) = F(:,2);
end

