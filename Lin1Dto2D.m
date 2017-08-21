function [ F ] = Lin1Dto2D( L )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

    length = size(L, 1);
    F = zeros(length/2, 2);
    A = 1:2:length;
    B = 2:2:length;
    F(:,1) = L(A);
    F(:,2) = L(B);

end

