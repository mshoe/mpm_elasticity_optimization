function [ shape_val ] = BiLinGrad( x, y, ind_x, ind_y, h )
%BILINGRAD Summary of this function goes here
%   Detailed explanation goes here
dx = x - ind_x*h;
if (-1 <= dx) && (dx < 0)
    shape_valx = 1 + dx/h;
    shape_valdx = 1/h;
elseif (0 <= dx) && (dx < 1)
    shape_valx = 1 - dx/h;
    shape_valdx = -1/h;
else
    shape_valx = 0;
    shape_valdx = 0;
end

dy = y - ind_y*h;
if (-1 <= dy) && (dy < 0)
    shape_valy = 1 + dy/h;
    shape_valdy = 1/h;
elseif (0 <= dy) && (dy < 1)
    shape_valy = 1 - dy/h;
    shape_valdy = -1/h;
else
    shape_valy = 0;
    shape_valdy = 0;
end

shape_val = [shape_valdx * shape_valy, shape_valdy * shape_valx];
end

