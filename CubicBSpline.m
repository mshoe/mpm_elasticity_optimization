function [ shape_val ] = CubicBSpline( x, y, ind_x, ind_y, h )
% Shape function

dx = x - ind_x*h;
if (-2*h <= dx) && (dx < -h)
    shape_valx = 1/(6*h^3)*dx^3 + 1/(h^2)*dx^2 + 2/h*dx + 4/3;
elseif (-h <= dx) && (dx < 0)
    shape_valx = -1/(2*h^3)*dx^3 - 1/(h^2)*dx^2 + 4/3; 
elseif (0 <= dx) && (dx < h)
    shape_valx = 1/(2*h^3)*dx^3 - 1/(h^2)*dx^2 + 4/3;
elseif (h <= dx) && (dx < 2*h);
    shape_valx = -1/(6*h^3)*dx^3 + 1/(h^2)*dx^2 - 2/h*dx + 4/3;
else
    shape_valx = 0;
end

dy = y - ind_y*h;
if (-2*h <= dy) && (dy < -h)
    shape_valy = 1/(6*h^3)*dy^3 + 1/(h^2)*dy^2 + 2/h*dy + 4/3;
elseif (-h <= dy) && (dy < 0)
    shape_valy = -1/(2*h^3)*dy^3 - 1/(h^2)*dy^2 + 4/3; 
elseif (0 <= dy) && (dy < h)
    shape_valy = 1/(2*h^3)*dy^3 - 1/(h^2)*dy^2 + 4/3;
elseif (h <= dy) && (dy < 2*h);
    shape_valy = -1/(6*h^3)*dy^3 + 1/(h^2)*dy^2 - 2/h*dy + 4/3;
else
    shape_valy = 0;
end

shape_val = shape_valx * shape_valy;
end

