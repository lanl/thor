function [ val_Hx ] = Hx_value( x, y, z, t)
    val_Hx = -sin(pi*x)*cos(pi*y)*sin(pi*t);
end