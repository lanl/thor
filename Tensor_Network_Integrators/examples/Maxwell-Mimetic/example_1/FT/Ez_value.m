function [ val_Ez ] = Ez_value( x, y, z, t )

val_Ez = sin(pi*x)*sin(pi*y)*cos(pi*t);
end