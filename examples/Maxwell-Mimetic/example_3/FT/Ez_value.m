function [ val_Ez ] = Ez_value( x, y, z, t )
K =5;
val_Ez = 0;
for k = 1:K
  val_Ez = val_Ez + sin(k*pi*x)*sin(k*pi*y);
end
  val_Ez = val_Ez*cos(pi*t);
end