function [ val_Hx ] = Hx_value( x, y, z, t)
K = 5;
val_Hx = 0;
for k = 1:K
  val_Hx = val_Hx -k*sin(k*pi*x)*cos(k*pi*y);
end
val_Hx = val_Hx*sin(pi*t);
end