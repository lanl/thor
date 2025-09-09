function [ val_Hy ] = Hy_value( x, y, z, t )
K = 5;
val_Hy = 0;
for k = 1:K
  val_Hy = val_Hy + k*cos(k*pi*x)*sin(k*pi*y);
end
val_Hy = val_Hy*sin(pi*t);

end