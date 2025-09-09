function [ val_Hy ] = tt_Hy_value( xgrid, ygrid, zgrid, t, tol)
K=5;
val_Hy = 0;
for k = 1:K
G = {k*cos(k*pi.*xgrid), sin(k*pi.*ygrid), ones(1,numel(zgrid))};
val_Hy = round(val_Hy + cell2core(tt_tensor, G), tol);
end
val_Hy = round(val_Hy*sin(pi*t), tol);
end