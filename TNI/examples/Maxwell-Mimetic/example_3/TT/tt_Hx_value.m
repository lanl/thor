function [ val_Hx ] = tt_Hx_value( xgrid, ygrid, zgrid, t, tol)
K=5;
val_Hx = 0;
for k = 1:K
G = {-k*sin(k*pi.*xgrid), cos(k*pi.*ygrid), ones(1,numel(zgrid))};
val_Hx = round(val_Hx + cell2core(tt_tensor,G), tol);
end
val_Hx = round(val_Hx*sin(pi*t),tol);
end