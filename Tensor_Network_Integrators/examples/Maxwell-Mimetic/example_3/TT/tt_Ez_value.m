function [ val_Ez ] = tt_Ez_value( xgrid, ygrid, zgrid, t, tol)
K=5;
val_Ez = 0;
for k = 1:K
  mz = numel(zgrid);
  val_Ez = round(val_Ez + cell2core(tt_tensor,{sin(k*pi.*xgrid),...
    sin(k*pi.*ygrid), ones(1,mz)}), tol);
end
val_Ez = cos(pi*t)*val_Ez;
end