function [ val_Ey ] = tt_Ey_value( xgrid, ygrid, zgrid, t )

mx = numel(xgrid);
my = numel(ygrid);
mz = numel(zgrid);
val_Ey = cell2core(tt_tensor,{zeros(1,mx),zeros(1,my),zeros(1,mz)});

end