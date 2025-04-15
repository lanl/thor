function [ val_Ez ] = tt_Ez_value( xgrid, ygrid, zgrid, t )

mz = numel(zgrid);
val_Ez = cos(pi*t)*cell2core(tt_tensor,{sin(pi.*xgrid),...
      sin(pi.*ygrid), ones(1,mz)});
end