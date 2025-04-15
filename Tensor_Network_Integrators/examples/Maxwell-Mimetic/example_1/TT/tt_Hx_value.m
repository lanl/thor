function [ val_Hx ] = tt_Hx_value( xgrid, ygrid, zgrid, t)

G = {-sin(pi.*xgrid), cos(pi.*ygrid), ones(1,numel(zgrid)).*sin(pi.*t)};
val_Hx = cell2core(tt_tensor,G);

end