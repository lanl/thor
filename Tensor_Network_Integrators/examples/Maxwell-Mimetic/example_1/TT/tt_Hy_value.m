function [ val_Hy ] = tt_Hy_value( xgrid, ygrid, zgrid, t )

G = {cos(pi.*xgrid), sin(pi.*ygrid), ones(1,numel(zgrid)).*sin(pi.*t)};
val_Hy = cell2core(tt_tensor, G);

end