function [ Jx, Jy, Jz ] = tt_setup_Jc( nx, ny, nz )


Jx = cell2core(tt_tensor,{zeros(1,nx-1),zeros(1,ny), zeros(1,nz)});
Jy = cell2core(tt_tensor,{zeros(1,nx),zeros(1,ny-1), zeros(1,nz)});
Jz = cell2core(tt_tensor,{zeros(1,nx),zeros(1,ny), zeros(1,nz-1)});
% setup
% Jx = zeros(nx-1,ny,nz);
% Jy = zeros(nx,ny-1,nz);
% Jz = zeros(nx,ny,nz-1);

end