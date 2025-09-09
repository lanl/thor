function [ Ex, Ey, Ez ] = tt_setup_E(dx, dy, dz, nx, ny, nz, X0 )
%% pulse in the box problem -- initial condition
% Ez = exp(-beta*(x^2 + y^2)) = exp(-beta*x^2)*exp(-beta*y^2)
beta = 25;
% setup Ex,
% xgrid = X0 + dx*([1:nx-1] -1/2);
% ygrid = X0 + dy*([1:ny] - 1);
% zgrid = X0 + dz*([1:nz]- 1);
% Ex = tt_Ex_value(xgrid, ygrid, zgrid, time);
% Ex = eps2*tt_ones2({nx-1,ny,nz});
Ex = tt_zeros2({nx-1,ny,nz});
% Ex = cell2core(tt_tensor,{exp(-beta.*(xgrid.^2)), exp(-beta.*(ygrid.^2)),...
%   exp(-beta.*(zgrid.^2))});

%setup Ey
% xgrid = dx*([1:nx] -1);
% ygrid = dy*([1:ny-1] - 1/2);
% zgrid = dz*([1:nz]- 1);
% Ey = tt_Ey_value(xgrid, ygrid, zgrid, time);
% Ey = eps2*tt_ones2({nx,ny-1,nz});
Ey = tt_zeros2({nx,ny-1,nz});
%setup Ez
xgrid = X0 + dx*([1:nx] -1);
ygrid = X0 + dy*([1:ny] - 1);
zgrid = X0 + dz*([1:nz-1]- 1/2);
% Ez = tt_zeros2({nx,ny,nz-1});
% Ez = cell2core(tt_tensor,{exp(-beta.*(xgrid.^2)), exp(-beta.*(ygrid.^2)),...
%   exp(-beta.*(zgrid.^2))});
Ez = cell2core(tt_tensor,{exp(-beta.*(xgrid.^2)), exp(-beta.*(ygrid.^2)),...
  ones(1,numel(zgrid))});



end