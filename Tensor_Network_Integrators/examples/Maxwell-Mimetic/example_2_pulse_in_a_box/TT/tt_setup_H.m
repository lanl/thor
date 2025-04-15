function [ Hx, Hy, Hz ] = tt_setup_H(dx, dy, dz, nx, ny, nz, X0)

% setup Hx,Hy,Hz in tt format

% setup Hx,
% xgrid = dx*([1:nx] -1);
% ygrid = dy*([1:ny-1] - 1/2);
% zgrid = dz*([1:nz-1]- 1/2);
% Hx = tt_Hx_value(xgrid, ygrid, zgrid, time);
% Hx = eps2*tt_ones2({nx,ny-1,nz-1});
Hx = tt_zeros2({nx,ny-1,nz-1});

% setup Hy
% xgrid = dx*([1:nx-1] -1/2);
% ygrid = dy*([1:ny] - 1);
% zgrid = dz*([1:nz-1]- 1/2);
% Hy = tt_Hy_value2(xgrid, ygrid, zgrid, time);
% Hy = eps2*tt_ones2({nx-1,ny,nz-1});
Hy = tt_zeros2({nx-1,ny,nz-1});
% setup Hz

% xgrid = dx*([1:nx-1] -1/2);
% ygrid = dy*([1:ny-1] - 1/2);
% zgrid = dz*([1:nz]- 1);
% Hz = tt_Hz_value(xgrid, ygrid, zgrid, time);
% Hz = eps2*tt_ones2({nx-1,ny-1,nz});
Hz = tt_zeros2({nx-1,ny-1,nz});
end 