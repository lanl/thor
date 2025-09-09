function [ tEz, tEx ] = tt_update_HBC_Y( tEz, tEx, dx, dy, dz, nx, ny, nz, time, tol)

% -----------
% HBC at "y0"
% -----------
%% update tEx at jy = 1,ny
tEx = tt_1D_set_zeros(tEx,2,[1,ny]);

%% update tEz at jy=1,ny
tEz = tt_1D_set_zeros(tEz,2,[1,ny]);

end