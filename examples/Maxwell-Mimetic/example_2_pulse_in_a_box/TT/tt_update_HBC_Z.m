function [ tEx, tEy ] = tt_update_HBC_Z( tEx, tEy, dx, dy, dz, nx, ny, nz, time, tol )

% --------------------
% HBC at "z0" (bottom)
% --------------------
%% update Ex at kz =1
tEx = tt_1D_set_zeros(tEx,3,[1,nz]);

%% update Ey at kz = 1
tEy = tt_1D_set_zeros(tEy,3,[1,nz]);

end
