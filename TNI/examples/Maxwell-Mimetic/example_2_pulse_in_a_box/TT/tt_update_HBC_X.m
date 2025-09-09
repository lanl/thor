function [ tEy, tEz ] = tt_update_HBC_X( tEy, tEz, dx, dy, dz, nx, ny, nz, time, tol)

%% update tEy at ix = 1,nx
tEy = tt_1D_set_zeros(tEy,1,[1,nx]);

%% update tEz at ix = 1,nx
tEz = tt_1D_set_zeros(tEz,1,[1,nx]);


end