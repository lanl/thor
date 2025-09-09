function [ tHx, tHy, tHz ] = tt_update_HBC_H(tHx, tHy, tHz, dx, dy, dz, nx, ny, nz, time, tol)

%% update tHx at ix = [1,nx] -- zero boundary
tHx = tt_1D_set_zeros(tHx,1,[1,nx]);
%% update tHy at jy = [1, ny]
tHy = tt_1D_set_zeros(tHy,2,[1,ny]);
%% update tHz at kz = [1, ny]
tHz = tt_1D_set_zeros(tHz,3,[1,nz]);

end