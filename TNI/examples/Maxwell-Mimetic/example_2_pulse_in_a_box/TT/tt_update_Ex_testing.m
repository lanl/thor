function [ tEx ] = tt_update_Ex_testing( tEx, tHy, tHz, dy, dz, dt, nx, ny, nz, tol)

%% update tEx
% tEx = round(tEx + dt*tJx, tol);
%Enlarge and finite diffrance the tHx and tHz
dHzdy = (dt/dy)*tt_FD(tHz,2);
dHzdy = tt_zeropad_fn(dHzdy,2,[1,1]);

dHydz = (dt/dz)*tt_FD(tHy,3);
dHydz = tt_zeropad_fn(dHydz,3,[1,1]);

tupd = round((dHzdy - dHydz),tol);

%Cut and pad the tupd terms
% adding = [0,2,2];
% dims = {(1:nx-1),(2:ny-1),(2:nz-1)};
% tupd = cut_and_pad(tupd,dims,adding);

%Update
tEx = round(tEx + tupd,tol);

end