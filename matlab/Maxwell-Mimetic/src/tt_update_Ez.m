function [ tEz ] = tt_update_Ez( tEz, tJz, tHx, tHy, dx, dy, dt, nx, ny, nz, tol)

%% update tEz
tEz = round(tEz + dt*tJz, tol);

%Enlarge and finite diffrance the tHy and tHx
dHydx = (dt/dx)*tt_FD(tHy,1);
dHydx = tt_zeropad_fn(dHydx,1,[1,1]);
dHxdy = (dt/dy)*tt_FD(tHx,2);
dHxdy = tt_zeropad_fn(dHxdy,2,[1,1]);
tupd = round(dHydx - dHxdy,tol);

%Cut and pad the tupd terms
% adding = [2,2,0];
% dims = {(2:nx-1),(2:ny-1),(1:nz-1)};
% tupd = cut_and_pad(tupd,dims,adding);

%Updating tEz
tEz = round(tEz + tupd, tol);


end

