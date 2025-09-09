function [ tEy ] = tt_update_Ey( tEy, tJy, tHz, tHx, dz, dx, dt, nx, ny, nz, tol )

tEy = round(tEy + dt*tJy,tol);

%Enlarge and finite diffrance the tHx and tHz
dHxdz = (dt/dz)*tt_FD(tHx,3);
dHxdz = tt_zeropad_fn(dHxdz,3,[1,1]); 

dHzdx = (dt/dx)*tt_FD(tHz,1);
dHzdx = tt_zeropad_fn(dHzdx,1,[1,1]);
tupd = round(dHxdz - dHzdx,tol);


%Cut and pad the tupd terms
% adding = [2,0,2];
% dims = {(2:nx-1),(1:ny-1),(2:nz-1)};
% tupd = cut_and_pad(tupd,dims,adding);

%Updating tEy
tEy = round(tEy + tupd, tol);

end