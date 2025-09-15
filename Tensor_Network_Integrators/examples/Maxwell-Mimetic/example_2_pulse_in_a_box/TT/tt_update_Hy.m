function [ tHy ] = tt_update_Hy( tHy, tEz, tEx, dz, dx, dt, tol )

%compute dEzdx
dEzdx = tt_FD(tEz,1);
dEzdx = round(dEzdx,tol);
%compute dExdz
dExdz = tt_FD(tEx,3);
dExdz = round(dExdz,tol);

Dterm = round((dt/dx)*dEzdx - (dt/dz)*dExdz, tol);
tHy = round(tHy + Dterm,tol);
% tHy = round(tHy, tol);


% for ix=1:nx-1
%   for jy=1:ny
%     for kz=1:nz-1
%       Hy(ix, jy, kz) = Hy(ix, jy, kz) + dt/(dz*dx) * (...
%         + dz*(Ez(ix+1, jy, kz) - Ez(ix, jy, kz))...
%         - dx*(Ex(ix, jy, kz+1) - Ex(ix, jy, kz)) ) ;
%     end
%   end
% end

end
