function [ tHx ] = tt_update_Hx(tHx, tEy, tEz, dy, dz, dt, tol)

%compute dEy/dz
dEydz = tt_FD(tEy,3);
dEydz = round(dEydz,tol);
%compute dEz/dy
dEzdy = tt_FD(tEz,2);
dEzdy = round(dEzdy,tol);
Dterm = round((dt/dz)*dEydz - (dt/dy)*dEzdy, tol);
tHx = round(tHx + Dterm,tol);
% tHx = round(tHx, tol);

% for ix=1:nx
%   for jy=1:ny-1
%     for kz=1:nz-1
%       Hx(ix, jy, kz) = Hx(ix, jy, kz) + dt/(dy*dz) * (...
%         + dy*(Ey(ix, jy, kz+1) - Ey(ix, jy, kz))...
%         - dz*(Ez(ix, jy+1, kz) - Ez(ix, jy, kz)) ) ;
% 
%     end
%   end
% end

end