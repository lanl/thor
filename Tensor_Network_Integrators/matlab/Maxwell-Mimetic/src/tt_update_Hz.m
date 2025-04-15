function [ tHz ] = tt_update_Hz( tHz, tEx, tEy, dx, dy, dt, tol)

% compute dExdy
dExdy = tt_FD(tEx,2);
% dExdy = round(dExdy,tol);
%compute dEydx
dEydx = tt_FD(tEy,1);
% dEydx = round(dEydx,tol);
% tHz = tHz + dt/(dx*dy)*round(dx*dExdy - dy*dEydx,tol);
tHz = tHz + round((dt/dy)*dExdy - (dt/dx)*dEydx,tol);
tHz= round(tHz, tol);

% 
% for ix=1:nx-1
%   for jy=1:ny-1
%     for kz=1:nz
%       Hz(ix, jy, kz) = Hz(ix, jy, kz) + dt/(dx*dy) * (...
%         + dx*(Ex(ix, jy+1, kz) - Ex(ix, jy, kz))...
%         - dy*(Ey(ix+1, jy, kz) - Ey(ix, jy, kz)) ) ;
%     end
%   end
% end

end