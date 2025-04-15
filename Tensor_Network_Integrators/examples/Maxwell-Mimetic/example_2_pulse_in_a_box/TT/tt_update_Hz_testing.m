function [ tHz] = tt_update_Hz_testing( tHz, tEx, tEy, dx, dy, dt, tol)

% compute dExdy
dExdy = tt_FD(tEx,2);
dExdy = round(dExdy,tol);
% ll(1) = fprintf('max rank of dExdy = %d \n', max(dExdy.r));
rExdy = max(dExdy.r);

%compute dEydx
dEydx = tt_FD(tEy,1);
dEydx = round(dEydx,tol);
% ll(2) = fprintf('max rank of dEydx = %d \n', max(dEydx.r));

% 
Dterm = round((dt/dy)*dExdy - (dt/dx)*dEydx, tol);
% ll(3) = fprintf('max rank of Dterm = %d \n', max(Dterm.r));
% ll(4) = fprintf('max rank of dHz before update = %d \n', max(tHz.r));
tHz = round( tHz + Dterm,tol);
% ll(5) = fprintf('max rank of dHz after update = %d \n', max(tHz.r));
% lls = sum(ll);
% tHz= round(tHz, tol);

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