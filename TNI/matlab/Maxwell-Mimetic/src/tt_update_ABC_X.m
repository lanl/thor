function [ tEy, tEz, tEy_x0, tEy_x1, tEz_x0, tEz_x1] = ...
  tt_update_ABC_X( tEy, tEz, tEy_x0, tEy_x1, tEz_x0, tEz_x1, nx, ny, nz, abccoef,tol)


% -----------
% ABC at "x0"
% -----------

ix = 1;
[tEy, tEy_x0] = update_E_slice(tEy, tEy_x0,ix, nx, ny-1, nz, abccoef,tol);
[tEz, tEz_x0] = update_E_slice(tEz, tEz_x0,ix, nx, ny, nz-1, abccoef,tol);
% % Ey
% for jy=1:ny-1
%   for kz=1:nz
%     Ey(ix, jy, kz) = Ey_x0(jy, kz) + abccoef * ( Ey(ix+1, jy, kz) - Ey(ix, jy, kz) );
%     Ey_x0 (jy, kz) = Ey(ix+1, jy, kz);
%   end
% end
% 
% % Ez
% for jy=1:ny
%   for kz=1:nz-1
%     Ez(ix, jy, kz) = Ez_x0(jy, kz) + abccoef * ( Ez(ix+1, jy, kz) - Ez(ix, jy, kz) );
%     Ez_x0 (jy, kz) = Ez(ix+1, jy, kz);
%   end
% end

% -----------
% ABC at "x1"
% -----------


ix = nx;
[tEy, tEy_x1] = update_E_slice(tEy, tEy_x1,ix, nx, ny-1, nz, abccoef,tol);
[tEz, tEz_x1] = update_E_slice(tEz, tEz_x1,ix, nx, ny, nz-1, abccoef,tol);

% ix = nx;
% 
% % Ey
% for jy=1:ny-1
%   for kz=1:nz
%     Ey(ix, jy, kz) = Ey_x1(jy, kz) + abccoef * (Ey(ix-1, jy, kz) - Ey(ix, jy, kz));
%     Ey_x1 (jy, kz) = Ey(ix-1, jy, kz);
%   end
% end
% 
% %Ez
% for jy=1:ny
%   for kz=1:nz-1
%     Ez(ix, jy, kz) = Ez_x1(jy, kz) + abccoef * (Ez(ix-1, jy, kz) - Ez(ix, jy, kz));
%     Ez_x1 (jy, kz) = Ez(ix-1, jy, kz);
%   end
% end
% 
% end
end
%% subroutines
function [tEy, tEy_x0] = update_E_slice(tEy, tEy_x0,ix, nx, ny, nz, abccoef,tol)

tEy_x0 = expand3d(tEy_x0,[1,ix, nx]);
% get Ey(ix+1)
if ix == 1
  FD = get_tt_slice(tEy,[1,ix+1,ix]);
elseif ix==nx
  FD = get_tt_slice(tEy,[1,ix-1,ix]);
end
%gett Ey(ix)
BD = get_tt_slice(tEy,[1, ix, ix]);

tEy = round(tEy-BD,tol) + tEy_x0  + abccoef * round(FD - BD,tol);
tEy = round(tEy,tol);
tEy_x0 = FD;
%Is needed here only for testing purposes at the moment 
tEy_x0 = project_to_2d(tEy_x0,ix,nx,ny,nz);
end















