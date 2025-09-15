function [ tHx, tHy, tHz ] = tt_update_DBC_H(tHx, tHy, tHz, dx, dy, dz, nx, ny, nz, time, tol)

% -----------
% HBC at "x0"
% -----------

%% update tHx at ix = 1
ix = 1;
x  = dx*(ix-1);

%compute the Hx_value_slice
xgrid = [x];
ygrid = dy*([1:ny-1]-1/2);
zgrid = dz*([1:nz-1]-1/2);
tempnew = tt_Hx_value(xgrid,ygrid,zgrid,time);

%expand the slice - the new boundary
tempnew = expand3d(tempnew,[1,ix,nx]);

% remove the old boundary from tHx
tHx{1}(:,ix,:) = 0;
% 
% % get tHx old slice
% tempold = get_tt_slice(tHx,[1,ix,ix]);

% update tHx
tHx = round(tHx + tempnew, tol); 
%%

% -----------
% HBC at "x1"
% -----------


%% update tHx at ix = nx
ix = nx;
x = dx*(ix-1);
%compute the Hx_value_slice
xgrid = [x];
tempnew = tt_Hx_value(xgrid,ygrid,zgrid,time);
%expand the slice
tempnew = expand3d(tempnew,[1,ix,nx]);
% get tHx old slice
tempold = get_tt_slice(tHx,[1,ix,ix]);
% update tHx
tHx = round(tHx - tempold + tempnew, tol); 

%%
% -----------
% HBC at "y0"
% -----------

%% update tHy at jy = 1
jy = 1;
y  = dy*(jy-1);

%compute the Hx_value_slice
xgrid = dx*([1:nx-1]-1/2);
ygrid = [y];
zgrid = dz*([1:nz-1]-1/2);
tempnew = tt_Hy_value(xgrid,ygrid,zgrid,time);
%expand the slice
tempnew = expand3d(tempnew,[2,jy,ny]);
% get tHx old slice
tempold = get_tt_slice(tHy,[2,jy,jy]);
% update tHx
tHy = round(tHy - tempold + tempnew, tol); 

%%
% -----------
% HBC at "y1"
% -----------
%% update tHy at jy = ny
jy = ny;
y  = dy*(jy-1);
tol = 1e-15;
%compute the Hx_value_slice
ygrid = [y];
tempnew = tt_Hy_value(xgrid,ygrid,zgrid,time);
%expand the slice
tempnew = expand3d(tempnew,[2,jy,ny]);
% get tHx old slice
tempold = get_tt_slice(tHy,[2,jy,jy]);
% update tHx
tHy = round(tHy - tempold + tempnew, tol); 
%%
% --------------------
% HBC at "z0" (bottom)
% --------------------
%% update tHz at kz = 1
kz = 1;
z  = dz*(kz-1);

%compute the Hx_value_slice
xgrid = dx*([1:nx-1]-1/2);
ygrid = dy*([1:ny-1]-1/2);
zgrid = [z];
tempnew = tt_Hz_value(xgrid,ygrid,zgrid,time);

%expand the slice
tempnew = expand3d(tempnew,[3,kz,nz]);
% get tHx old slice
tempold = get_tt_slice(tHz,[3,kz,kz]);
% update tHx
tHz = round(tHz - tempold + tempnew, tol); 
%%
% -----------------
% HBC at "z1" (top)
% -----------------

%% update tHz at kz=nz
kz = nz;
z  = dz*(kz-1);
%compute the Hx_value_slice
zgrid = [z];
tempnew = tt_Hz_value(xgrid,ygrid,zgrid,time);
%expand the slice
tempnew = expand3d(tempnew,[3,kz,nz]);
% get tHx old slice
tempold = get_tt_slice(tHz,[3,kz,kz]);
% update tHx
tHz = round(tHz - tempold + tempnew, tol);     

end