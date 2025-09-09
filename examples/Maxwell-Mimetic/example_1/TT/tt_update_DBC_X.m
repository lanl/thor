function [ tEy, tEz ] = tt_update_DBC_X( tEy, tEz, dx, dy, dz, nx, ny, nz, time, tol)

% -----------
% HBC at "x0"
% -----------

%% update tEy at ix = 1,nx
ix = [1 nx];
xgrid = [dx*(ix-1)];
ygrid = dy*([1:ny-1]-1/2);
zgrid = dz*([1:nz]-1);

tempnew = tt_Ey_value(xgrid,ygrid,zgrid,time);
%expand the boundary
r = tempnew.r;
Gi = zeros(r(1),nx,r(2));
Gi(:,1,:) = tempnew{1}(:,1,:);%slice 1
Gi(:,nx,:) = tempnew{1}(:,2,:);%slice nx
tempnew{1} = Gi;

% set the old boundary of tHx to zeros
tEy = tt_1D_set_zeros(tEy,1,[1,nx]);

% update tHx
tEy = round(tEy + tempnew, tol); 

%% update tEz at ix = 1,nx

xgrid = [dx*(ix-1)];
ygrid = dy*([1:ny]-1);
zgrid = dz*([1:nz-1]-1/2);

tempnew = tt_Ez_value(xgrid,ygrid,zgrid,time);

%expand the boundary
r = tempnew.r;
Gi = zeros(r(1),nx,r(2));
Gi(:,1,:) = tempnew{1}(:,1,:);%slice 1
Gi(:,nx,:) = tempnew{1}(:,2,:);%slice nx
tempnew{1} = Gi;

% set the old boundary of tHx to zeros
tEz = tt_1D_set_zeros(tEz,1,[1,nx]);

% update tHx
tEz = round(tEz + tempnew, tol); 


end