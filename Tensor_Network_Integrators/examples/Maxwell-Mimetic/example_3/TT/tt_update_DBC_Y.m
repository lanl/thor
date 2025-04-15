function [ tEz, tEx ] = tt_update_DBC_Y( tEz, tEx, dx, dy, dz, nx, ny, nz, time, tol)

% -----------
% HBC at "y0"
% -----------
%% update tEx at jy = 1
jy = [1 ny];
xgrid = dx*([1:(nx-1)]-1/2);
ygrid = dy*(jy-1);
zgrid = dz*([1:nz]-1);

tempnew = tt_Ex_value(xgrid,ygrid,zgrid,time);

%expand the boundary
r = tempnew.r;
Gi = zeros(r(1),ny,r(2));
Gi(:,1,:) = tempnew{2}(:,1,:);%slice 1
Gi(:,ny,:) = tempnew{2}(:,2,:);%slice nx
tempnew{2} = Gi;

% set the old boundary of tHx to zeros
tEx = tt_1D_set_zeros(tEx,2,[1,ny]);

% update tHx
tEx = round(tEx + tempnew, tol); 


%% update tEz at jy=1

xgrid = dx*([1:nx]-1);
ygrid = dy*(jy-1);
zgrid = dz*([1:nz-1]-1/2);

tempnew = tt_Ez_value(xgrid,ygrid,zgrid,time, tol);
%expand the boundary
r = tempnew.r;
Gi = zeros(r(1),ny,r(2));
Gi(:,1,:) = tempnew{2}(:,1,:);%slice 1
Gi(:,ny,:) = tempnew{2}(:,2,:);%slice nx
tempnew{2} = Gi;

% set the old boundary of tHx to zeros
tEz = tt_1D_set_zeros(tEz,2,[1,ny]);

% update tHx
tEz = round(tEz + tempnew, tol); 

end