function [ tEx, tEy ] = tt_update_DBC_Z( tEx, tEy, dx, dy, dz, nx, ny, nz, time, tol )

% --------------------
% HBC at "z0" (bottom)
% --------------------
%% update Ex at kz =1
kz = [1,nz];
xgrid = dx*([1:(nx-1)]-1/2);
ygrid = dy*([1:ny]-1);
zgrid = dz*(kz-1);

tempnew = tt_Ex_value(xgrid,ygrid,zgrid,time);

%expand the boundary
r = tempnew.r;
Gi = zeros(r(1),nz,r(2));
Gi(:,1,:) = tempnew{3}(:,1,:);%slice 1
Gi(:,nz,:) = tempnew{3}(:,2,:);%slice nx
tempnew{3} = Gi;

% set the old boundary of tHx to zeros
tEx = tt_1D_set_zeros(tEx,3,[1,nz]);

% update tHx
tEx = round(tEx + tempnew, tol); 


%% update Ey at kz = 1
xgrid = dx*([1:nx]-1);
ygrid = dy*([1:ny-1]-1/2);
zgrid = dz*(kz-1);

tempnew = tt_Ey_value(xgrid,ygrid,zgrid,time);

%expand the boundary
r = tempnew.r;
Gi = zeros(r(1),nz,r(2));
Gi(:,1,:) = tempnew{3}(:,1,:);%slice 1
Gi(:,nz,:) = tempnew{3}(:,2,:);%slice nx
tempnew{3} = Gi;

% set the old boundary of tHx to zeros
tEy = tt_1D_set_zeros(tEy,3,[1,nz]);

% update tHx
tEy = round(tEy + tempnew, tol); 

end
