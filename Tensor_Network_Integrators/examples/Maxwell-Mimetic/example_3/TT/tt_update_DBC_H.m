function [ tHx, tHy, tHz ] = tt_update_DBC_H(tHx, tHy, tHz, dx, dy, dz, nx, ny, nz, time, tol)

%% update tHx at ix = [1,nx]

%compute the Hx_value_slice
xgrid = dx*([1,nx]-1);
ygrid = dy*([1:ny-1]-1/2);
zgrid = dz*([1:nz-1]-1/2);
tempnew = tt_Hx_value(xgrid,ygrid,zgrid,time, tol);

%expand the boundary
r = tempnew.r;
Gi = zeros(r(1),nx,r(2));
Gi(:,1,:) = tempnew{1}(:,1,:);%slice 1
Gi(:,nx,:) = tempnew{1}(:,2,:);%slice nx
tempnew{1} = Gi;

% set the old boundary of tHx to zeros
tHx = tt_1D_set_zeros(tHx,1,[1,nx]);

% update tHx
tHx = round(tHx + tempnew, tol); 


%% update tHy at jy = [1, ny]

%compute the Hx_value_slice
xgrid = dx*([1:nx-1]-1/2);
ygrid = dy*([1 ny]-1);
zgrid = dz*([1:nz-1]-1/2);
tempnew = tt_Hy_value(xgrid,ygrid,zgrid,time, tol);
%expand the slice
r = tempnew.r;
Gi = zeros(r(2),ny,r(3)); %second dimension
Gi(:,1,:) = tempnew{2}(:,1,:);%slice 1
Gi(:,ny,:) = tempnew{2}(:,2,:);%slice nx
tempnew{2} = Gi;

%set tHy boundary in the y dim to zeros
tHy = tt_1D_set_zeros(tHy,2,[1,ny]);

% update tHx
tHy = round(tHy + tempnew, tol); 


% --------------------
% HBC at "z0" (bottom)
% --------------------
%% update tHz at kz = [1, ny]

%compute the Hx_value_slice
xgrid = dx*([1:nx-1]-1/2);
ygrid = dy*([1:ny-1]-1/2);
zgrid = dz*([1, nz]-1);
tempnew = tt_Hz_value(xgrid,ygrid,zgrid,time);

%expand the slice
r = tempnew.r;
Gi = zeros(r(3),nz,r(4)); %second dimension
Gi(:,1,:) = tempnew{3}(:,1,:);%slice 1
Gi(:,nz,:) = tempnew{3}(:,2,:);%slice nz
tempnew{3} = Gi;

%set boundary to zeros
tHz = tt_1D_set_zeros(tHz,3,[1,nz]);

% update tHx
tHz = round(tHz + tempnew, tol); 

end