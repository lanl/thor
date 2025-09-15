function [err_E, nrm_E] =  tt_error_E(tEx, tEy, tEz, dx, dy, dz, nx, ny, nz, time, tol)

%set time
t = time;

%% tHx
xgrid = dx*([1:nx-1]-1/2);
ygrid = dy*([1:ny]-1);
zgrid = dz*([1:nz]-1);

exact_tEx = tt_Ex_value(xgrid,ygrid,zgrid,t);
%absolute error
err_Ex = sum(round((tEx-exact_tEx).^2,tol));
err_Ex = sqrt( err_Ex/((nx-1)*ny*nz));

%true solution norm
nrm_Ex = sum(round(exact_tEx.^2,tol));
nrm_Ex = sqrt( nrm_Ex / ( (nx-1)*ny*nz ) );

%% tHy
xgrid = dx*([1:nx]-1);
ygrid = dy*([1:ny-1]-1/2);
zgrid = dz*([1:nz]-1);

exact_tEy = tt_Ey_value(xgrid,ygrid,zgrid,t);
%absolute error
err_Ey = sum(round((tEy-exact_tEy).^2,tol));
err_Ey = sqrt( err_Ey / ( nx*(ny-1)*nz ) );


%true solution norm
nrm_Ey = sum(round(exact_tEy.^2,tol));
nrm_Ey = sqrt( nrm_Ey / ( nx*(ny-1)*nz ) );


%% tHz
xgrid = dx*([1:nx]-1);
ygrid = dy*([1:ny]-1);
zgrid = dz*([1:nz-1]-1/2);

exact_tEz = tt_Ez_value(xgrid,ygrid,zgrid,t);
%absolute error
err_Ez = sum(round((tEz-exact_tEz).^2,tol));
err_Ez = sqrt( err_Ez / ( nx*ny*(nz-1) ) );

%true solution norm
nrm_Ez = sum(round(exact_tEz.^2,tol));
nrm_Ez = sqrt( nrm_Ez / ( nx*ny*(nz-1) ) );

%% final error
err_E = err_Ex + err_Ey + err_Ez;
nrm_E = nrm_Ex + nrm_Ey + nrm_Ez;
