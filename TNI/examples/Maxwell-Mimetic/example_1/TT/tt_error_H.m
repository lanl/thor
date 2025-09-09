function [err_H, nrm_H] =  tt_error_H(tHx, tHy, tHz, dx, dy, dz, nx, ny, nz, time, tol)

%set time
t = time;

%% tHx
xgrid = dx*([1:nx]-1);
ygrid = dy*([1:ny-1]-1/2);
zgrid = dz*([1:nz-1]-1/2);

exact_tHx = tt_Hx_value(xgrid,ygrid,zgrid,t);
%absolute error
err_Hx = sum(round((tHx-exact_tHx).^2,tol));
err_Hx = sqrt( err_Hx / ( nx*(ny-1)*(nz-1) ) );

%true solution norm
nrm_Hx = sum(round(exact_tHx.^2,tol));
nrm_Hx = sqrt( nrm_Hx / ( nx*(ny-1)*(nz-1) ) );

%% tHy
xgrid = dx*([1:nx-1]-1/2);
ygrid = dy*([1:ny] -1);
zgrid = dz*([1:nz-1]-1/2);

exact_tHy = tt_Hy_value(xgrid,ygrid,zgrid,t);
%absolute error
err_Hy = sum(round((tHy-exact_tHy).^2,tol));
err_Hy = sqrt( err_Hy / ( (nx-1)*ny*(nz-1) ) );


%true solution norm
nrm_Hy = sum(round(exact_tHy.^2,tol));
nrm_Hy = sqrt( nrm_Hy / ( (nx-1)*ny*(nz-1) ) );

%% tHz
xgrid = dx*([1:nx-1]-1/2);
ygrid = dy*([1:ny-1]-1/2);
zgrid = dz*([1:nz]-1);

exact_tHz = tt_Hz_value(xgrid,ygrid,zgrid,t);
%absolute error
err_Hz = sum(round((tHz-exact_tHz).^2,tol));
err_Hz = sqrt( err_Hz / ( (nx-1)*(ny-1)*nz ) );

%true solution norm
nrm_Hz = sum(round(exact_tHz.^2,tol));
nrm_Hz = sqrt( nrm_Hz / ( (nx-1)*(ny-1)*nz ) );

%% final error
err_H = err_Hx + err_Hy + err_Hz;
nrm_H = nrm_Hx + nrm_Hy + nrm_Hz;
