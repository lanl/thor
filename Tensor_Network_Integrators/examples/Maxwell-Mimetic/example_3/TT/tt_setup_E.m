function [ Ex, Ey, Ez ] = tt_setup_E( time, dx, dy, dz, nx, ny, nz, tol )

% setup Ex,
xgrid = dx*([1:nx-1] -1/2);
ygrid = dy*([1:ny] - 1);
zgrid = dz*([1:nz]- 1);
Ex = tt_Ex_value(xgrid, ygrid, zgrid, time);


%setup Ey
xgrid = dx*([1:nx] -1);
ygrid = dy*([1:ny-1] - 1/2);
zgrid = dz*([1:nz]- 1);
Ey = tt_Ey_value(xgrid, ygrid, zgrid, time);

%setup Ez
xgrid = dx*([1:nx] -1);
ygrid = dy*([1:ny] - 1);
zgrid = dz*([1:nz-1]- 1/2);
Ez = tt_Ez_value(xgrid, ygrid, zgrid, time, tol);


end