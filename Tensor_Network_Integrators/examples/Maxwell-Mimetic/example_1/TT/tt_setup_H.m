function [ Hx, Hy, Hz ] = tt_setup_H(time, dx, dy, dz, nx, ny, nz )

% setup Hx,Hy,Hz in tt format


% setup Hx,
xgrid = dx*([1:nx] -1);
ygrid = dy*([1:ny-1] - 1/2);
zgrid = dz*([1:nz-1]- 1/2);
Hx = tt_Hx_value(xgrid, ygrid, zgrid, time);


% setup Hy
xgrid = dx*([1:nx-1] -1/2);
ygrid = dy*([1:ny] - 1);
zgrid = dz*([1:nz-1]- 1/2);
Hy = tt_Hy_value(xgrid, ygrid, zgrid, time);

% setup Hz

xgrid = dx*([1:nx-1] -1/2);
ygrid = dy*([1:ny-1] - 1/2);
zgrid = dz*([1:nz]- 1);
Hz = tt_Hz_value(xgrid, ygrid, zgrid, time);
end 