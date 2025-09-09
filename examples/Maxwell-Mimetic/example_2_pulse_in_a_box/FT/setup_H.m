function [ Hx, Hy, Hz ] = setup_H( time, dx, dy, dz, nx, ny, nz )

% setup
Hx   = zeros(nx,ny-1,nz-1);
Hy   = zeros(nx-1,ny,nz-1);
Hz   = zeros(nx-1,ny-1,nz);

end