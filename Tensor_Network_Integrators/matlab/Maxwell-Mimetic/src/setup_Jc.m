function [ Jx, Jy, Jz ] = setup_Jc( nx, ny, nz )

% setup
Jx = zeros(nx-1,ny,nz);
Jy = zeros(nx,ny-1,nz);
Jz = zeros(nx,ny,nz-1);

end