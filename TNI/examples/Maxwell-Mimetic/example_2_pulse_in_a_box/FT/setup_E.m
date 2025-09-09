function [ Ex, Ey, Ez ] = setup_E( time, dx, dy, dz, nx, ny, nz, X0 )

% setup
Ex = zeros(nx-1,ny,nz);
Ey = zeros(nx,ny-1,nz);
Ez = zeros(nx,ny,nz-1);
b = 25; %beta value
% setEz
for ix=1:nx
  for jy=1:ny
    for kz=1:nz-1
      x = X0 + dx*(ix-1);
      y = X0 + dy*(jy-1);
      z = X0 + dz*(kz-1)+dz/2;
      Ez(ix, jy, kz) = exp(-b*(x^2+y^2));
    end
  end
end

end