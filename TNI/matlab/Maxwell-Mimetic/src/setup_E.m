function [ Ex, Ey, Ez ] = setup_E( time, dx, dy, dz, nx, ny, nz )

% setup
Ex = zeros(nx-1,ny,nz);
Ey = zeros(nx,ny-1,nz);
Ez = zeros(nx,ny,nz-1);

% set Ex
for ix=1:nx-1
  for jy=1:ny
    for kz=1:nz
      x = dx*(ix-1)+dx/2;
      y = dy*(jy-1);
      z = dz*(kz-1);
      Ex(ix, jy, kz) = Ex_value( x, y, z, time );
    end
  end
end

% set Ey
for ix=1:nx
  for jy=1:ny-1
    for kz=1:nz
      x = dx*(ix-1);
      y = dy*(jy-1)+dy/2;
      z = dz*(kz-1);
      Ey(ix, jy, kz) = Ey_value( x, y, z, time );
    end
  end
end
% setEz
for ix=1:nx
  for jy=1:ny
    for kz=1:nz-1
      x = dx*(ix-1);
      y = dy*(jy-1);
      z = dz*(kz-1)+dz/2;
      Ez(ix, jy, kz) = Ez_value( x, y, z, time );
    end
  end
end

end