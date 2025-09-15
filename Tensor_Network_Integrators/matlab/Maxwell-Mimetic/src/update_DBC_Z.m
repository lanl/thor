function [ Ex, Ey ] = update_DBC_Z( Ex, Ey, dx, dy, dz, nx, ny, nz, time )

% --------------------
% HBC at "z0" (bottom)
% --------------------

kz = 1;

% Ex
for ix=1:nx-1
  for jy=1:ny
    x = dx*(ix-1)+dx/2;
    y = dy*(jy-1);
    z = dz*(kz-1);
    Ex  (ix, jy, kz) = Ex_value( x, y, z, time );
  end
end

% Ey
for ix=1:nx
  for jy=1:ny-1
    x = dx*(ix-1);
    y = dy*(jy-1)+dy/2;
    z = dz*(kz-1);
    Ey  (ix, jy, kz) = Ey_value( x, y, z, time );
  end
end

% -----------------
% HBC at "z1" (top)
% -----------------

kz = nz;

% Ex
for ix=1:nx-1
  for jy=1:ny
    x = dx*(ix-1)+dx/2;
    y = dy*(jy-1);
    z = dz*(kz-1);
    Ex  (ix, jy, kz) = Ex_value( x, y, z, time );
  end
end

% Ey
for ix=1:nx
  for jy=1:ny-1
    x = dx*(ix-1);
    y = dy*(jy-1)+dy/2;
    z = dz*(kz-1);
    Ey  (ix, jy, kz) = Ey_value( x, y, z, time );
  end
end

end
