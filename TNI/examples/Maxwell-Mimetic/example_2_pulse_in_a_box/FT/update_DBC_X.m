function [ Ey, Ez ] = update_DBC_X( Ey, Ez, dx, dy, dz, nx, ny, nz, time )

% -----------
% HBC at "x0"
% -----------

ix = 1;

% Ey
for jy=1:ny-1
  for kz=1:nz
    x = dx*(ix-1);
    y = dy*(jy-1)+dy/2;
    z = dz*(kz-1);
    Ey(ix, jy, kz) = 0;
  end
end

% Ez
for jy=1:ny
  for kz=1:nz-1
    x = dx*(ix-1);
    y = dy*(jy-1);
    z = dz*(kz-1)+dz/2;
    Ez(ix, jy, kz) = 0;
  end
end

% -----------
% HBC at "x1"
% -----------

ix = nx;

% Ey
for jy=1:ny-1
  for kz=1:nz
    x = dx*(ix-1);
    y = dy*(jy-1)+dy/2;
    z = dz*(kz-1);
    Ey  (ix, jy, kz) = 0;
  end
end

% Ez
for jy=1:ny
  for kz=1:nz-1
    x = dx*(ix-1);
    y = dy*(jy-1);
    z = dz*(kz-1)+dz/2;
    Ez  (ix, jy, kz) = 0;
  end
end

end