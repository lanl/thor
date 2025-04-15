function [ Ez, Ex ] = update_DBC_Y( Ez, Ex, dx, dy, dz, nx, ny, nz, time )

% -----------
% HBC at "y0"
% -----------

jy = 1;

% Ex
for ix=1:nx-1
  for kz=1:nz
    x = dx*(ix-1)+dx/2;
    y = dy*(jy-1);
    z = dz*(kz-1);
    Ex(ix, jy, kz) = 0;
  end
end

% Ez
for ix=1:nx
  for kz=1:nz-1
    x = dx*(ix-1);
    y = dy*(jy-1);
    z = dz*(kz-1)+dz/2;
    Ez  (ix, jy, kz) = 0;
  end
end

% -----------
% HBC at "y1"
% -----------

jy = ny;

% Ex
for ix=1:nx-1
  for kz=1:nz
    x = dx*(ix-1)+dx/2;
    y = dy*(jy-1);
    z = dz*(kz-1);
    Ex(ix, jy, kz) = 0;
  end
end

% Ez
for ix=1:nx
  for kz=1:nz-1
    x = dx*(ix-1);
    y = dy*(jy-1);
    z = dz*(kz-1)+dz/2;
    Ez  (ix, jy, kz) = 0;
  end
end

end