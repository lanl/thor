function [ Hx, Hy, Hz ] = update_DBC_H( Hx, Hy, Hz, dx, dy, dz, nx, ny, nz, time )

% -----------
% HBC at "x0"
% -----------

ix = 1;
x  = dx*(ix-1);
    
% Hx
for jy=1:ny-1
  for kz=1:nz-1
    y = dy*(jy-1)+dy/2;
    z = dz*(kz-1)+dz/2;
    Hx(ix, jy, kz) = Hx_value( x, y, z, time );
  end
end

% -----------
% HBC at "x1"
% -----------

ix = nx;
x = dx*(ix-1);

% Hx
for jy=1:ny-1
  for kz=1:nz-1
    y = dy*(jy-1)+dy/2;
    z = dz*(kz-1)+dz/2;
    Hx(ix, jy, kz) = Hx_value( x, y, z, time );
  end
end

% -----------
% HBC at "y0"
% -----------

jy = 1;
y  = dy*(jy-1);
    
% Hy
for ix=1:nx-1
  for kz=1:nz-1
    x = dx*(ix-1)+dx/2;
    z = dz*(kz-1)+dz/2;
    Hy  (ix, jy, kz) = Hy_value( x, y, z, time );
  end
end

% -----------
% HBC at "y1"
% -----------

jy = ny;
y  = dy*(jy-1);

% Hy
for ix=1:nx-1
  for kz=1:nz-1
    x = dx*(ix-1)+dx/2;
    z = dz*(kz-1)+dz/2;
    Hy  (ix, jy, kz) = Hy_value( x, y, z, time );
  end
end


% --------------------
% HBC at "z0" (bottom)
% --------------------

kz = 1;
z  = dz*(kz-1);
    
% Hz
for ix=1:nx-1
  for jy=1:ny-1
    x = dx*(ix-1)+dx/2;
    y = dy*(jy-1)+dy/2;
    Hz  (ix, jy, kz) = Hz_value( x, y, z, time );
  end
end

% -----------------
% HBC at "z1" (top)
% -----------------

kz = nz;
z  = dz*(kz-1);
    
% Hz
for ix=1:nx-1
  for jy=1:ny-1
    x = dx*(ix-1)+dx/2;
    y = dy*(jy-1)+dy/2;
    Hz  (ix, jy, kz) = Hz_value( x, y, z, time );
  end
end

end