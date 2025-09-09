function [ Hx, Hy, Hz ] = update_DBC_H( Hx, Hy, Hz, dx, dy, dz, nx, ny, nz, time )

% -----------
% HBC at "x0"
% -----------

ix = 1;

    
% Hx
for jy=1:ny-1
  for kz=1:nz-1
    Hx(ix, jy, kz) = 0;
  end
end

% -----------
% HBC at "x1"
% -----------

ix = nx;


% Hx
for jy=1:ny-1
  for kz=1:nz-1
    Hx(ix, jy, kz) = 0;
  end
end

% -----------
% HBC at "y0"
% -----------

jy = 1;
    
% Hy
for ix=1:nx-1
  for kz=1:nz-1
    Hy  (ix, jy, kz) = 0;
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
    Hy  (ix, jy, kz) = 0;
  end
end


% --------------------
% HBC at "z0" (bottom)
% --------------------

kz = 1;

% Hz
for ix=1:nx-1
  for jy=1:ny-1
    Hz  (ix, jy, kz) = 0;
  end
end

% -----------------
% HBC at "z1" (top)
% -----------------

kz = nz;
    
% Hz
for ix=1:nx-1
  for jy=1:ny-1
    Hz  (ix, jy, kz) = 0;
  end
end

end