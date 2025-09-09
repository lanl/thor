function [ Hx, Hy, Hz ] = setup_H( time, dx, dy, dz, nx, ny, nz )

% setup
Hx   = zeros(nx,ny-1,nz-1);
Hy   = zeros(nx-1,ny,nz-1);
Hz   = zeros(nx-1,ny-1,nz);

% set magnetic coeffs
for ix=1:nx
  for jy=1:ny-1
    for kz=1:nz-1
      x = dx*(ix-1);
      y = dy*(jy-1)+dy/2;
      z = dz*(kz-1)+dz/2;
      Hx  (ix, jy, kz) = Hx_value( x, y, z, time );
    end
  end
end

% set magnetic coeffs
for ix=1:nx-1
  for jy=1:ny
    for kz=1:nz-1
      x = dx*(ix-1)+dx/2;
      y = dy*(jy-1);
      z = dz*(kz-1)+dz/2;
      Hy  (ix, jy, kz) = Hy_value( x, y, z, time );
    end
  end
end


% set magnetic coeffs
for ix=1:nx-1
  for jy=1:ny-1
    for kz=1:nz
      x = dx*(ix-1)+dx/2;
      y = dy*(jy-1)+dy/2;
      z = dz*(kz-1);
      Hz  (ix, jy, kz) = Hz_value( x, y, z, time ); 
    end
  end
end


end 