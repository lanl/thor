function [ Hx ] = update_Hx( Hx, Ey, Ez, dy, dz, dt, nx, ny, nz )

for ix=1:nx
  for jy=1:ny-1
    for kz=1:nz-1
      Hx(ix, jy, kz) = Hx(ix, jy, kz) + dt/(dy*dz) * (...
        + dy*(Ey(ix, jy, kz+1) - Ey(ix, jy, kz))...
        - dz*(Ez(ix, jy+1, kz) - Ez(ix, jy, kz)) ) ;
%           Hx(ix, jy, kz) = Hx(ix, jy, kz) + (...
%         + dt/dz*(Ey(ix, jy, kz+1) - Ey(ix, jy, kz))...
%         - dt/dy*(Ez(ix, jy+1, kz) - Ez(ix, jy, kz)) ) ;
    end
  end
end

end