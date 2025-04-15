function [ Hz ] = update_Hz( Hz, Ex, Ey, dx, dy, dt, nx, ny, nz )

for ix=1:nx-1
  for jy=1:ny-1
    for kz=1:nz
      Hz(ix, jy, kz) = Hz(ix, jy, kz) + dt/(dx*dy) * (...
        + dx*(Ex(ix, jy+1, kz) - Ex(ix, jy, kz))...
        - dy*(Ey(ix+1, jy, kz) - Ey(ix, jy, kz)) ) ;
    end
  end
end

end