function [ Ez ] = update_Ez( Ez, Jz, Hx, Hy, dx, dy, dt, nx, ny, nz )

for ix=2:nx-1
  for jy=2:ny-1
    for kz=1:nz-1
      Ez(ix, jy, kz) = Ez(ix, jy, kz) + dt/(dx*dy) * (...
        + dy*(Hy(ix, jy, kz) - Hy(ix-1, jy, kz))...
        - dx*(Hx(ix, jy, kz) - Hx(ix, jy-1, kz))...
        ) + dt * Jz(ix, jy, kz);
    end
  end
end

end