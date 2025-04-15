function [ Ey ] = update_Ey( Ey, Jy, Hz, Hx, dz, dx, dt, nx, ny, nz )

for ix=2:nx-1
  for jy=1:ny-1
    for kz=2:nz-1
      Ey(ix, jy, kz) = Ey(ix, jy, kz) + dt/(dz*dx) * (...
        + dx*(Hx(ix, jy, kz) - Hx(ix, jy, kz-1))...
        - dz*(Hz(ix, jy, kz) - Hz(ix-1, jy, kz))...
        ) + dt * Jy(ix, jy, kz);
    end
  end
end

end