function [ Ex ] = update_Ex( Ex, Jx, Hy, Hz, dy, dz, dt, nx, ny, nz )

for ix=1:nx-1
  for jy=2:ny-1
    for kz=2:nz-1
      Ex(ix, jy, kz) = Ex(ix, jy, kz) + dt/(dy*dz) * (...
        + dz*(Hz(ix, jy, kz) - Hz(ix, jy-1, kz))...
        - dy*(Hy(ix, jy, kz) - Hy(ix, jy, kz-1))...
        ) + dt * Jx(ix, jy, kz) ;
    end
  end
end

end