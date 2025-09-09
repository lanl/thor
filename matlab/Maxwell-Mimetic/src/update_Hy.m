function [ Hy ] = update_Hy( Hy, Ez, Ex, dz, dx, dt, nx, ny, nz  )

for ix=1:nx-1
  for jy=1:ny
    for kz=1:nz-1
      Hy(ix, jy, kz) = Hy(ix, jy, kz) + dt/(dz*dx) * (...
        + dz*(Ez(ix+1, jy, kz) - Ez(ix, jy, kz))...
        - dx*(Ex(ix, jy, kz+1) - Ex(ix, jy, kz)) ) ;
    end
  end
end

end
