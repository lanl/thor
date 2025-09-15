function[ Jx, Jy, Jz ] = update_Jc( Jx, Jy, Jz, dx, dy, dz, dt, nx, ny, nz, time )

% set time
tm=time;
te=time+dt/2;

%% set Jx
for ix=1:nx-1
  for jy=1:ny
    for kz=1:nz
      
      x = dx*(ix-1);
      y = dy*(jy-1);
      z = dz*(kz-1);
      
      Ex_1 = Ex_value( x+dx/2, y,         z,         te    );
      Ex_0 = Ex_value( x+dx/2, y,         z,         te-dt );
      
      Hz_1 = Hz_value( x+dx/2, y+dy/2,    z,         tm    );
      Hz_0 = Hz_value( x+dx/2, y+dy/2-dy, z,         tm    );
      
      Hy_1 = Hy_value( x+dx/2, y,         z+dz/2,    tm    );
      Hy_0 = Hy_value( x+dx/2, y,         z+dz/2-dz, tm    );
      
      Jx(ix, jy, kz) = ( Ex_1 - Ex_0 )/dt - 1/(dy*dz) * ( dz*(Hz_1 - Hz_0) - dy*(Hy_1 - Hy_0) );
      
      %  Ex(ix, jy, kz) = Ex(ix, jy, kz) + dt/(dy*dz) * (...
      %         + dz*(Hz(ix, jy, kz) - Hz(ix, jy-1, kz))...
      %         - dy*(Hy(ix, jy, kz) - Hy(ix, jy, kz-1))...
      %  )
    end
  end
end

%% set Jy
for ix=1:nx
  for jy=1:ny-1
    for kz=1:nz
      
      x = dx*(ix-1);
      y = dy*(jy-1);
      z = dz*(kz-1);
      
      Ey_1 = Ey_value( x,         y+dy/2, z,         te    );
      Ey_0 = Ey_value( x,         y+dy/2, z,         te-dt );
      
      Hx_1 = Hx_value( x,         y+dy/2, z+dz/2,    tm    );
      Hx_0 = Hx_value( x,         y+dy/2, z+dz/2-dz, tm    );
      
      Hz_1 = Hz_value( x+dx/2,    y+dy/2, z,         tm    );
      Hz_0 = Hz_value( x+dx/2-dx, y+dy/2, z,         tm    );
      
      Jy(ix, jy, kz) = ( Ey_1 - Ey_0 ) / dt - 1/(dz*dx) * ( dx*(Hx_1 - Hx_0) - dz*(Hz_1 - Hz_0) );
      
      %  Ey(ix, jy, kz) = Ey(ix, jy, kz) + dt/(dz*dx) * (...
      %         + dx*(Hx(ix, jy, kz) - Hx(ix, jy, kz-1))...
      %         - dz*(Hz(ix, jy, kz) - Hz(ix-1, jy, kz))...
      %  ) + dt * Jy(ix, jy, kz);
    end
  end
end

% set Jz
for ix=1:nx
  for jy=1:ny
    for kz=1:nz-1
      
      x = dx*(ix-1);
      y = dy*(jy-1);
      z = dz*(kz-1);
      
      Ez_1 = Ez_value( x,         y,         z+dz/2, te    );
      Ez_0 = Ez_value( x,         y,         z+dz/2, te-dt );
      
      Hy_1 = Hy_value( x+dx/2,    y,         z+dz/2, tm    );
      Hy_0 = Hy_value( x+dx/2-dx, y,         z+dz/2, tm    );
      
      Hx_1 = Hx_value( x,         y+dy/2,    z+dz/2, tm    );
      Hx_0 = Hx_value( x,         y+dy/2-dy, z+dz/2, tm    );
      
      Jz(ix, jy, kz) = ( Ez_1 - Ez_0 ) / dt - 1/(dx*dy) * ( dy*(Hy_1 - Hy_0) - dx*(Hx_1 - Hx_0) );
      
      %  Ez(ix, jy, kz) = Ez(ix, jy, kz) + dt/(dx*dy) * (...
      %         + dy*(Hy(ix, jy, kz) - Hy(ix-1, jy, kz))...
      %         - dx*(Hx(ix, jy, kz) - Hx(ix, jy-1, kz))...
      %  ) + dt * Jz(ix, jy, kz);
    end
  end
end

end