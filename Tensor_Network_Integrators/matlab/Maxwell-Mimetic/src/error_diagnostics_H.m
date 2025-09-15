function [ err_H, nrm_H ] = error_diagnostics_H( Hx, Hy, Hz, dx, dy, dz, nx, ny, nz, time )

% set time
t = time;

% Hx
err_Hx=0;
nrm_Hx=0;
for ix=1:nx
  for jy=1:ny-1
    for kz=1:nz-1
      x = dx*(ix-1);
      y = dy*(jy-1)+dy/2;
      z = dz*(kz-1)+dz/2;
      err_Hx = err_Hx + (Hx(ix, jy, kz) - Hx_value( x, y, z, t ))^2; 
      nrm_Hx = nrm_Hx + Hx_value( x, y, z, t )^2;
      %fprintf("(x,y,z,t) = (%f %f %f %f) -->> Hx_num = %14.7e  Hx_exa = %14.7e\n",x,y,z,t,Hx(ix,jy,kz),Hx_value(x,y,z,t));
    end
  end
end
err_Hx = sqrt( err_Hx / ( nx*(ny-1)*(nz-1) ) );
nrm_Hx = sqrt( nrm_Hx / ( nx*(ny-1)*(nz-1) ) );

% Hy
err_Hy=0;
nrm_Hy=0;
for ix=1:nx-1
  for jy=1:ny
    for kz=1:nz-1
      x = dx*(ix-1)+dx/2;
      y = dy*(jy-1);
      z = dz*(kz-1)+dz/2;
      err_Hy = err_Hy + (Hy(ix, jy, kz) - Hy_value( x, y, z, t ))^2; 
      nrm_Hy = nrm_Hy + Hy_value( x, y, z, t )^2;
      % fprintf("(x,y,z,t) = (%f %f %f %f) -->> Hy_num = %14.7e  Hy_exa = %14.7e\n",x,y,z,t,Hy(ix,jy,kz),Hy_value(x,y,z,t));
    end
  end
end
err_Hy = sqrt( err_Hy / ( (nx-1)*ny*(nz-1) ) );
nrm_Hy = sqrt( nrm_Hy / ( (nx-1)*ny*(nz-1) ) );

% Hz
err_Hz=0;
nrm_Hz=0;
for ix=1:nx-1
  for jy=1:ny-1
    for kz=1:nz
      x = dx*(ix-1)+dx/2;
      y = dy*(jy-1)+dy/2;
      z = dz*(kz-1);
      err_Hz = err_Hz + (Hz(ix, jy, kz) - Hz_value( x, y, z, t ))^2; 
      nrm_Hz = nrm_Hz + Hz_value( x, y, z, t )^2; 
      % fprintf("(x,y,z,t) = (%f %f %f %f) -->> Hz_num = %14.7e  Hz_exa = %14.7e\n",x,y,z,t,Hz(ix,jy,kz),Hz_value(x,y,z,t));
    end
  end
end
err_Hz = sqrt( err_Hz / ( (nx-1)*(ny-1)*nz ) );
nrm_Hz = sqrt( nrm_Hz / ( (nx-1)*(ny-1)*nz ) );

% final
err_H = err_Hx + err_Hy + err_Hz;
nrm_H = nrm_Hx + nrm_Hy + nrm_Hz;

% fprintf("pause here!\n");
% fprintf("\n");

end