function [ err_E, nrm_E ] = error_diagnostics_E( Ex, Ey, Ez, dx, dy, dz, nx, ny, nz, time )

% set time
t = time;

% Ex
err_Ex=0;
nrm_Ex=0;
for ix=1:nx-1
  for jy=1:ny
    for kz=1:nz
      x = dx*(ix-1)+dx/2;
      y = dy*(jy-1);
      z = dz*(kz-1);
      err_Ex = err_Ex + (Ex(ix, jy, kz)-Ex_value( x, y, z, t ))^2;
      nrm_Ex = nrm_Ex + Ex_value( x, y, z, t )^2;
      % fprintf("(x,y,z,t) = (%f %f %f %f) -->> Ex_num = %14.7e  Ex_exa = %14.7e\n",x,y,z,t,Ex(ix,jy,kz),Ex_value(x,y,z,t));
    end
  end
end
err_Ex = sqrt( err_Ex / ( (nx-1)*ny*nz ) );
nrm_Ex = sqrt( nrm_Ex / ( (nx-1)*ny*nz ) );

% Ey
err_Ey=0;
nrm_Ey=0;
for ix=1:nx
  for jy=1:ny-1
    for kz=1:nz
      x = dx*(ix-1);
      y = dy*(jy-1)+dy/2;
      z = dz*(kz-1);
      err_Ey = err_Ey + (Ey(ix, jy, kz) - Ey_value( x, y, z, t ))^2;
      nrm_Ey = nrm_Ey + Ey_value( x, y, z, t )^2;
      % fprintf("(x,y,z,t) = (%f %f %f %f) -->> Ey_num = %14.7e  Ey_exa = %14.7e\n",x,y,z,t,Ey(ix,jy,kz),Ey_value(x,y,z,t));
    end
  end
end
err_Ey = sqrt( err_Ey / ( nx*(ny-1)*nz ) );
nrm_Ey = sqrt( nrm_Ey / ( nx*(ny-1)*nz ) );

% Ez
err_Ez=0;
nrm_Ez=0;
for ix=1:nx
  for jy=1:ny
    for kz=1:nz-1
      x = dx*(ix-1);
      y = dy*(jy-1);
      z = dz*(kz-1)+dz/2;
      err_Ez = err_Ez + (Ez(ix, jy, kz) - Ez_value( x, y, z, t ))^2;
      nrm_Ez = nrm_Ez + Ez_value( x, y, z, t )^2;
      % fprintf("(x,y,z,t) = (%f %f %f %f) -->> Ez_num = %14.7e  Ez_exa = %14.7e\n",x,y,z,t,Ez(ix,jy,kz),Ez_value(x,y,z,t));
    end
  end
end
err_Ez = sqrt( err_Ez / ( nx*ny*(nz-1) ) );
nrm_Ez = sqrt( nrm_Ez / ( nx*ny*(nz-1) ) );

% final
err_E = err_Ex + err_Ey + err_Ez;
nrm_E = nrm_Ex + nrm_Ey + nrm_Ez;

% fprintf("pause here!\n");
% fprintf("\n");

end