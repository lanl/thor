addpath(genpath('../../../../matlab/Maxwell-Mimetic/src/'))
addpath(genpath('../../../../matlab/utils/chebfun/'))
addpath(genpath('../../../../matlab/utils/tt-toolbox/'))
addpath(genpath('../../../../matlab/utils/ttfunc/'))

clc; clear; close all;

divflag = 1;
energyflag = 1;

%% setup run parameters 
starttime = datetime;
% pblm size (points)
nlev = 3;      % 0 -->> 2x2x2
nref = 2^nlev;

% be careful here: (nx,ny,nz) are the number of partitions
nx =  2*nref+1;
ny =  2*nref+1;
nz =  2*nref+1;

X0 = -1; 
X1 = 1;

% mesh size
dx = (X1-X0)./(nx-1);
dy = (X1-X0)./(ny-1);
dz = (X1-X0)./(nz-1);

% check
assert(nx>1);
assert(ny>1);
assert(nz>1);

% CFL factor
CFL = 0.5; 
% CFL = 1/sqrt(3);

% run parameters
dt = CFL * min([dx,dy,dz]);

% start & final time
T0 = 0;
T1 = 20;
nt = floor((T1-T0)/dt);

% add a time cycle to get the final integration time
if T0+nt*dt<T1
  nt = nt+1;
end

% print before starting the solver
fprintf('nlev = %d \n', nlev);
fprintf("number of cells      = %d x %d x %d\n",(nx-1), (ny-1), (nz-1));
fprintf("number of time loops = %d\n",nt);
fprintf("time step  = %e\n",dt);

%% setup magnetic fields and magnetic coeffs
time = T0;
[ Hx, Hy, Hz ] = setup_H( time, dx, dy, dz, nx, ny, nz );

%% setup electric field
time = T0+dt/2;
[ Ex, Ey, Ez ] = setup_E( time, dx, dy, dz, nx, ny, nz, X0 );

%% setup initial current
[ Jx, Jy, Jz ] = setup_Jc( nx, ny, nz );

divH = nan(1,nt); %compute the divergence of H at every time step
%% time loop
time = T0;
ll = 0;
for it=1:nt
  % print time at the end of the integration interval
  if time+dt<T1
    % first nt time steps
    time = time+dt;
    if mod(it,5)==0
      fprintf(repmat('\b',1,ll));
      ll = fprintf("-- %.2f seconds --- it = %d\n",...
        seconds(datetime-starttime),it);
    end
  else
    % last time step to match the final time
    dt    = T1-time;
    time  = time + dt;
    fprintf(repmat('\b',1,ll));
    ll = fprintf("-- %.2f seconds --- it = %d\n",...
        seconds(datetime-starttime),it);
  end

  % update magnetic field components
  Hx = update_Hx( Hx, Ey, Ez, dy, dz, dt, nx, ny, nz );
  Hy = update_Hy( Hy, Ez, Ex, dz, dx, dt, nx, ny, nz );
  Hz = update_Hz( Hz, Ex, Ey, dx, dy, dt, nx, ny, nz );

  % update DBC of H
  [ Hx, Hy, Hz ] = update_DBC_H( Hx, Hy, Hz, dx, dy, dz, nx, ny, nz, time );
 
  % update electric field components
  Ex = update_Ex( Ex, Jx, Hy, Hz, dy, dz, dt, nx, ny, nz );
  Ey = update_Ey( Ey, Jy, Hz, Hx, dz, dx, dt, nx, ny, nz );
  Ez = update_Ez( Ez, Jz, Hx, Hy, dx, dy, dt, nx, ny, nz );
  
  % update DBC of E
  [ Ey, Ez ] = update_DBC_X( Ey, Ez, dx, dy, dz, nx, ny, nz, time+dt/2 );
  [ Ez, Ex ] = update_DBC_Y( Ez, Ex, dx, dy, dz, nx, ny, nz, time+dt/2 );
  [ Ex, Ey ] = update_DBC_Z( Ex, Ey, dx, dy, dz, nx, ny, nz, time+dt/2 );
  
  %% compute the divergence at each time step
  if divflag
    temp = 0;
    temp = temp + (Hx(3:nx-1,2:ny-2,2:nz-2) - Hx(2:nx-2,2:ny-2,2:ny-2))/dx;
    temp = temp + (Hy(2:ny-2,3:ny-1,2:ny-2) - Hy(2:ny-2,2:ny-2,2:ny-2))/dy;
    temp = temp + (Hz(2:ny-2,2:ny-2,3:nz-1) - Hz(2:ny-2,2:ny-2,2:nz-2))/dz;
    divH(1,it) = norm(temp(:));
  end
  %% energy conservation
  if energyflag
    % interpolate the cell center values
    Hxc = 0.5*(Hx(2:end,:,:)+Hx(1:end-1,:,:));
    Hyc = 0.5*(Hy(:,2:end,:)+Hy(:,1:end-1,:));
    Hzc = 0.5*(Hz(:,:,2:end)+Hz(:,:,1:end-1));
    Exc = 0.25*(Ex(:,1:end-1,1:end-1) ...
    + Ex(:,1:end-1,2:end) + Ex(:,2:end,1:end-1) + Ex(:,2:end,2:end));
    Eyc = 0.25*(Ey(1:end-1,:,1:end-1) ...
    + Ey(1:end-1,:,2:end) + Ey(2:end,:,1:end-1) + Ey(2:end,:,2:end));
    Ezc = 0.25*(Ez(1:end-1,1:end-1,:) ...
    + Ez(1:end-1,2:end,:) + Ez(2:end,1:end-1,:) + Ez(2:end,2:end,:));
    % compute H energy
    Energy_H(it) = (dx*dy*dz)*sum((Hxc.^2 + Hyc.^2 + Hzc.^2),'all');
    Energy_E(it) = (dx*dy*dz)*sum((Exc.^2 + Eyc.^2 + Ezc.^2),'all');
  end
end
save('../plot_data/fg_results.mat','T0','T1','nt','dt','divH','Energy_H','Energy_E');