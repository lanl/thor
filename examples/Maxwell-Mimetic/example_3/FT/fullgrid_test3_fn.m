function [Err, Elapsed_time] =  fullgrid_test3_fn(nlev)
% what to output
% Error and norm of E and H
% Elapsed time

starttime = datetime;
% nlev = 2;      % 0 -->> 2x2x2
nref = 2^nlev;

% be careful here: (nx,ny,nz) are the number of partitions
nx =  2*nref+1;
ny =  2*nref+1;
nz =  2*nref+1;

% mesh size
dx = 1./(nx-1);
dy = 1./(ny-1);
dz = 1./(nz-1);

% check
assert(nx>1);
assert(ny>1);
assert(nz>1);

% CFL factor
CFL = 0.5; 

% run parameters
dt = CFL * min([dx,dy,dz]);

% start & final time
T0 = 0.2;
T1 = 1.2;
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
fprintf('Expected error = %.5e \n', 7.58843/(4^(nlev-1)));

%% setup magnetic fields and magnetic coeffs
time = T0;
[ Hx, Hy, Hz ] = setup_H( time, dx, dy, dz, nx, ny, nz );

%% setup electric field
time = T0+dt/2;
[ Ex, Ey, Ez ] = setup_E( time, dx, dy, dz, nx, ny, nz );

%% setup initial current
[ Jx, Jy, Jz ] = setup_Jc( nx, ny, nz );

%% time loop
time = T0;
ll=0;
for it=1:nt
  % print time at the end of the integration interval
  if time+dt<T1
    % first nt time steps
    time = time+dt;
    if mod(it,5)==0
      fprintf(repmat('\b',1,ll));
      ll = fprintf("-- %.2f seconds --- it = %d",...
        seconds(datetime-starttime),it);
    end
  else
    % last time step to match the final time
    dt    = T1-time;
    time  = time + dt;
    fprintf(repmat('\b',1,ll));
    ll = fprintf("-- %.2f seconds --- it = %d",...
        seconds(datetime-starttime),it);
  end

  % update magnetic field components
  Hx = update_Hx( Hx, Ey, Ez, dy, dz, dt, nx, ny, nz );
  Hy = update_Hy( Hy, Ez, Ex, dz, dx, dt, nx, ny, nz );
  Hz = update_Hz( Hz, Ex, Ey, dx, dy, dt, nx, ny, nz );

  % update DBC of H
  [ Hx, Hy, Hz ] = update_DBC_H( Hx, Hy, Hz, dx, dy, dz, nx, ny, nz, time );
  
  % setup initial current (running time is set inside)
  [ Jx, Jy, Jz ] = update_Jc( Jx, Jy, Jz, dx, dy, dz, dt, nx, ny, nz, time );
  
  % update electric field components
  Ex = update_Ex( Ex, Jx, Hy, Hz, dy, dz, dt, nx, ny, nz );
  Ey = update_Ey( Ey, Jy, Hz, Hx, dz, dx, dt, nx, ny, nz );
  Ez = update_Ez( Ez, Jz, Hx, Hy, dx, dy, dt, nx, ny, nz );
  
  % update DBC of E
  [ Ey, Ez ] = update_DBC_X( Ey, Ez, dx, dy, dz, nx, ny, nz, time+dt/2 );
  [ Ez, Ex ] = update_DBC_Y( Ez, Ex, dx, dy, dz, nx, ny, nz, time+dt/2 );
  [ Ex, Ey ] = update_DBC_Z( Ex, Ey, dx, dy, dz, nx, ny, nz, time+dt/2 );
  
end
Elapsed_time = seconds(datetime - starttime);
%% error diagnostics: check magnetic field at t = T1
[err_H, nrm_H ] = error_diagnostics_H( Hx, Hy, Hz, dx, dy, dz, nx, ny, nz, time );

%% error diagnostics: check the electric field at t = T1+dt/2 
time = time+dt/2;
[err_E, nrm_E ] = error_diagnostics_E( Ex, Ey, Ez, dx, dy, dz, nx, ny, nz, time );
fprintf('---> Elapsed time = %.2fs \n', Elapsed_time)
fprintf('E_err + H_err = %.5e \n', err_E+ err_H)
Err = [err_H, nrm_H, err_E, nrm_E];
