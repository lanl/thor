function [errH, errE, T] = tt_3D_manufactured_fn(test,nlev,tol,T0,T1)
% 0 = convergence run; 1,2,3,4,5... = patch tests
% (1-4 test exact solution, error must be zero)
% (5   test BCs for quadratic convergence)
tic()
global patch_test;
patch_test = test;

% pblm size (points)
% nlev = 2;      % 0 -->> 2x2x2
nref = 2^nlev;

% tol = 1e-15;

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
CFL = 0.2;
% CFL = 1/sqrt(3);

% run parameters
dt = CFL * min([dx,dy,dz]);

% start & final time
% T0 = 0.25;
% T1 = 1.25;
nt = floor((T1-T0)/dt);

% add a time cycle to get the final integration time
if T0+nt*dt<T1
  nt = nt+1;
end

% print before starting the solver
fprintf("number of cells      = %d\n",(nx-1)*(ny-1)*(nz-1));
fprintf("number of time loops = %d\n",nt);
fprintf("initial   time step  = %e\n",dt);

%% setup magnetic fields and magnetic coeffs
time = T0;

[ tHx, tHy, tHz ] = tt_setup_H(time, dx, dy, dz, nx, ny, nz );

%% setup electric field
time = T0+dt/2;

[ tEx, tEy, tEz ] = tt_setup_E( time, dx, dy, dz, nx, ny, nz);

%% setup initial current
% [ tJx, tJy, tJz ] = tt_setup_Jc( nx, ny, nz );

%% time loop
time = T0;
for it=1:nt
  
  % print time at the end of the integration interval
  if time+dt<T1
    % first nt time steps
    time = time+dt;
    if mod(it,5)==0 || true
      fprintf("it = %3d -->> time = %14.7e\n",it,time);
    end
  else
    % last time step to match the final time
    dt    = T1-time;
    time  = time + dt;
    fprintf("it = %3d -->> time = %14.7e (final time step = %e)\n",it,time,dt );
  end
  
  %% update magnetic field components

  tHx = tt_update_Hx(tHx,tEy,tEz,dy,dz,dt,tol);
  tHy = tt_update_Hy(tHy,tEz,tEx,dz,dx,dt,tol);
  tHz = tt_update_Hz(tHz,tEx,tEy,dx,dy,dt,tol);
  
  %% update DBC of H
  [ tHx, tHy, tHz ] = tt_update_DBC_H(tHx, tHy, tHz, dx, dy, dz, nx, ny, nz, time, tol);
  
  %% update Jc
  [ tJx, tJy, tJz ] = tt_update_Jc(dx, dy, dz, dt, nx, ny, nz, time,tol);

  %% update electric field components 
  
  tEx = tt_update_Ex( tEx, tJx, tHy, tHz, dy, dz, dt, nx, ny, nz, tol);
  tEy = tt_update_Ey( tEy, tJy, tHz, tHx, dz, dx, dt, nx, ny, nz, tol);
  tEz = tt_update_Ez( tEz, tJz, tHx, tHy, dx, dy, dt, nx, ny, nz, tol);
  
  %% update DBC of E
  
  [ tEy, tEz ] = tt_update_DBC_X( tEy, tEz, dx, dy, dz, nx, ny, nz, time+dt/2, tol);
  [ tEz, tEx ] = tt_update_DBC_Y( tEz, tEx, dx, dy, dz, nx, ny, nz, time+dt/2,tol );
  [ tEx, tEy ] = tt_update_DBC_Z( tEx, tEy, dx, dy, dz, nx, ny, nz, time+dt/2,tol );
end
%assert( abs(time-T1)<1E-12 );
fprintf("\n");
T = toc;
%% error diagnostics: check magnetic field at t = T1
[ terr_H, tnrm_H ] = error_diagnostics_H( tHx, tHy, tHz, dx, dy, dz, nx, ny, nz, time );
fprintf("error diagnostics -->> final time = %14.7e\n",time);
fprintf("terr_H = %14.7e  tnrm_H = %14.7e  terr_H/tnrm_H = %14.7e\n\n",terr_H,tnrm_H,terr_H/tnrm_H);
errH = terr_H/tnrm_H;
%% error diagnostics: check the electric field at t = T1+dt/2
time = time+dt/2;
[ terr_E, tnrm_E ] = error_diagnostics_E( tEx, tEy, tEz, dx, dy, dz, nx, ny, nz, time );
fprintf("error diagnostics -->> final time = %14.7e\n",time);
fprintf("terr_E = %14.7e  tnrm_E = %14.7e  terr_E/tnrm_E = %14.7e\n\n",terr_E,tnrm_E,terr_E/tnrm_E);
errE = terr_E/tnrm_E;
