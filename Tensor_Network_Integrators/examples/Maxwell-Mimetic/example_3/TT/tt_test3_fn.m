function [Err, Elapsed_time, tt_numel_E, full_numel_E] = tt_test3_fn(nlev,tol)

starttime = datetime;
% pblm size (points)
% nlev = 10;      % 0 -->> 2x2x2
nref = 2^nlev;
% be careful here: (nx,ny,nz) are the number of partitions
nx =  2*nref+1;
ny =  2*nref+1;
nz =  2*nref+1;

% mesh size
dx = 1./(nx-1);
dy = 1./(ny-1);
dz = 1./(nz-1);

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
fprintf('nlev = %d ---- tol = %.2e \n', nlev, tol);
fprintf("number of cells      = %d x %d x %d\n",(nx-1), (ny-1), (nz-1));
fprintf("number of time loops = %d\n",nt);
fprintf("time step  = %e\n",dt);
fprintf('Expected error = %.5e \n', (7.5)/(4^(nlev-1)));

%% setup magnetic fields and magnetic coeffs
time = T0;

[ tHx, tHy, tHz ] = tt_setup_H(time, dx, dy, dz, nx, ny, nz, tol);

%% setup electric field
time = T0+dt/2;

[ tEx, tEy, tEz ] = tt_setup_E( time, dx, dy, dz, nx, ny, nz, tol);

%% setup initial current
% [ tJx, tJy, tJz ] = tt_setup_Jc( nx, ny, nz );

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
Elapsed_time = seconds(datetime - starttime);
%% compute storage 
[~, tt_numel_E, full_numel_E] = compress_ratio_tt(tEz);
%% error diagnostics: check magnetic field at t = T1
[ terr_H, tnrm_H ] = tt_error_H( tHx, tHy, tHz, dx, dy, dz, nx, ny, nz, time, tol);
% %% error diagnostics: check the electric field at t = T1+dt/2
time = time+dt/2;
[ terr_E, tnrm_E ] = tt_error_E( tEx, tEy, tEz, dx, dy, dz, nx, ny, nz, time, tol);
fprintf('---> Elapsed time = %.2fs \n',seconds(datetime - starttime))
fprintf('E_err + H_err = %.5e \n', terr_E+ terr_H)
fprintf('tt Ez numel = %.2e \n', tt_numel_E)
fprintf('full Ez numel = %.2e \n', full_numel_E)
Err = [terr_H, tnrm_H, terr_E, tnrm_E];
