addpath(genpath('../../../../matlab/Maxwell-Mimetic/src/'))
addpath(genpath('../../../../matlab/utils/chebfun/'))
addpath(genpath('../../../../matlab/utils/tt-toolbox/'))
addpath(genpath('../../../../matlab/utils/ttfunc/'))


close all; clear; clc;

% this code computes the divH for the tt version
%%
ts = datetime;
% pblm size (points)
nlev = 3;      % 0 -->> 2x2x2
nref = 2^nlev;
tol = 1e-12;

% be careful here: (nx,ny,nz) are the number of partitions
nx =  2*nref+1;
ny =  2*nref+1;
nz =  2*nref+1;

X0 = -1;
X1 = 1;

% mesh size --- the domain now is [-1,1]
dx = (X1-X0)./(nx-1);
dy = (X1-X0)./(ny-1);
dz = (X1-X0)./(nz-1);

% create grid

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
fprintf('grid size %d x %d x %d \n', nx-1, ny-1, nz-1);
fprintf('tol = %.2e \n', tol)
% fprintf("number of cells      = %d\n",(nx-1)*(ny-1)*(nz-1));
fprintf("number of time loops = %d\n",nt);
fprintf("time step size  = %.5e\n",dt);
fprintf('End time = %.2fs \n', T1)
%% setup magnetic fields and magnetic coeffs
time = T0;
[ tHx, tHy, tHz ] = tt_setup_H(dx, dy, dz, nx, ny, nz, X0);
%% setup electric field
time = T0+dt/2;
[ tEx, tEy, tEz ] = tt_setup_E(dx, dy, dz, nx, ny, nz, X0);

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
      ll=fprintf("%.2f s --- it = %3d \n",seconds(datetime-ts), it);
    end
  else
    % last time step to match the final time
    dt    = T1-time;
    time  = time + dt;
    fprintf(repmat('\b',1,ll));
    %     ll=fprintf("it = %3d -->> time = %.3fs (final time step = %.3f)\n",it,time,dt );
    ll=fprintf("%.2f s --- it = %3d \n",seconds(datetime-ts), it);
  end

  %% update magnetic field components
  rHz = max(tHz.r);
  tHx = tt_update_Hx(tHx,tEy,tEz,dy,dz,dt,tol);
  tHy = tt_update_Hy(tHy,tEz,tEx,dz,dx,dt,tol);
  tHz = tt_update_Hz(tHz,tEx,tEy,dx,dy,dt,tol);
  %% update DBC of H
  [ tHx, tHy, tHz ] = tt_update_HBC_H(tHx, tHy, tHz, dx, dy, dz, nx, ny, nz, time, tol);

  %% update electric field components

  tEx = tt_update_Ex( tEx, tHy, tHz, dy, dz, dt, nx, ny, nz, tol);
  tEy = tt_update_Ey( tEy, tHz, tHx, dz, dx, dt, nx, ny, nz, tol);
  tEz = tt_update_Ez( tEz, tHx, tHy, dx, dy, dt, nx, ny, nz, tol);
  if mod(it,5)==10
    ll1=fprintf('tEx max rank = %d \n',max(tEx.r));
    ll2=fprintf('tEy max rank = %d \n',max(tEy.r));
    ll3=fprintf('tEz max rank = %d \n',max(tEz.r));
    ll = ll + ll1 + ll2 + ll3;
  end
  %% update DBC of E

  [ tEy, tEz ] = tt_update_HBC_X( tEy, tEz, dx, dy, dz, nx, ny, nz, time+dt/2, tol);
  [ tEz, tEx ] = tt_update_HBC_Y( tEz, tEx, dx, dy, dz, nx, ny, nz, time+dt/2, tol );
  [ tEx, tEy ] = tt_update_HBC_Z( tEx, tEy, dx, dy, dz, nx, ny, nz, time+dt/2, tol );

  %% compute divergence
  if 1
    temp = 0;
    temp = temp + (tHx(3:nx-1,2:nx-2,2:nx-2) - tHx(2:nx-2,2:ny-2,2:nz-2))/dx;
    temp = temp + (tHy(2:nx-2,3:ny-1,2:nx-2) - tHy(2:nx-2,2:ny-2,2:nz-2))/dy;
    temp = temp + (tHz(2:nx-2,2:ny-2,3:nz-1) - tHz(2:nx-2,2:ny-2,2:ny-2))/dz;
    divH(1,it) = norm(temp(:));
  end
  %% compute the energy
  if 1
    %interpolation
    tHxc = tt_intp_Hc_fn(tHx,1, tol);
    tHyc = tt_intp_Hc_fn(tHy,2, tol);
    tHzc = tt_intp_Hc_fn(tHz,3, tol);

    tExc = tt_intp_Ec_fn(tEx,[2,3], tol);
    tEyc = tt_intp_Ec_fn(tEy,[1,3], tol);
    tEzc = tt_intp_Ec_fn(tEz,[1,2], tol);

    %compute energy
    Energy_H(it) = (dx*dy*dz)*sum((tHxc.^2 + tHyc.^2 + tHzc.^2));
    Energy_E(it) = (dx*dy*dz)*sum((tExc.^2 + tEyc.^2 + tEzc.^2));

  end
end
Elapsed_time = seconds(datetime-ts);
fprintf('Elapsed time = %.2fs \n', Elapsed_time);

save('../plot_data/tt_results.mat','T0','T1','nt','dt','divH','Energy_H','Energy_E');






