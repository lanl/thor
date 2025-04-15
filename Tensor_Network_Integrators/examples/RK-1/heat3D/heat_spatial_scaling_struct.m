function heat_spatial_scaling_struct(input,savename,runFG)

addpath(genpath('../../../matlab/RK-1/'))
addpath(genpath('../../../matlab/utils/tt-toolbox/'))
addpath(genpath('../../../matlab/utils/ttfunc/'))

% problem inputs
% tspan - time interval
% nvals - vector of spatial grid sizes
% rtol - relative accuracy tolerance 
% tt_tol - rounding tolerance 
% mname - name of explicit Runge Kutta method
% uexact - anonymous function 
% 
% Presets 
% spatial domain is [-1,1].^3
% set up of right hand side assumes zero boundary conditions so chosen 
% exact solution must have zero boundaries on [-1,1].^3
% Finite difference solution assumed 
% Uniform grid spacing in space 
% Adaptive in time 

%% Unpack input struct 
if contains(path, 'struct2vars')
    struct2vars(input);
else
    tspan = input.tspan; 
    nvals = input.nvals; 
    rtol = input.rtol; 
    atol = input.atol;
    tt_tol = input.tt_tol;
    mname = input.mname; 
    uexact = input.uexact;
    diff = input.diff;
end

%% Parameters  
xl = -1; xr = 1; yl = xl; zl = xl; yr = xr; zr = xr; % spatial domain
prob = generate_test_general(uexact,diff); % generate manufactured problem
R = cell(1,length(nvals));

for j = 1:length(nvals)
    n = nvals(j); 

    %% Spatial Grids
    nx = n; ny = n; nz = n;
    dx = (xr-xl)/(nx-1); 
    dy = (yr-yl)/(ny-1); 
    dz = (zr-zl)/(nz-1);
    fprintf('\n------------------------------------\n');
    fprintf('Grid size = %d \n',n); 

    xp = linspace(xl,xr,nx); 
    yp = linspace(yl,yr,ny);
    zp = linspace(zl,zr,nz);
    
    %fg
    [XP,YP,ZP] = ndgrid(xp,yp,zp);
    [Ix, Iy, Iz] = ndgrid(1:nx, 1:ny, 1:nz);     % index arrays
    boundaryMask = (Ix == 1 | Ix == nx | ...
                    Iy == 1 | Iy == ny | ...
                    Iz == 1 | Iz == nz); % index-based logical mask, true on the boundary 
    % interiorMask = ~boundaryMask;
    
    %tt
    Itt = {ones(1,nx),ones(1,ny),ones(1,nz)};
    temp = Itt; 
    temp{1} = xp; 
    Ctt{1} = cell2core(tt_tensor,temp); 
    temp = Itt; 
    temp{2} = yp; 
    Ctt{2} = cell2core(tt_tensor,temp); 
    temp = Itt; 
    temp{3} = zp; 
    Ctt{3} = cell2core(tt_tensor,temp);
    
    %% Time
    dtinit = 0.01;
    % Stability condition for explicit method
    if dtinit > min([dx, dy, dz])^2 / diff
        disp('Time step size is too large for stability. Reducing dt.');
        dtinit = min([dx, dy, dz])^2 / diff;
        fprintf('dtinit = %g\n',dtinit)
    end
    dtmin = 1e-2*dtinit;  dtmax = 1e2*dtinit; 

    fprintf('dtinit = %g, dtmin = %g, dtmax = %g \n',dtinit,dtmin,dtmax);

    % RK method
    B = butcher(mname);  s = numel(B(1,:))-1;
    
    %% Store problem info
    prob.N = [nx,ny,nz];  
    prob.diffcoeff = diff; 
    prob.domain = [xl,xr,yl,yr,zl,zr];
    prob.deltas = [dx,dy,dz];
    prob.discretepoints = [xp;yp;zp];
    prob.forcingfn_xyz = @(t) prob.forcingfn(t,XP,YP,ZP);
    prob.uexact_xyz = @(t) uexact(t,XP,YP,ZP);
    prob.mask = boundaryMask;
    prob.Ctt = Ctt;
    % disp(prob)
    
    %% Initial conditions and RHS functions 
    uinitial = @(x,y,z) uexact(tspan(1),x,y,z);
    
    % fg
    U0 = uinitial(XP,YP,ZP);
    oderhs = @(t,U) fg_oderhs(t,U,prob);

    % tt
    U0tt = amen_cross_zero(Ctt,@(x) cross_fun_nD(x,uinitial),tt_tol,'verb',0);
    ttrhs = @(t,Utt,ett) tt_oderhs(t,Utt,prob,ett);
    
    % tt_assert_error(U0,U0tt,tt_tol,'U0');
    atolfg = atol*ones(size(U0));
    atoltt = atol*tt_ones(U0tt.n);

    par.B = B;
    par.rtol = rtol;
    par.dtmin = dtmin;
    par.dtmax = dtmax;
    par.dtinit = dtinit;
    par.atol = atolfg;
    

    %% Exact solution 
    Uexact = uexact(tspan(2),XP,YP,ZP);

    %% FG ERK
    if runFG
      fprintf('\nRunning with ERK integrator: %s (order = %i)\n',mname,B(s+1,1))

      tstart = tic;
      [~,Urk,nsteps,fgoutput] = solve_ERK_struct(oderhs,tspan,U0,par);
      FGtime = toc(tstart);  
  
      Urk_end = Urk{end};
      Urk_end(boundaryMask) = Uexact(boundaryMask);

      rkerror = Uexact - Urk_end;
      fgabserror = norm(rkerror,'fro');
      fgrelerror = fgabserror/norm(Uexact,'fro');

      fprintf('FG nsteps = %d\n', nsteps);
      fprintf('FG Elapsed time = %.5f seconds\n', FGtime);
      fprintf('FG Elapsed time per step = %.5f seconds\n', FGtime/nsteps);
      fprintf('FG relative error (interior) = %.5e \n',fgrelerror);
    end

    %% TT ERK
    fprintf('\nRunning with ERK-TT integrator: %s (order = %i)\n',mname,B(s+1,1))
    par.atol = atoltt;
    par.ettinit = tt_tol;
    
    tstart = tic;
    [~,UttRK,nstepstt,ttoutput] = solve_ERKtt_struct(ttrhs,tspan,U0tt,par);
    TTtime = toc(tstart);
    
    UttRK_end = UttRK{end};
    UttRK2full = reshape(full(UttRK_end),nx,ny,nz);
    UttRK2full(boundaryMask) = Uexact(boundaryMask);

    tterror = UttRK2full - Uexact;
    ttabserror = norm(tterror,'fro');
    ttrelerror = ttabserror/norm(Uexact,'fro');
    
    fprintf('TT nsteps = %d\n', nstepstt);
    fprintf('TT Elapsed time = %.5f seconds\n', TTtime);
    fprintf('TT Elapsed time per step = %.5f seconds\n', TTtime/nstepstt);
    fprintf('TT error against exact (interior) = %.5e \n',ttrelerror);

    %% Save data 
    if runFG
      c.timefg  = FGtime;
      c.fgerror_abs = fgabserror;
      c.fgerror_rel = fgrelerror;
      c.fgerkoutput = fgoutput;
      c.fg_steps = nsteps;
    end

    c.input = input; 
    c.timett  = TTtime;
    c.tterror_abs = ttabserror;
    c.tterror_rel = ttrelerror;
    c.dtinit = dtinit;
    c.tterkoutput = ttoutput;
    c.tt_steps = nstepstt;
    R{j} = c;
    save([savename, '.mat'],'R','-v7.3');
    fprintf('output is saved !!!\n')
    
end
