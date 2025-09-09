function data = initMetaData(data, SWEtype, Tfinal, dtOption, CFL, C_eps, Nx, Ny, Lx, Ly, H,...
                             gacc, Coriolis_f, Q0, BC, Recon, Ref, isManuf, isConsErr)
    %
    % set the number of equations for SWEs
    data.tt.Neq = 3;

    % set function defitions
    data.SWEtype   = SWEtype;
    data.isConsErr = isConsErr;
    data.Q0        = Q0;
    data.BC        = BC;
    data.ref       = Ref;
    data.isManuf   = isManuf;
    data.dtOption  = dtOption;
    data.ReconType = Recon;
    %
    if(SWEtype == "linear")
        data.F   = @F_linearSWE;
        data.Eig = @Eig_linearSWE;
    elseif(SWEtype == "nonlinear")
        data.F   = @F_nonlinearSWE;
        data.Eig = @Eig_nonlinearSWE;
    else
        error('Invalid SWEtype: %s. Valid options are "linear" or "nonlinear".',SWEtype)
    end
    %
    % calculate other characteristic scales from available ref conditions
    %
    data.ref.g = (data.ref.U.^2)/data.ref.H;
    data.ref.t = data.ref.L/data.ref.U;
    data.ref.f = 1/data.ref.t;
    %
    % nondimensionalize inputs
    %
    data.Lx     = Lx/data.ref.L;
    data.Ly     = Ly/data.ref.L;
    data.H      = H/data.ref.H;
    data.gacc   = gacc/data.ref.g;
    data.f      = Coriolis_f/data.ref.f;
    data.Tfinal = Tfinal/data.ref.t;  
    %
    % determine the reference scales for conserved variables
    %
    data.Qref = [data.ref.H; data.ref.U; data.ref.U];
    if(data.SWEtype == "nonlinear")
        data.Qref(2:3) = data.ref.H*data.Qref(2:3);
    end
    %
    % In 1-D, the mesh will look like:
    %    |--------|--------|--------| BC |---------|.......|--------| BC |--------|--------|--------|
    %    |   i=1  |   i=2  |   i=3  | BC |   i=4   |       | i=Nx+3 | BC | i=Nx+4 | i=Nx+5 | i=Nx+6 |
    %    |  Ghost |  Ghost |  Ghost | BC |          Interior        | BC |  Ghost |  Ghost |  Ghost |
    %    |--------|--------|--------| BC |---------|.......|--------| BC |--------|--------|--------|
    
    % Calculate cell size in each direction
    dx = data.Lx/Nx;
    dy = data.Ly/Ny;

    % Create ghost cells
    Nx_total = Nx + 6;
    Ny_total = Ny + 6;

    % save grid size info
    data.Nx = Nx; 
    data.Ny = Ny; 
    %
    data.Nx_total = Nx_total; 
    data.Ny_total = Ny_total; 

    % Create cell centers coordinates
    x = [dx:dx:data.Lx]' - dx/2;
    y = [dy:dy:data.Ly]' - dy/2;
    %
    xg = [-2*dx:dx:data.Lx+3*dx] - dx/2;
    yg = [-2*dy:dy:data.Ly+3*dy] - dy/2;
    
    % Calculate cell face areas and volumes
    data.total_vol = data.Lx*data.Ly;
    data.vol       = dx*dy;
    data.area      = [dy dx];
    data.h         = [dx dy];
    data.hmin      = min([dx dy]);
    data.dx        = dx;
    data.dy        = dy;
    data.x         = x;
    data.y         = y;
    data.xg        = xg;
    data.yg        = yg;
    
    % set parameters needed by tt computations
    data.tt.N   = [Nx_total Ny_total];
    data.tt.d   = 2;
    data.tt.xy  = tt_meshgrid_vert(tt_tensor(x),tt_tensor(y));

    %
    data.tt.eps    = 1e-14;                % only as an initial value 
    data.tt.eps_cr = 1e-14;                % only as an initial value 
    data.tt.eps_rk = [1e-14,1e-14,1e-14];  % only as an initial value 
    data.tt.C_eps  = C_eps;

    % initialize explicit solver parameters
    %
    data.rk_time   = 0;                       % time at rk stage 
    data.crk       = [0 1;3/4 1/4;1/3 2/3];   % rk3 coefficients
    data.cdt       = [0; 1; 0.5];             % rk3 dt values
    %
    data.curr_time = 0;                       % current time of the simulations
    data.time_step = 0;                       % time step index
    data.dt        = 0;                       % time step
    data.max_eig   = 0;                       % max eigenvalue for all cells
    data.cfl       = CFL;                     % cfl number
    %
    % High-order method settings
    %
    if(Recon=="Upwind3")
        %
        data.Recon       = @recon3ideal;       % function handle to set the high-order reconstruction method
        data.recon5fun   = [];                 % function handle for reconstruction using cross interpolation (not used for upwind)
        data.recon5GQfun = [];                 % function handle for reconstruction at quadrature points using cross interpolation (not used for upwind)
        data.weno.weps   = [];                 % small number used in the denominator of the WENO smoothness indicator  (not used for upwind)
        data.wenok       = 2;                  % number of reconstruction polynomials in a candidate stencil
        data.hp          = 3;                  % order of accuracy for the spatial discretization
        data.rkmax       = 3;                  % maximum number of Runge-Kutta stages
        Nq               = 2;                  % number of quadrature points
        %
    elseif(Recon=="Upwind5")
        %
        data.Recon       = @recon5ideal;       % function handle to set the high-order reconstruction method
        data.recon5fun   = [];                 % function handle for reconstruction using cross interpolation (not used for upwind)
        data.recon5GQfun = [];                 % function handle for reconstruction at quadrature points using cross interpolation (not used for upwind)
        data.weno.weps   = [];                 % small number used in the denominator of the WENO smoothness indicator  (not used for upwind)
        data.wenok       = 3;                  % number of reconstruction polynomials in a candidate stencil
        data.hp          = 5;                  % order of accuracy for the spatial discretization
        data.rkmax       = 3;                  % maximum number of Runge-Kutta stages
        Nq               = 3;                  % number of quadrature points
        %
    elseif(Recon=="WENO5")
        %
        data.Recon       = @recon5;          % function handle to set the high-order reconstruction method
        data.recon5fun   = @fun_weno5;       % function handle for reconstruction using cross interpolation
        data.recon5GQfun = @fun_weno5GQ;     % function handle for reconstruction at quadrature points using cross interpolation
        data.weno.weps   = min(data.h)^2;    % small number used in the denominator of the WENO smoothness indicator
        data.wenok       = 3;                % number of reconstruction polynomials in a candidate stencil
        data.hp          = 5;                % order of accuracy for the spatial discretization
        data.rkmax       = 3;                % maximum number of Runge-Kutta stages
        Nq               = 3;                % number of quadrature points
        %
    else
        error("%s is an unknown reconstruction method",Recon);
    end

    % create 2D Gauss quadrature rule
    [data.quad.r,data.quad.w] = createQuadRule(Nq);
    % save the number of quad points
    data.quad.n = Nq;

    % create 1-d quadrature points
    xq = zeros(Nx_total*Nq,1);
    yq = zeros(Ny_total*Nq,1);

    % calculate quadratre points in x-dim
    for i=1:Nx_total
        j1 = 1+(i-1)*Nq;
        j2 = i*Nq;
        xq(j1:j2) = (data.quad.r + i - 1 - 3)*data.dx;
    end

    % calculate quadratre points in y-dim
    for i=1:Ny_total
        j1 = 1+(i-1)*Nq;
        j2 = i*Nq;
        yq(j1:j2) = (data.quad.r + i - 1 - 3)*data.dy;
    end
    % calculate the mesh with quadrature points in tt format
    data.tt.q_xy = tt_meshgrid_vert(tt_tensor(xq),tt_tensor(yq));

    %
    % allocate memory for linear weights 
    if(mod(data.quad.n,2)==0)
        data.weno.gamma = zeros(data.wenok,data.quad.n);
    else
        data.weno.gamma = zeros(data.wenok,data.quad.n+1);
    end
    % allocate mem: high-order interpolation coefficients for cell averages
    data.weno.C = cell(1,size(data.weno.gamma,2));
    %
    for i=1:data.quad.n
        data.weno.C{i} = wenokCoeff(data.quad.r(i),data.wenok);
    end
    % calculate linear weno weights
    for i=1:data.quad.n
        if(data.hp==3)
            data.weno.gamma(:,i) = optimalLinearWeno3Weights(data.quad.r(i));
        elseif(data.hp==5)
            data.weno.gamma(:,i) = optimalLinearWeno5Weights(data.quad.r(i));
        else
            data.weno.gamma(:,i)=1;
        end
    end
    % correct interp coefficients and negative weights
    data.weno.imid = -999;
    if(mod(data.quad.n,2)==1)
        %
        j=(data.quad.n-1)/2+1;
        %
        data.weno.imid = j;
        %
        data.weno.C{data.quad.n+1} = data.weno.C{j};
        %
        gamma = data.weno.gamma(:,j);
        %
        data.weno.gamma(:,j) = 0;
        %
        gammatp = 0.5*(3*abs(gamma)+gamma);
        gammatm = 0.5*(3*abs(gamma)-gamma);
        %
        data.weno.sigma_p = sum(gammatp); 
        data.weno.sigma_m = sum(gammatm); 
        %
        data.weno.gamma(:,j)             = gammatp/data.weno.sigma_p;
        data.weno.gamma(:,data.quad.n+1) = gammatm/data.weno.sigma_m;
        %
    end
    %
    % Limited error checks (more to be added later) 
    %
    if(data.isManuf)
        %
        if(~exist("manuf.m","file"))
            error("For problems with a manufactured solution, the relevant 'manuf.m' file is needed in the path");
        end
        %
    end
    %
end