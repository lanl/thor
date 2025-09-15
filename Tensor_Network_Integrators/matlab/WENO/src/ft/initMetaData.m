function data = initMetaData(data, Tfinal, dtOption, CFL, Neq, Nx, Ny, Nz, Lx, Ly, Lz,...
                             F, Eig, Q0, BC, isSource, gam, isConsErr)

    % set function defitions
    data.isConsErr = isConsErr;
    data.F         = F;
    data.Eig       = Eig;
    data.Q0        = Q0;
    data.BC        = BC;
    data.isSource  = isSource;
    data.dtOption  = dtOption;
    data.gam       = gam;
    %
    % In 1-D, the mesh will look like:
    %    |--------|--------|--------| BC |---------|.......|--------| BC |--------|--------|--------|
    %    |   i=1  |   i=2  |   i=3  | BC |   i=4   |       | i=Nx+3 | BC | i=Nx+4 | i=Nx+5 | i=Nx+6 |
    %    |  Ghost |  Ghost |  Ghost | BC |          Interior        | BC |  Ghost |  Ghost |  Ghost |
    %    |--------|--------|--------| BC |---------|.......|--------| BC |--------|--------|--------|
    
    % Calculate cell size in each direction
    dx = Lx/Nx;
    dy = Ly/Ny;
    dz = Lz/Nz;

    % Create ghost cells
    Nx_total = Nx + 6;
    Ny_total = Ny + 6;
    Nz_total = Nz + 6;

    % save grid size info
    data.Nx = Nx; 
    data.Ny = Ny; 
    data.Nz = Nz;
    %
    data.Nx_total = Nx_total; 
    data.Ny_total = Ny_total; 
    data.Nz_total = Nz_total; 
    %
    data.Neq = Neq;

    % Create cell centers coordinates
    x = [dx:dx:Lx]' - dx/2;
    y = [dy:dy:Ly]' - dy/2;
    z = [dz:dz:Lz]' - dz/2;
    %
    xg = [-2*dx:dx:Lx+3*dx] - dx/2;
    yg = [-2*dy:dy:Ly+3*dy] - dy/2;
    zg = [-2*dz:dz:Lz+3*dz] - dz/2;
    
    % Calculate cell face areas and volumes
    data.vol  = dx*dy*dz;
    data.area = [dy*dz; dx*dz; dx*dy];
    data.h    = [dx dy dz];
    data.dx   = dx;
    data.dy   = dy;
    data.dz   = dz;
    data.x    = x;
    data.y    = y;
    data.z    = z;
    data.xg   = xg;
    data.yg   = yg;
    data.zg   = zg;
    %
    data.X = zeros(length(xg),length(yg),length(zg));
    data.Y = data.X;
    data.Z = data.X;
    %
    for k=1:length(zg)
        for j=1:length(yg)
            for i=1:length(xg)
                 data.X(i,j,k) = xg(i);
                 data.Y(i,j,k) = yg(j);
                 data.Z(i,j,k) = zg(k);
            end
        end
    end
    %
    % initialize explicit solver parameters
    %
    data.rk_time   = 0;                       % time at rk stage 
    data.crk       = [0 1;3/4 1/4;1/3 2/3];   % rk3 coefficients
    data.cdt       = [0; 1; 0.5];             % rk3 dt values
    data.rkmax     = 3;
    data.hp        = 5;
    %
    data.curr_time = 0;       % current time of the simulations
    data.time_step = 0;       % time step index
    data.dt        = 0;       % time step
    data.Tfinal    = Tfinal;  % final time
    data.max_eig   = 0;       % max eigenvalue for all cells
    data.cfl       = CFL;     % cfl number
    %
    % save the number of equations
    %
    data.weno.weps   = min(data.h)^2;
    %
    % Limited error checks (more to be added later) 
    %
    if(data.isSource)
        %
        if(~exist("sourceTerm.m","file"))
            error("For problems with a source term, the relevant 'sourceTerm.m' file is needed in the path");
        end
        %
    end
    %
end