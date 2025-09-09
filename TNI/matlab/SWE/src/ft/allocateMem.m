function data = allocateMem(data)
    % set array sizes
    Neq      = data.Neq;
    Nx       = data.Nx;
    Ny       = data.Ny;
    Nx_total = data.Nx_total;
    Ny_total = data.Ny_total;
    %
    % solution related variables
    %
    data.Q   = zeros(Neq,Nx_total,Ny_total);                   % conserved variables at cell centers
    data.RHS = zeros(Neq,Nx,Ny);                               % right-hand side vector
    data.FL  = zeros(Neq,Nx_total,Ny_total);                   % integrated flux values at the left boundary of the cells
    data.FR  = zeros(Neq,Nx_total,Ny_total);                   % integrated flux values at the right boundary of the cells
    data.QL  = zeros(Neq,Nx_total,Ny_total,data.quad.n);       % reconstructed conserved variable values at the left boundary of the cell on quad points
    data.QR  = zeros(Neq,Nx_total,Ny_total,data.quad.n);       % reconstructed conserved variable values at the right boundary of the cells on quad points
    data.FLq = zeros(Neq,Nx_total,Ny_total,data.quad.n);       % fluxes at the left boundary of the cell on quad points
    data.FRq = zeros(Neq,Nx_total,Ny_total,data.quad.n);       % fluxes at the right boundary of the cell on quad points
    data.eig = zeros(Neq,Nx_total,Ny_total,data.quad.n);       % eigenvalues of flux jacobians
    % allocate source term if needed 
    if (data.isManuf)
        % source term
        data.source = zeros(size(data.RHS));
    end
    %
end