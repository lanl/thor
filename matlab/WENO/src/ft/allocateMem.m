function data = allocateMem(data)
    % set array sizes
    Neq      = data.Neq;
    Nx       = data.Nx;
    Ny       = data.Ny;
    Nz       = data.Nz;
    Nx_total = data.Nx_total;
    Ny_total = data.Ny_total;
    Nz_total = data.Nz_total;
    %
    % solution related variables
    %
    data.Q     = zeros(Neq,Nx_total,Ny_total,Nz_total);        % conserved variables at cell centers
    data.RHS   = zeros(Neq,Nx,Ny,Nz);                          % right-hand side vector
    data.eig   = zeros(Neq,Nx_total,Ny_total,Nz_total);        % eigenvalues of flux jacobians
    data.FL    = zeros(Neq,Nx_total,Ny_total,Nz_total);        % reconstructed flux values at the left boundary of the cells
    data.FR    = zeros(Neq,Nx_total,Ny_total,Nz_total);        % reconstructed flux values at the right boundary of the cells
    data.invFL = zeros(Neq,Nx_total,Ny_total,Nz_total);        % split flux values at the left boundary of the cells
    data.invFR = zeros(Neq,Nx_total,Ny_total,Nz_total);        % split Flux values at the right boundary of the cells
    % allocate source term if needed  
    if(data.isSource)
        data.source = zeros(size(data.RHS));
    end
    %
end
