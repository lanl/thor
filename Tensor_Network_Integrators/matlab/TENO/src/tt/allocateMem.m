function data = allocateMem(data)
    %
    data.tt.zero = tt_zeros(data.tt.N, data.tt.d);  
    %
    % solution related variables
    %
    data.tt.Q      = cell(1,data.tt.Neq);                % conserved variables at cell centers
    data.tt.RHS    = cell(1,data.tt.Neq);                % right-hand side vector
    data.tt.FL     = cell(1,data.tt.Neq);                % reconstructed flux values at the left boundary of the cells
    data.tt.FR     = cell(1,data.tt.Neq);                % reconstructed flux values at the right boundary of the cells
    data.tt.invFL  = cell(1,data.tt.Neq);                % split flux values at the left boundary of the cells
    data.tt.invFR  = cell(1,data.tt.Neq);                % split Flux values at the right boundary of the cells
    %
    data.tt.F0L = cell(3,1);                             % reserved for initial guess for flux cross interpolation
    data.tt.F0R = cell(3,1);                             % reserved for initial guess for flux cross interpolation
    %
    for i=1:data.tt.Neq
        data.tt.Q{i}   = data.tt.zero; 
        data.tt.RHS{i} = data.tt.zero; 
        data.tt.FL{i}  = data.tt.zero;
        data.tt.FR{i}  = data.tt.zero;
    end
    % allocate source term if needed
    if (data.isSource)
        data.tt.source = cell(1,data.tt.Neq);
        for i=1:data.tt.Neq
            data.tt.source{i} = data.tt.zero;
        end
    end
    %
    data.tt.eps_rk = cell(data.tt.Neq,1);
    %
end