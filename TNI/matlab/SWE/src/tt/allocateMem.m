function data = allocateMem(data)
    %
    data.tt.zero = tt_zeros(data.tt.N, data.tt.d);  
    %
    % solution related variables
    data.tt.Q   = cell(1,data.tt.Neq);                % conserved variables at cell centers
    data.tt.RHS = cell(1,data.tt.Neq);                % right-hand side vector of forward euler method
    % reconstructed variables
    data.tt.QL  = cell(1,data.tt.Neq);                % reconstructed conserved variable at the left bondary of the cell
    data.tt.QR  = cell(1,data.tt.Neq);                % reconstructed conserved variable at the right boundary of the cell
    data.tt.FL  = cell(1,data.tt.Neq);                % numerical flux at the left boundary of the cell
    data.tt.FR  = cell(1,data.tt.Neq);                % numerical flux at the right boundary of the cell
    %
    for i=1:data.tt.Neq
        data.tt.Q{i}   = data.tt.zero; 
        data.tt.RHS{i} = data.tt.zero; 
        data.tt.QL{i}  = data.tt.zero;     % reconstructed conserved variable at the left bondary of the cell
        data.tt.QR{i}  = data.tt.zero;     % reconstructed conserved variable at the right boundary of the cell
        data.tt.FL{i}  = data.tt.zero;     % numerical flux at the left boundary of the cell
        data.tt.FR{i}  = data.tt.zero;     % numerical flux at the right boundary of the cell
    end  

    % allocate source term if needed
    if (data.isManuf)
        data.tt.source  = cell(1,data.tt.Neq);
        data.tt.sourceQ = cell(1,data.tt.Neq);
        data.tt.Q0      = cell(1,data.tt.Neq);
        for i=1:data.tt.Neq
            data.tt.source{i} = data.tt.zero;
        end
    end
    %
    data.tt.Qq  = cell(1,data.tt.Neq);
    data.Qcore  = cell(1,2);
    %
end