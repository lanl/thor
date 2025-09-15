function data = prepareOutput(data)
    % remove ghost cell values from the tt
    data.tt.Q = remove_ghost(data.tt.Q,data.tt.Neq);

    % convert solution from tt to ft
    data.ft.Q = zeros(data.tt.Neq,data.Nx,data.Ny,data.Nz);
    for i=1:data.tt.Neq
        data.ft.Q(i,:,:,:) = full(data.tt.Q{i},data.tt.Q{i}.n');
    end
    
    % set cell center coordinates for post-processing
    data.X = zeros(data.Nx,data.Ny,data.Nz);
    data.Y = zeros(data.Nx,data.Ny,data.Nz);
    data.Z = zeros(data.Nx,data.Ny,data.Nz);
    %
    for k=1:data.Nz
       for j=1:data.Ny
           for i=1:data.Nx
               data.X(i,j,k) = data.x(i);
               data.Y(i,j,k) = data.y(j);
               data.Z(i,j,k) = data.z(k);
           end
       end
    end
    %
end
