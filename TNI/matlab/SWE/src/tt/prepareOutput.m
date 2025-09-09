function data = prepareOutput(data)
    %
    data.Q = remove_ghost(data.tt.Q,data.tt.Neq);
    %
    for i=1:data.tt.Neq
        data.Q{i} = data.Q{i}*data.Qref(i);
    end
    %
    data.x = data.x*data.ref.L;
    data.y = data.y*data.ref.L;
    %
    data.X = zeros(data.Nx,data.Ny);
    data.Y = zeros(data.Nx,data.Ny);
    %
    for i=1:data.Nx
        for j=1:data.Ny
            data.X(i,j) = data.x(i);
            data.Y(i,j) = data.y(j);
        end
    end
    %
end