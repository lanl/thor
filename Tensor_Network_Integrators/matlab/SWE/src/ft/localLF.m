function data=localLF(data,d)
    % get size
    [Neq,Nx_total,Ny_total]=size(data.Q);
    
    % Left side
    i = [4:Nx_total-3+(d==1)];
    j = [4:Ny_total-3+(d==2)];
    %
    data.FLq(:,i,j,:) = 0.5*(data.F(data.QL(:,i,j,:),d,data.gacc,data.H) - data.eig(:,i,j,:).*data.QL(:,i,j,:));
    % apply quadrature rule
    FLq            = reshape(data.FLq,Neq*Nx_total*Ny_total,data.quad.n);
    data.FL(:,:,:) = reshape(FLq*data.quad.w,[Neq,Nx_total,Ny_total]);

    % Right side
    i = [4-(d==1):Nx_total-3];
    j = [4-(d==2):Ny_total-3];
    %
    ip1 = i + (1==d);
    jp1 = j + (2==d);
    %
    data.FRq(:,i,j,:) = 0.5*(data.F(data.QR(:,i,j,:),d,data.gacc,data.H) + data.eig(:,ip1,jp1,:).*data.QR(:,i,j,:));
    % apply quadrature rule
    FRq            = reshape(data.FRq,Neq*Nx_total*Ny_total,data.quad.n);
    data.FR(:,:,:) = reshape(FRq*data.quad.w,[Neq,Nx_total,Ny_total]);
    %
end