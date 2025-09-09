function data = BC(data,t)
    %
    Neq = data.Neq;
    Nx  = data.Nx_total;
    Ny  = data.Ny_total;
    %
    dx  = data.dx;
    dy  = data.dy;
    %
    wij = data.quad.wij;
    rx  = data.quad.rx;
    ry  = data.quad.ry;
    %
    I = zeros(Neq,6,Ny);
    %
    ii = [1:3 Nx-2:Nx];
    %
    for j = 1:Ny
        for cnt=1:length(ii)
            i = ii(cnt);
            %
            qx = (rx + i - 1 - 3)*dx;
            qy = (ry + j - 1 - 3)*dy;
            %
            Int = solExact(qx,qy,t,data);
            %
            for eq=1:Neq
                I(eq,cnt,j) = sum(wij(:)'.*Int(eq,:),"all");
            end
        end
    end
    %
    % non-periodic in x (exact)
    %
    data.Q(:,1:3,:)      = I(:,1:3,:);      
    data.Q(:,end-2:end,:)= I(:,4:6,:);
    % periodic in y
    data.Q(:,:,1:3)      = data.Q(:,:,end-5:end-3);
    data.Q(:,:,end-2:end)= data.Q(:,:,4:6);
    %
end
  