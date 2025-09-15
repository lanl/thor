function data = Q0(data)
    %
    uint = calc_uint(data);
    %
    hb = calc_hb(data,uint); 
    %
    u = calc_U(data); 
    %
    % setup index arrays
    %
    i = [4:data.Nx_total-3];
    j = [4:data.Ny_total-3];
    %
    Nx = data.Nx_total;
    Ny = data.Ny_total;
    %
    data.Q(:,i,j) = integrate_q(data,u,hb);
    %
    data.Q(:,1:3,:)     = data.Q(:,Nx-5:Nx-3,:);
    data.Q(:,Nx-2:Nx,:) = data.Q(:,4:6,:);
    %
    data.Q(:,:,1:3)     = data.Q(:,:,[4 4 4]);
    data.Q(:,:,Ny-2:Ny) = data.Q(:,:,[Ny-3 Ny-3 Ny-3]);
    %
    data.Qe = data.Q(:,4:end-3,4:end-3);
    %
end
%
function uint = calc_uint(data)
    %
    phi0 = -pi/7;
    phi1 =  pi/7;
    %
    Ny  = data.Ny;
    %
    r = data.quad.r;
    w = data.quad.w;
    n = data.quad.n;
    %
    % integrate u(y)dy on fine grid cells using gauss quadratures
    %
    uint  = zeros(Ny*(n+1),1);
    %
    dr = [0;r;1]; dr = dr(2:end)-dr(1:end-1);
    %
    idx=1;
    %
    for j=1:Ny
        %
        yl = (j-1)*data.dy; % left side of the coarse cell
        %
        ylf = yl; % left side of the fine cell
        %
        for i1=1:n+1
            %
            dy = dr(i1)*data.dy;
            %
            for i2=1:n
                %
                y = ylf + r(i2)*dy; % quadrature point on the fine cell
                %
                phi = 2*pi*y/data.Lx - pi/2;
                %
                term1 = 4/(phi1-phi0)^2;
                %
                if(phi > phi0 && phi < phi1)
                    %
                    term2 = 1/(phi-phi1)/(phi-phi0);
                    %
                    Uq = 80/data.ref.U*exp(term1 + term2);
                    %
                else
                    %
                    Uq = 0;
                    %
                end
                %
                uint(idx) = uint(idx) + w(i2)*Uq*dy;
                %
            end
            %
            ylf = ylf + dy;
            %
            idx = idx+1;
            %
        end
        %
    end
    %
end
%
function hb = calc_hb(data,uint)
    %
    Ny  = data.Ny;
    %
    r = data.quad.r;
    w = data.quad.w;
    n = data.quad.n;
    %
    q = zeros(Ny*n,1);
    %
    idx1 = 1;
    idx2 = 1;
    %
    int = 0;
    %
    for j=1:Ny
        %
        for qj=1:n
            %
            int = int + uint(idx2);
            %
            q(idx1) = int;
            %
            idx1 = idx1 + 1;
            idx2 = idx2 + 1;
            %
        end
        % add the remaining part in the coarse cell after the last actual quad point
        int  = int + uint(idx2);
        idx2 = idx2 + 1;
        %
    end
    %
    temp = data.f/data.gacc*q;
    %
    h0 = 0;
    %
    idx1=1;
    for j=1:Ny
        %
        for qj=1:n
            %
            h0   = h0 + w(qj)*temp(idx1)*data.dy;
            idx1 = idx1 + 1;
            %
        end
        %
    end
    %
    h0 = 10000/data.ref.H + h0/data.Ly;
    %
    hb = h0 - temp;
    %
    h0 = 0;
    %
    idx1=1;
    for j=1:Ny
        %
        for qj=1:n
            %
            h0   = h0 + w(qj)*hb(idx1)*data.dy;
            idx1 = idx1 + 1;
            %
        end
        %
    end
    %
end
%
function u = calc_U(data)
    %
    phi0 = -pi/7;
    phi1 =  pi/7;
    %
    Ny  = data.Ny;
    %
    r = data.quad.r;
    w = data.quad.w;
    n = data.quad.n;
    %
    u  = zeros(Ny*n,1);
    %
    idx=1;
    %
    for j=1:Ny
        %
        yl = (j-1)*data.dy; % 
        %
        for i1=1:n
            %
            y = yl + r(i1)*data.dy;
            %
            phi = 2*pi*y/data.Lx  - pi/2;
            %
            term1 = 4/(phi1-phi0)^2;
            %
            if(phi > phi0 && phi < phi1)
                %
                term2 = 1/(phi-phi1)/(phi-phi0);
                %
                u(idx) = 80/data.ref.U*exp(term1 + term2);
                %
            else
                %
                u(idx) = 0;
                %
            end
            %
            idx = idx+1;
            %
        end
        %
    end
    %
end
%
function Q = integrate_q(data,u,hb)
    %
    Nx  = data.Nx;
    Ny  = data.Ny;
    %
    Q = zeros([3,Nx,Ny]);
    %
    r = data.quad.r;
    w = data.quad.w;
    n = data.quad.n;
    %
    idx=1;
    %
    for j=1:Ny
        %
        h  = 0;
        hU = 0;
        for i1=1:n
            %
            h  = h  + w(i1)*hb(idx);
            hU = hU + w(i1)*hb(idx)*u(idx);
            %
            idx = idx+1;
            %
        end
        %
        Q(1,:,j) = h;
        Q(2,:,j) = hU;
        %
    end
    %
end