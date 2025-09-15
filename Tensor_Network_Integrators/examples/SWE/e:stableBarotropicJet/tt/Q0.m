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
    nx = data.Nx_total*data.quad.n;
    ny = data.Ny_total*data.quad.n;
    %
    data.tt.Q0{1} = tt_zeros([nx ny]);
    data.tt.Q0{2} = tt_zeros([nx ny]);
    data.tt.Q0{3} = tt_zeros([nx ny]);
    %
    j = [data.quad.n*3+1:ny-data.quad.n*3];
    %
    core11=data.tt.Q0{1}{1};
    core12=data.tt.Q0{1}{2};
    %
    core21=data.tt.Q0{2}{1};
    core22=data.tt.Q0{2}{2};
    %
    core11(:,:,:) = 1;
    core12(:,j,:) = hb;
    %
    data.tt.Q0{1}{1} = core11;
    data.tt.Q0{1}{2} = core12;
    %
    core21(:,:,:) = 1;
    core22(:,j,:) = hb.*u;
    %
    data.tt.Q0{2}{1} = core21;
    data.tt.Q0{2}{2} = core22;
    %
    Q_quad = data.tt.Q0;
    %
    for i=1:data.tt.Neq
        data.tt.Q0{i}    = applyQuadrature(Q_quad{i},data.quad,data.h);
        data.tt.Q0{i}    = round(data.tt.Q0{i},data.tt.eps_rk(i));
        data.tt.Q0{i}{1} = data.tt.Q0{i}{1}/data.vol;
    end
    %
    %
    data.tt.Q = data.tt.Q0;
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