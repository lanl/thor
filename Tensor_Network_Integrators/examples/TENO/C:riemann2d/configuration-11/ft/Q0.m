function data=Q0(data)
    %
    gam = data.gam;
    %
    rho1 = 1;
    u1   = 0.1;
    v1   = 0;
    p1   = 1;
    %
    rho2 = 0.5313;
    u2   = 0.8276;
    v2   = 0;
    p2   = 0.4;
    %
    rho3 = 0.8;
    u3   = 0.1;
    v3   = 0;
    p3   = 0.4;
    %
    rho4 = 0.5313;
    u4   = 0.1;
    v4   = 0.7276;
    p4   = 0.4;
    %
    idx1 = data.X >= 0.5 & data.Y >= 0.5;
    idx2 = data.X  < 0.5 & data.Y >= 0.5;
    idx3 = data.X  < 0.5 & data.Y  < 0.5;
    idx4 = data.X >= 0.5 & data.Y  < 0.5;
    %
    rho = zeros(size(data.X));
    u   = zeros(size(data.X));
    v   = zeros(size(data.X));
    p   = zeros(size(data.X));
    %
    rho(idx1) = rho1;
    u(idx1)   = u1;
    v(idx1)   = v1;
    p(idx1)   = p1;
    %
    rho(idx2) = rho2;
    u(idx2)   = u2;
    v(idx2)   = v2;
    p(idx2)   = p2;
    %
    rho(idx3) = rho3;
    u(idx3)   = u3;
    v(idx3)   = v3;
    p(idx3)   = p3;
    %
    rho(idx4) = rho4;
    u(idx4)   = u4;
    v(idx4)   = v4;
    p(idx4)   = p4;
    %
    data.Q(1,:,:,:) = rho;
    data.Q(2,:,:,:) = rho.*u;
    data.Q(3,:,:,:) = rho.*v;
    data.Q(4,:,:,:) = 0;
    data.Q(5,:,:,:) = p/(gam-1) + 0.5*rho.*(u.^2 + v.^2);
    %
end