function data=Q0(data)
    %
    gam = data.gam;
    t   = 0;
    %
    th   = pi/3;
    %
    rhoL =  8;
    uL   =  8.25*sin(th);
    vL   = -8.25*cos(th);
    pL   =  116.5; 
    %
    rhoR = 1.4;
    uR   = 0;
    vR   = 0;
    pR   = 1; 
    %
    xs = 10*t/sin(th) + data.Y*cot(th) + 1/6;
    %
    idx_R = xs < data.X;
    idx_L = xs >= data.X;
    %
    rho = zeros(size(data.X));
    u   = zeros(size(data.X));
    v   = zeros(size(data.X));
    p   = zeros(size(data.X));
    %
    rho(idx_R) = rhoR;
    u(idx_R)   = uR;
    v(idx_R)   = vR;
    p(idx_R)   = pR;
    %
    rho(idx_L) = rhoL;
    u(idx_L)   = uL;
    v(idx_L)   = vL;
    p(idx_L)   = pL;
    %
    %
    data.Q(1,:,:,:) = rho;
    data.Q(2,:,:,:) = rho.*u;
    data.Q(3,:,:,:) = rho.*v;
    data.Q(4,:,:,:) = 0;
    data.Q(5,:,:,:) = p/(gam-1) + 0.5*rho.*(u.^2 + v.^2);
    %
end