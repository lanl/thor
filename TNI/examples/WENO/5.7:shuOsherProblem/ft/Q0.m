function data=Q0(data)
    %
    gam = data.gam;
    %
    rho = zeros(size(data.X));
    u   = zeros(size(data.X));
    p   = zeros(size(data.X));
    %
    xs = data.X - 5;
    %
    idxL = xs<-4;
    idxR = xs>=-4;
    %
    rho(idxL) = 27/7;
    u(idxL)   = 4*sqrt(35)/9;
    p(idxL)   = 31/3;
    %
    rho(idxR) = 1 + 0.2*sin(5*xs(idxR));
    u(idxR)   = 0;
    p(idxR)   = 1;
    %
    data.Q(1,:,:,:) = rho; 
    data.Q(2,:,:,:) = rho.*u;
    data.Q(5,:,:,:) = p/(gam-1) + 0.5*rho.*u.^2;
    %
end