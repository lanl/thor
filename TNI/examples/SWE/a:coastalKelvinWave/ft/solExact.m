function Q = solExact(x,y,t,data)
    %
    %  Q: nondimensional conserved variables
    %  x,y,t: nondimensional x/y coordinates and time 
    % 
    % first dimensionalize all nondimensional variables
    %
    g    = data.gacc*data.ref.g;
    Lref = data.ref.L;
    Href = data.ref.H;
    Uref = data.ref.U;
    tref = data.ref.t;
    %
    x = x*Lref;
    y = y*Lref;
    t = t*tref;
    %
    Lx = data.Lx * data.ref.L;
    Ly = data.Ly * data.ref.L;
    %
    H = data.H*data.ref.H;
    f = data.f*data.ref.f;
    %
    eta1 = 1e-4;
    eta2 = 2*eta1;
    %
    ky1 = 2*pi/Ly;
    ky2 = 2*ky1;
    %
    c = sqrt(g*H);
    R = c/f;
    %
    ypct = y + c*t;
    %
    Fy = eta1*sin(ky1*ypct) +  eta2*sin(ky2*ypct);
    %
    eta = -H*Fy.*exp(-x/R);
    v   = c*Fy.*exp(-x/R);
    %
    % conserved variables
    %
    Q = zeros([3,size(x)]);
    %
    Q(1,:) = eta(:)/data.Qref(1);
    Q(3,:) = v(:)/data.Qref(3);
    %
end