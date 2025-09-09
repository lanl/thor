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
    eta1 = 1e-1;
    eta2 = 2*eta1;
    %
    kx1 = 2*pi/Lx;
    kx2 = 2*kx1;
    %
    ky1 = 2*pi/Ly;
    ky2 = 2*ky1;
    %
    k1 = kx1^2 + ky1^2;
    k2 = kx2^2 + ky2^2;
    %
    c = sqrt(g*H);
    R = c/f;
    %
    om1 = sqrt(c^2*k1 + f^2);
    om2 = sqrt(c^2*k2 + f^2);
    %
    phi1 = kx1*x + ky1*y - om1*t;
    phi2 = kx2*x + ky2*y - om2*t;
    %
    [eta1,u1,v1]=calculateVars(g,f,eta1,phi1,om1,kx1,ky1);
    [eta2,u2,v2]=calculateVars(g,f,eta2,phi2,om2,kx2,ky2);
    %
    % conserved variables
    %
    Q = zeros([3,size(x)]);
    %
    Q(1,:) = (eta1(:) + eta2(:))/Href;
    Q(2,:) = (u1(:)   + u2(:))/Uref;
    Q(3,:) = (v1(:)   + v2(:))/Uref;
    %
end
%
function [eta,u,v]=calculateVars(g,f,etah,phi,om,kx,ky)
    %
    coeff = g*etah/(om^2-f^2);
    %
    cosp = cos(phi);
    sinp = sin(phi);
    %
    eta = etah*cosp;
    u   = coeff*(om*kx*cosp - f*ky*sinp);
    v   = coeff*(om*ky*cosp + f*kx*sinp);
    %
end