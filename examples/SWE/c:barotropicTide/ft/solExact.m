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
    etah1 = 0.2;
    etah2 = 2*etah1;
    %
    k1 = 2.5*pi/Lx;
    k2 = 4.5*pi/Lx;
    %
    [eta1,u1,v1]=calculateVars(g,f,H,etah1,k1,x,t);
    [eta2,u2,v2]=calculateVars(g,f,H,etah2,k2,x,t);
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
function [eta,u,v]=calculateVars(g,f,H,etah,k,x,t)
    %
    om = sqrt(g*H*k^2 + f^2);
    %
    coeff = g*etah*om*k/(om^2-f^2);
    %
    coskx = cos(k*x);
    sinkx = sin(k*x);
    %
    cosot = cos(om*t);
    sinot = sin(om*t);
    %
    eta =  etah*coskx.*cosot;
    u   = coeff*sinkx.*sinot;
    v   = coeff*f/om*sinkx.*cosot;
    %
end
