function Q=solExact(x,y,t,data)
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
    x{1} = x{1}*Lref;
    y{2} = y{2}*Lref;
    t    = t*tref;
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
    [eta1,u1,v1]=calculateVars(g,f,eta1,x,y,t,om1,kx1,ky1,data.tt.eps_rk);
    [eta2,u2,v2]=calculateVars(g,f,eta2,x,y,t,om2,kx2,ky2,data.tt.eps_rk);
    %
    Q      = cell(1,3);
    Q_quad = cell(1,3);
    %
    Q_quad{1} = round((eta1 + eta2)/data.Qref(1), data.tt.eps_rk(1));
    Q_quad{2} = round((u1   + u2  )/data.Qref(2), data.tt.eps_rk(2));
    Q_quad{3} = round((v1   + v2  )/data.Qref(3), data.tt.eps_rk(3));
    %
    for i=1:data.tt.Neq
        Q{i}    = applyQuadrature(Q_quad{i},data.quad,data.h);
        Q{i}    = round(Q{i},data.tt.eps_rk(i));
        Q{i}{1} = Q{i}{1}/data.vol;
    end
    %
end
%
function [eta,u,v]=calculateVars(g,f,etah,x,y,t,om,kx,ky,eps_rk)
    %
    coeff = g*etah/(om^2-f^2);
    %
    cosp = cosaxbyct(x,y,t,kx,ky,-om,min(eps_rk));
    sinp = sinaxbyct(x,y,t,kx,ky,-om,min(eps_rk));
    %
    eta = etah*cosp;
    u   = coeff*round(om*kx*cosp - f*ky*sinp,eps_rk(2));
    v   = coeff*round(om*ky*cosp + f*kx*sinp,eps_rk(3));
    %
end