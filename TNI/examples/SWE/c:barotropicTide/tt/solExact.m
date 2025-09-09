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
    %
    etah1 = 0.2;
    etah2 = 2*etah1;
    %
    k1 = 2.5*pi/Lx;
    k2 = 4.5*pi/Lx;
    %
    [eta1,u1,v1]=calculateVars(g,f,H,etah1,x,t,k1);
    [eta2,u2,v2]=calculateVars(g,f,H,etah2,x,t,k2);
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
function [eta,u,v]=calculateVars(g,f,H,etah,x,t,k)
    %
    x{1} = k*x{1};
    %
    om    = sqrt(g*H*k^2 + f^2);
    %
    coeff = g*etah*om*k/(om^2-f^2);
    %
    coskx = costt(x,1);
    sinkx = sintt(x,1);
    %
    cosot = cos(om*t);
    sinot = sin(om*t);
    %
    eta = etah*cosot*coskx;
    u   = coeff*sinot*sinkx;
    v   = (coeff*f/om*cosot)*sinkx;
    %
end