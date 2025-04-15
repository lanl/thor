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
    eta1 = 1e-4;
    eta2 = 2*eta1;
    %
    ky1 = 2*pi/Ly;
    ky2 = 2*ky1;
    %
    c = sqrt(g*H);
    R = c/f;
    %
    x.core{1} = -x.core{1}/R;
    y.core{2} = y.core{2} + c*t;
    %
    ypct1 = y;
    ypct2 = y;
    %
    ypct1.core{2} = ky1*ypct1.core{2};
    ypct2.core{2} = ky2*ypct2.core{2};
    %
    Fy = eta1*sintt(ypct1,2) + eta2*sintt(ypct2,2);
    %
    Q      = cell(1,3);
    Q_quad = cell(1,3);
    %
    Q_quad{1} = -round(H*Fy.*exptt(x,1)/data.Qref(1),data.tt.eps_rk(1));
    Q_quad{2} =  tt_zeros(data.tt.N*data.quad.n);
    Q_quad{3} =  round(c*Fy.*exptt(x,1)/data.Qref(3),data.tt.eps_rk(3));
    %
    for i=1:data.tt.Neq
        Q{i}    = applyQuadrature(Q_quad{i},data.quad,data.h);
        Q{i}    = round(Q{i},data.tt.eps_rk(i));
        Q{i}{1} = Q{i}{1}/data.vol;
    end
    %
end