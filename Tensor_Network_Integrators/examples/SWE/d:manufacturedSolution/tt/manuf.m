function data=manuf(x,y,t,data)
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
    etahat = 0.01;
    %
    kx     = 2*pi/Lx;
    ky     = 2*pi/Ly;
    c      = sqrt(g*H);
    omega  = c*sqrt(kx^2+ky^2);
    %
    %phi = kx*x + ky*y - omega*t;
    %
    cosp = etahat*cosaxbyct(x,y,t,kx,ky,-omega,2e-16);
    sinp = etahat*sinaxbyct(x,y,t,kx,ky,-omega,2e-16);
    %
    cosp_phi = -sinp;
    sinp_phi =  cosp;
    %
    phi_t = -omega;
    phi_x = kx;
    phi_y = ky;
    %
    cosp_t = cosp_phi*phi_t;
    cosp_x = cosp_phi*phi_x;
    cosp_y = cosp_phi*phi_y;
    %
    sinp_t = sinp_phi*phi_t;
    sinp_x = sinp_phi*phi_x;
    sinp_y = sinp_phi*phi_y;
    %
    U   = cosp;
    U_t = cosp_t;
    U_x = cosp_x;
    U_y = cosp_y;
    %
    V   = tt_zeros(x.n);
    V_t = V;
    V_x = V;
    V_y = V;
    %
    h   =  sinp + H;
    h_t =  sinp_t;
    h_x =  sinp_x;
    h_y =  sinp_y;
    %
    hU   = h.*U;
    hU_t = h_t.*U + h.*U_t;
    hU_x = h_x.*U + h.*U_x;
    hU_y = h_y.*U + h.*U_y;
    %
    hV   = h.*V;
    hV_t = h_t.*V + h.*V_t;
    hV_x = h_x.*V + h.*V_x;
    hV_y = h_y.*V + h.*V_y;
    %
    F_x = hU_x.*U + hU.*U_x + g*h.*h_x;
    F_y = hU_y.*V + hU.*V_y;
    %
    G_x = hU_x.*V + hU.*V_x;
    G_y = hV_y.*V + hV.*V_y + g*h.*h_y;
    %
    % calculate the time derivative of the conserved variable vector
    %
    data.tt.Q0{1} = round( h/data.Qref(1),data.tt.eps_rk(1));
    data.tt.Q0{2} = round(hU/data.Qref(2),data.tt.eps_rk(2));
    data.tt.Q0{3} = round(hV/data.Qref(3),data.tt.eps_rk(3));
    %
    data.tt.sourceQ{1} = round((h_t + hU_x + hV_y      )*(tref)/(Href)     , data.tt.eps);
    data.tt.sourceQ{2} = round((hU_t + F_x + F_y - f*hV)*(tref)/(Href*Uref), data.tt.eps);
    data.tt.sourceQ{3} = round((hV_t + G_x + G_y + f*hU)*(tref)/(Href*Uref), data.tt.eps);
    %
    %
    for i=1:data.tt.Neq
        %
        data.tt.source{i} = applyQuadrature(data.tt.sourceQ{i},data.quad,data.h);
        %
    end
    %
end