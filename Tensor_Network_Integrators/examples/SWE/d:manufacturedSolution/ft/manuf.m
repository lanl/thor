function [S,Q] = manuf(x,y,t,data)
    %
    %  [S,Q]: nondimensional source term and conserved variables
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
    etahat = 0.01;
    %
    kx     = 2*pi/Lx;
    ky     = 2*pi/Ly;
    c      = sqrt(g*H);
    omega  = c*sqrt(kx^2+ky^2);
    %
    phi = kx*x + ky*y - omega*t;
    %
    cosp = cos(phi);
    sinp = sin(phi);
    %
    cosp_phi = -sinp;
    sinp_phi =  cosp;
    %
    phi_t = -omega;
    phi_x = kx;
    phi_y = ky;
    %
    cosp_t = cosp_phi.*phi_t;
    cosp_x = cosp_phi.*phi_x;
    cosp_y = cosp_phi.*phi_y;
    %
    sinp_t = sinp_phi.*phi_t;
    sinp_x = sinp_phi.*phi_x;
    sinp_y = sinp_phi.*phi_y;
    %
    U   = etahat*cosp;
    U_t = etahat*cosp_t;
    U_x = etahat*cosp_x;
    U_y = etahat*cosp_y;
    %
    V   = zeros(size(U));
    V_t = V;
    V_x = V;
    V_y = V;
    %
    h   =  etahat*sinp + H;
    h_t =  etahat*sinp_t;
    h_x =  etahat*sinp_x;
    h_y =  etahat*sinp_y;
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
    S = zeros([3,size(x)]);
    Q = zeros([3,size(x)]);
    %
    Q(1,:) = h(:)/Href;
    Q(2,:) = hU(:)/(Href*Uref);
    Q(3,:) = hV(:)/(Href*Uref);
    %
    S(1,:) = (h_t(:) + hU_x(:) + hV_y(:))*(tref)/(Href);
    %
    S(2,:) = (hU_t(:) + F_x(:) + F_y(:) - f*hV(:))*(tref)/(Href*Uref);
    %
    S(3,:) = (hV_t(:) + G_x(:) + G_y(:) + f*hU(:))*(tref)/(Href*Uref);
    %
end