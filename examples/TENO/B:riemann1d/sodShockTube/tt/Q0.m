function data=Q0(data)
    %
    gam = data.gam;
    %
    x0    = 0.5;    % initial discontinuity
    p_r   = 0.1;    % right BC for pressure 
    u_r   = 0;      % right BC for velocity
    rho_r = 0.125;  % right BC for density
    p_l   = 1;      % left BC for pressure
    u_l   = 0;      % left BC for velocity
    rho_l = 1;      % left BC for density
    %
    n = [data.Nx_total;data.Ny_total;data.Nz_total];
    %
    tt_l = tt_ones(n);
    tt_r = tt_ones(n);
    %
    crX_l = tt_l{1};
    crX_r = tt_r{1};
    %
    crX_l(:,data.xg  >x0,:) = 0;
    crX_r(:,data.xg <=x0,:) = 0;
    %
    tt_l{1} = crX_l;
    tt_r{1} = crX_r;
    %
    data.tt.Q{1} = round(rho_l*tt_l + rho_r*tt_r,            data.tt.eps); 
    data.tt.Q{5} = round(p_l/(gam-1)*tt_l + p_r/(gam-1)*tt_r,data.tt.eps); 
    %
    data.tt.Q{2} = tt_zeros(n);
    data.tt.Q{3} = tt_zeros(n);
    data.tt.Q{4} = tt_zeros(n);
    %
end