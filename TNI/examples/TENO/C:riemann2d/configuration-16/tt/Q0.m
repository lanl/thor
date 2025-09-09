function data=Q0(data)
    %
    gam = data.gam;
    %
    rho1 = 0.5313;
    u1   = 0.1;
    v1   = 0.1;
    p1   = 0.4;
    %
    rho2 = 1.0222;
    u2   = -0.6179;
    v2   = 0.1;
    p2   = 1;
    %
    rho3 = 0.8;
    u3   = 0.1;
    v3   = 0.1;
    p3   = 1;
    %
    rho4 = 1;
    u4   = 0.1;
    v4   = 0.8276;
    p4   = 1;
    %
    d=3;
    n=[data.Nx_total;data.Ny_total;data.Nz_total];
    %
    y1   = tt_zeros(n);
    y2   = tt_zeros(n);
    y3   = tt_zeros(n);
    y4   = tt_zeros(n);
    %
    x = data.xg;
    y = data.yg;
    %
    idx1x = x >= 0.5;
    idx1y = y >= 0.5;
    %
    idx2x = x  < 0.5;
    idx2y = y >= 0.5;
    %
    idx3x = x  < 0.5;
    idx3y = y  < 0.5;
    %
    idx4x = x >= 0.5;
    idx4y = y  < 0.5;
    %
    core1x = y1{1};
    core1y = y1{2};
    %
    core2x = y2{1};
    core2y = y2{2};
    %
    core3x = y3{1};
    core3y = y3{2};
    %
    core4x = y4{1};
    core4y = y4{2};
    %
    core1x(:,idx1x,:) = 1;
    core1y(:,idx1y,:) = 1;
    %
    core2x(:,idx2x,:) = 1;
    core2y(:,idx2y,:) = 1;
    %
    core3x(:,idx3x,:) = 1;
    core3y(:,idx3y,:) = 1;
    %
    core4x(:,idx4x,:) = 1;
    core4y(:,idx4y,:) = 1;
    %
    y1{1} = core1x;
    y1{2} = core1y;
    %
    y2{1} = core2x;
    y2{2} = core2y;
    %
    y3{1} = core3x;
    y3{2} = core3y;
    %
    y4{1} = core4x;
    y4{2} = core4y;
    %
    corez1 = y1{3};
    corez2 = y2{3};
    corez3 = y3{3};
    corez4 = y4{3};
    %
    corez1(:) = 1;
    corez2(:) = 1;
    corez3(:) = 1;
    corez4(:) = 1;
    %
    y1{3} = corez1;
    y2{3} = corez2;
    y3{3} = corez3;
    y4{3} = corez4;
    %
    rho = y1*rho1 + y2*rho2 + y3*rho3 + y4*rho4; 
    u   = y1*u1   + y2*u2   + y3*u3   + y4*u4  ; 
    v   = y1*v1   + y2*v2   + y3*v3   + y4*v4  ; 
    p   = y1*p1   + y2*p2   + y3*p3   + y4*p4  ; 
    %
    rho = round(rho,data.tt.eps);
    u   = round(u,data.tt.eps);
    v   = round(v,data.tt.eps);
    p   = round(p,data.tt.eps);
    %
    rhou = round(rho.*u,data.tt.eps);
    rhov = round(rho.*v,data.tt.eps);
    %
    data.tt.Q{1} = rho;
    data.tt.Q{2} = rhou;
    data.tt.Q{3} = rhov;
    %
    data.tt.Q{4} = tt_zeros(n);
    %
    data.tt.Q{5} = round(round(0.5*rhou.*u,data.tt.eps) + round(0.5*rhov.*v,data.tt.eps),data.tt.eps);
    data.tt.Q{5} = round(data.tt.Q{5} + p/(gam-1),data.tt.eps);
    %
end