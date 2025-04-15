function data=Q0(data)
    %
    gam = data.gam;
    %
    X = data.tt.xyz{1};
    %
    n = data.tt.xyz{1}.n;
    %
    x = X{1} - 5;
    %
    rho = tt_ones(n);
    u   = tt_ones(n);
    p   = rho;
    %
    idxL = x<-4;
    idxR = x>=-4;
    %
    rho.core(idxL) = 27/7;
    u.core(idxL)   = 4*sqrt(35)/9;
    p.core(idxL)   = 31/3;
    %
    rho.core(idxR) = 1 + 0.2*sin(5*x(idxR));
    u.core(idxR)   = 0;
    p.core(idxR)   = 1;
    %
    rhou  = round(rho.*u,data.tt.eps);
    rhou2 = round(rhou.*u,data.tt.eps);
    %
    data.tt.Q{1} = rho;
    data.tt.Q{2} = rhou;
    data.tt.Q{3} = tt_zeros(n);
    data.tt.Q{4} = tt_zeros(n);
    data.tt.Q{5} = round(p/(gam-1) + 0.5*rhou2,data.tt.eps);
    %
    data.tt.Qbc = data.tt.Q; % needed to enforce boundary conditions
    %
end