function data=Q0(data)
    %    
    gam  = data.gam;
    th   = pi/3;
    %
    rhoL =  8;
    uL   =  8.25*sin(th);
    vL   = -8.25*cos(th);
    pL   =  116.5; 
    %
    rhoR = 1.4;
    uR   = 0;
    vR   = 0;
    pR   = 1; 
    %
    d=3;
    n=[data.Nx_total;data.Ny_total;data.Nz_total];
    %
    y0   = tt_zeros(n);
    left = y0;
    th   = pi/3;
    %
    k1 = (data.Nx_total + data.Ny_total) + 1;
    k2 = sum(n);
    %
    for j=1:left.n(2)
        %
        y0 = tt_zeros(n);
        %
        xs = 10*data.curr_time/sin(th) + data.yg(j)*cot(th) + 1/6;
        %    
        y0.core(xs >= data.xg) = 1;
        y0.core(j+y0.ps(2)-1)  = 1;
        %
        y0.core(k1:k2) = 1;
        %
        left = left + y0;
        %
        left = round(left, data.tt.eps);
        %
    end 
    %
    rho = round((rhoL - rhoR)*left + rhoR,data.tt.eps); 
    u   = round((  uL -   uR)*left +   uR,data.tt.eps); 
    v   = round((  vL -   vR)*left +   vR,data.tt.eps); 
    p   = round((  pL -   pR)*left +   pR,data.tt.eps); 
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
