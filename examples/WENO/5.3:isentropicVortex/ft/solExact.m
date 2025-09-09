function Q = solExact(x,y,z,t,gam)
    %
    beta   = 5;
    uinf   = 1;
    vinf   = 1;
    x0     = 5;
    y0     = 5;
    %
    xbar = x - (x0 + uinf*t); 
    ybar = y - (y0 + vinf*t); 
    %
    r2 = xbar.^2 + ybar.^2;
    %
    rho = (1 - 0.125*(gam-1)/gam*(beta^2)/(pi^2)*exp(1-r2)).^(1/(gam-1));
    u   = uinf - 0.5*beta/pi*exp(0.5*(1-r2)).*ybar;
    v   = vinf + 0.5*beta/pi*exp(0.5*(1-r2)).*xbar;
    p   = rho.^gam;
    %
    Q = zeros([5,size(x)]);
    Q(1,:,:,:) = rho;
    Q(2,:,:,:) = rho.*u;
    Q(3,:,:,:) = rho.*v;
    Q(5,:,:,:) = p/(gam-1) + 0.5*rho.*(u.^2 + v.^2);
    %
end
%