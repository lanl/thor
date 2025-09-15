function Q = solExact(data,t,Q0)
    %
    if(~exist("Q0","var"))
        Q0 = cell(data.tt.Neq,1); 
    end
    %
    Q = cell(1,data.tt.Neq);
    for i=1:data.tt.Neq
        Q{i} = cross_interpolation(data.tt.xyz, @(xxx) funSolExact(xxx,i,t,data.gam), data.tt.eps_cr, Q0{i}, 0);
    end
    %
    % Since amen_cross fails for identically zero vectors, w was set to 1 uniformly in cross_interpolation:
    % Set w-velocity to zero to fix this
    %
    Q{4}.core(:) = 0;  
    %
end
%
function Q = funSolExact(xx,i,t,gam)
    %
    % xx : xyz in tt format
    % i  : equation index
    % t  : current time 
    % gam: specific heat ratio
    %
    beta   = 5;
    uinf   = 1;
    vinf   = 1;
    x0     = 5;
    y0     = 5;
    %
    xbar = xx(:,1) - (x0 + uinf*t); 
    ybar = xx(:,2) - (y0 + vinf*t); 
    %
    r2 = xbar.^2 + ybar.^2;
    %
    exp1 =  exp(1-r2);
    %
    exp2    = exp(0.5*(1-r2));
    %
    rho_coeff   = 0.125*(gam-1)/gam*(beta^2)/(pi^2);
    rho_term    = 1 - rho_coeff*exp1;
    %
    ucoeff = 0.5*beta/pi;
    uterm  = ucoeff*exp2;
    %
    rho_pw = 1/(gam-1);
    rho    = rho_term.^rho_pw;
    %
    U = uinf - uterm.*ybar;
    V = vinf + uterm.*xbar;
    W = ones(size(rho)); % so that amen_cross does not fail
    %
    % calculate the conserved variable vector
    %
    if(i==1)
        % 
        Q = rho; 
        %
    elseif(i==2)
        %
        Q = rho.*U;  
        %
    elseif(i==3)
        %
        Q = rho.*V;
        %
    elseif(i==4)
        %
        Q = rho.*W;
        %
    elseif(i==5)
        %
        p = rho.^gam;
        %
        rhoe = p/(gam-1);
        %
        u2pv2 = U.^2 + V.^2; % (ignore w)
        %
        Et = rhoe + 0.5*rho.*u2pv2;
        %
        Q = Et;
        %
    end
    %
end