
function [y,err,k]=taylor_inverse(x,eps_tt,Nsteps)
    %
    a    = 1.00;
    xp   = tt_ones(x.n);
    %
    dof = prod(x.n);
    %
    xavg = sum(x)/dof;
    %
    dof = sqrt(dof);
    %
    x   = x/xavg;
    xt  = round(1-x,eps_tt);
    %
    y  = tt_ones(x.n);
    %
    for k=1:Nsteps
        %
        xp = round(xp.*xt,eps_tt);
        %
        err = norm(xp)/dof;
        %
        if(err<eps_tt)
            break;
        end
        %
        y = round(y + xp,eps_tt);
        %
    end
    %
    y = y/xavg;
    %
end