function data=Q0(data)
    %
    for i=1:data.tt.Neq
        data.tt.Q{i} = cross_interpolation(data.tt.xyz, @(xxx) funQ0(xxx,i,[]), data.tt.eps_cr, [], 0);
    end
    %
end
%
function u = funQ0(xx,i,gam)
    %
    % xx : xyz in tt format
    % i  : equation index
    % T  : current time 
    % gam: specific heat ratio
    %
    x=xx(:,1); y=xx(:,2); z=xx(:,3);
    %
    a   = 0.5;
    b   = 0.5;
    k   = 2*pi;
    %
    arg = k*x;
    % initial guess
    u = a+b*sin(arg); 
    %
    err      = 1;
    tol      = 1e-12;
    max_iter = 10000;
    %
    t  = 0;
    %
    for iter=1:max_iter
        %
        arg = k*(x - u*t);
        %
        u0  = a+b*sin(arg);
        du0 = b*k*cos(arg);
        %
        f  = u - u0;
        %
        err = max(abs(f(:)));
        %
        if(err<tol)
            break;
        end
        %
        fp = 1 + t*du0;
        u  = u - 0.5*f./fp;
    end
    %
    if(err>=tol)
        error("Newton solver could not converge. The error at the last iteration was %e",err);
    end
    %
end