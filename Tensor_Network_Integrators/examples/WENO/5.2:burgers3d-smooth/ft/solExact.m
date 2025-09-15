function u = solExact(x,y,z,t,gam)
    %
    a   = 0.5;
    b   = 0.5;
    k   = 2*pi;
    %
    arg = k*(x + y + z);
    % initial guess
    u = a+b*sin(arg); 
    %
    err      = 1;
    tol      = 1e-12;
    max_iter = 100000;
    %
    for iter=1:max_iter
        %
        arg = k*(x + y + z - 3*u*t);
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
        fp = 1 + 3*t*du0;
        %
        u(:,:,:) = u(:,:,:) - 0.5*f./fp;
    end
    %
    if(err>=tol)
        error("Newton solver could not converge. The error at the last iteration was %e",err);
    end
    %
    u = reshape(u,[1,size(u)]);
    %
end
%