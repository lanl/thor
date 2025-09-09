function data=Q0(data)
    %
    a   = 0.5;
    b   = 0.5;
    k   = 2*pi;
    %
    arg = k*data.X;
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
        arg = k*(data.X-u*t);
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
        %
        u(:,:,:) = u(:,:,:) - 0.5*f./fp;
    end
    %
    if(err>=tol)
        error("Newton solver could not converge. The error at the last iteration was %e",err);
    end
    %
    data.Q = reshape(u,[1,size(u)]);
    %
end
