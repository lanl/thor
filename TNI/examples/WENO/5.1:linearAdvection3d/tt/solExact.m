function u = solExact(data,t)
    %
    u = cell(1,data.tt.Neq);
    for i=1:data.tt.Neq
        u{i} = cross_interpolation(data.tt.xyz, @(xxx) funSolExact(xxx,i,t,[]), data.tt.eps_cr, [], 0);
    end
    %
end
%
function u = funSolExact(xx,i,t,gam)
    %
    % xx : xyz in tt format
    % i  : equation index
    % t  : current time 
    % gam: specific heat ratio
    %
    x=xx(:,1); y=xx(:,2); z=xx(:,3);
 
    n1 = 2*pi; n2 = 2*pi; n3=2*pi; w = 6*pi;
    %
    arg = n1*x + n2*y + n3*z - w*t;
    %
    u = sin(arg);
    %
end