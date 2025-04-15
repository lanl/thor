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
end
%
function Q = funSolExact(xx,i,t,gam)
    %
    % xx : xyz in tt format
    % i  : equation index
    % t  : current time 
    % gam: specific heat ratio
    %
    [~,Q] = solManuf(xx,i,t,gam);
    %
end