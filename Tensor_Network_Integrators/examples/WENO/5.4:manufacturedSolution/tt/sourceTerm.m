function data = sourceTerm(data)
    %
    for i=1:data.tt.Neq
        %
        data.tt.source{i} = cross_interpolation(data.tt.xyz, @(xxx) solManuf(xxx,i,data.rk_time,data.gam), data.tt.eps_cr, [], 0);
        %
        data.tt.RHS{i} = round(data.tt.RHS{i} + data.tt.source{i},data.tt.eps);
        %
    end
    %
end