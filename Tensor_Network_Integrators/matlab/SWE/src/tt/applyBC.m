function Q=applyBC(Q,data,t)
    %
    for d=1:data.tt.d
        Q=data.BC(Q,data,t,d);
    end
    %
    for i=1:data.tt.Neq
        Q{i}=round(Q{i},data.tt.eps_rk(i));
    end 
    %
end