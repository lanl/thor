function data=calculateEps(data)
    %
    Q = data.tt.Q;
    %
    Q=zero_ghost(Q,1);
    Q=zero_ghost(Q,2);
    Q=zero_ghost(Q,3);
    %
    nrm = 0;
    %
    h  = data.hmin;
    p  = max([1,data.hp-1.5]);
    hp = h^p;
    %
    for eq=1:data.tt.Neq
        %
        nrm_eq = norm(Q{eq});
        nrm    = max([nrm nrm_eq]);
        %
    end
    %
    data.tt.eps    = set_eps_tt(data.total_vol,data.tt.C_eps,hp,nrm);
    data.tt.eps_cr = data.tt.eps;
    %
end
%
function eps_tt = set_eps_tt(total_vol,C_eps,hp,nrm)
    %
    eps_tt    = min([1e-3,max([1e-14 sqrt(total_vol)*C_eps*hp/nrm])]);
    %
end