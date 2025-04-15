function Q=BC(Q,data,t,d)
    %
    if(d==1) % only apply BCs for once
        %
        Qbc = solExact(data.tt.q_xy{1},data.tt.q_xy{2},t,data);
        %
        for j=1:2
            Q = zero_ghost(Q,j);
        end
        %
        %
        Qinterior = Qbc;
        %
        for j=1:2
            %
            Qinterior = zero_ghost(Qinterior,j);
            %
        end
        %
        % loop over equations
        %
        for i=1:length(Qbc)
            %
            Qbc{i} = round(Qbc{i} - Qinterior{i},data.tt.eps_rk(i));
            %
        end
        %
        Qbc = zero_ghost(Qbc,2);
        %
        for i=1:length(Q)
            Q{i} = Q{i} + Qbc{i};
        end
        %
        Q=BCperiodic(Q,data,t,2);
        %
    end
    %
end