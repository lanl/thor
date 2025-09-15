function Q=BC(Q,data,t,d)
    if(d==1)
        % this is from the initialization step
        Qbc = data.tt.Qbc;
        %
        Qinterior = Qbc;
        % this is to apply BCs only in x
        Qinterior = zero_ghost(Qinterior,1); 
        %
        for j=1:3
            Q = zero_ghost(Q,j);
        end
        %
        % loop over equations to obtain BCs only in x (this will set ghost cells in y and z to zero)
        %
        for i=1:5
            Qbc{i} = round(Qbc{i} - Qinterior{i},eps);
        end
        %
        for i=1:5
            Q{i} = Q{i} + Qbc{i};
        end
    else
        % assuming periodic in y and z is OK
        Q=BCperiodic(Q,data,t,d);
    end
    %
end