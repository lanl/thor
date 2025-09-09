
function Q=remove_ghost(Q,Neq)
    %
    % remove ghost cells
    %
    for eq=1:Neq
        %
        tt = create_tt(Q{eq}.n-6,Q{eq}.r,Q{eq}.d);
        %
        for d=1:tt.d
            %
            core = Q{eq}{d};
            %
            tt{d} = core(:,4:Q{eq}.n(d)-3,:);
            %
        end
        %
        Q{eq} = tt;
    end
    %
end