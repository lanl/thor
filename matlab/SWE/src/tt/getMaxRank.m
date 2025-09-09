function [rmax,idx] = getMaxRank(Q)
    %
    rmax = 0;
    for eq=1:length(Q)
        %
        r = max(Q{eq}.r);
        %
        if(r>rmax)
            idx  = eq;
            rmax = r;
        end
        %
    end
    %
end