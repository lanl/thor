function F=F_linearSWE(Q,data,d)
    %
    H = data.H;
    g = data.gacc;
    %
    F = cell(1,3);
    %
    if(d==1)
        %
        F{1} = H*Q{2};
        F{2} = g*Q{1};
        F{3} = tt_zeros(Q{1}.n);
        %
    elseif(d==2)
        %
        F{1} = H*Q{3};
        F{2} = tt_zeros(Q{1}.n);
        F{3} = g*Q{1};
        %
    end
    %
end

