function Q=zero_ghost(Q,d)
    %
    for i=1:length(Q)
        %
        cr = Q{i}{d};
        %
        cr(:,[1:3],:)       = 0;
        cr(:,[end-2:end],:) = 0;
        %
        Q{i}{d} = cr;
        %
    end 
    %
end