function Q=BC(Q,data,t,d)
    % loop over equations
    for i=1:length(Q)
        %
        cr = Q{i}{d};
        %
        cr(:,1:3,:)       = cr(:,[4 4 4],:);
        cr(:,end-2:end,:) = cr(:,[end-3 end-3 end-3],:);
        %
        Q{i}{d} = cr;
        %
    end
    %
end