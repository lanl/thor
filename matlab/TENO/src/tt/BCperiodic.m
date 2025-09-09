function Q=BCperiodic(Q,data,t,d)
    %
    % apply periodic BCs in the d-direction
    %
    % loop over equations
    for i=1:length(Q)
        %
        cr = Q{i}{d};
        %
        cr(:,1:3,:)       = cr(:,end-5:end-3,:);
        cr(:,end-2:end,:) = cr(:,4:6,:);
        %
        Q{i}{d} = cr;
        %
    end
    %
end