function data = BC(data,t)
    % no flow in y
    data.Q(1:2,:,1:3)       = data.Q(1:2,:,[6 5 4]);
    data.Q(1:2,:,end-2:end) = data.Q(1:2,:,[end-3 end-4 end-5]);
    %
    data.Q(3,:,1:3)       = -data.Q(3,:,[6 5 4]);
    data.Q(3,:,end-2:end) = -data.Q(3,:,[end-3 end-4 end-5]);
    % periodic in x
    data.Q(:,1:3,:)      = data.Q(:,end-5:end-3,:);
    data.Q(:,end-2:end,:)= data.Q(:,4:6,:);
    %
end
  