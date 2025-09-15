function data = BC(data,t)
    % periodic in x
    data.Q(:,1:3,:)      = data.Q(:,end-5:end-3,:);
    data.Q(:,end-2:end,:)= data.Q(:,4:6,:);
    % periodic in y
    data.Q(:,:,1:3)      = data.Q(:,:,end-5:end-3);
    data.Q(:,:,end-2:end)= data.Q(:,:,4:6);
    %
end
  