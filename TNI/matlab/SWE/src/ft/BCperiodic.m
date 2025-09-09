function Q=BCperiodic(Q,x,y,t,gacc)
    %
    % apply periodic BC
    %
    Q(:,1:3,:)      = Q(:,end-5:end-3,:);
    Q(:,end-2:end,:)= Q(:,4:6,:);
    %
    Q(:,:,1:3)      = Q(:,:,end-5:end-3);
    Q(:,:,end-2:end)= Q(:,:,4:6);
    