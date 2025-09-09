function Q=BCisentropicVortex(Q,x,y,z,t,gam)
    %
    % Exact BCs in x and y
    %
    Qbc = solExact(x,y,z,t,gam);
    %
    Q(:,      1:3,:,:) = Qbc(:,      1:3,:,:);
    Q(:,end-2:end,:,:) = Qbc(:,end-2:end,:,:);
    Q(:,:,      1:3,:) = Qbc(:,:,      1:3,:);
    Q(:,:,end-2:end,:) = Qbc(:,:,end-2:end,:);
    %
    % Periodic in z
    %
    Q(:,:,:,      1:3) = Q(:,:,:,end-5:end-3);
    Q(:,:,:,end-2:end) = Q(:,:,:,        4:6);
    %
end