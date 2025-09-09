function Q=BC(Q,x,y,z,t,gam)
    %
    Q(:,1:3,:,:)      = Q(:,[4 4 4],:,:);
    Q(:,end-2:end,:,:)= Q(:,[end-3 end-3 end-3],:,:);
    %
    Q(:,:,1:3,:)       = Q(:,:,[4 4 4],:);
    Q(:,:,end-2:end,:) = Q(:,:,[end-3 end-3 end-3],:);
    %
    Q(:,:,:,1:3)       = Q(:,:,:,[4 4 4]);
    Q(:,:,:,end-2:end) = Q(:,:,:,[end-3 end-3 end-3]);
    %
end