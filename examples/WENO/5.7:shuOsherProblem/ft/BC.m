function Q=BC(Q,x,y,z,t,gam)
    %
    rho_l = 27/7*ones(1,3);
    u_l   = 4*sqrt(35)/9*ones(1,3);
    p_l   = 31/3*ones(1,3);
    %
    rho_r = 1 + 0.2*sin(5*x(end-2:end,1,1));
    u_r   = zeros(1,3);
    p_r   = ones(1,3);
    %
    rhoE_l = p_l/(gam-1) + 0.5*rho_l.*u_l.^2;
    rhoE_r = p_r/(gam-1) + 0.5*rho_r.*u_r.^2;
    %
    for i=1:3
        %
        Q(1,1:3,:,:)     = rho_l(i);
        Q(1,end-3+i,:,:) = rho_r(i);
        %
        Q(2,1:3,:,:)     = rho_l(i)*u_l(i);
        Q(2,end-3+i,:,:) = rho_r(i)*u_r(i);
        %
        Q(3:4,1:3,:,:)     = 0;
        Q(3:4,end-3+i,:,:) = 0;
        %
        Q(5,1:3,:,:)     = rhoE_l(i);
        Q(5,end-3+i,:,:) = rhoE_r(i);
        %
    end
    %
    Q(:,:,1:3,:)       = Q(:,:,[4 4 4],:);
    Q(:,:,end-2:end,:) = Q(:,:,[end-3 end-3 end-3],:);
    %
    Q(:,:,:,1:3)       = Q(:,:,:,[4 4 4]);
    Q(:,:,:,end-2:end) = Q(:,:,:,[end-3 end-3 end-3]);
    %
end