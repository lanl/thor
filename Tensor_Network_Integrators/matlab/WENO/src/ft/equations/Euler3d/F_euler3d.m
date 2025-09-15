function F=F_euler3d(Q,dir,gam)
    % allocate memory
    F = zeros(size(Q));
    %
    rho  = Q(1,:);
    rhoU = Q(2,:);
    rhoV = Q(3,:);
    rhoW = Q(4,:);
    rhoE = Q(5,:);
    %
    U = rhoU./rho;
    V = rhoV./rho;
    W = rhoW./rho;
    %
    rhoe = rhoE - 0.5*(rhoU.*U + rhoV.*V + rhoW.*W);
    p    = (gam-1)*rhoe;
    p(p<1e-13) = 1e-13; 
    %
    if(min(p(:))<0) 
        error("Negative pressure")
    end
    %
    if(min(rho(:))<0) 
        error("Negative density")
    end
    %
    if(dir==1)
        %
        F(1,:) = rhoU(:);
        F(2,:) = rhoU(:).*U(:) + p(:);
        F(3,:) = rhoV(:).*U(:);
        F(4,:) = rhoW(:).*U(:);
        F(5,:) = U(:).*(rhoE(:) + p(:));
        %
    elseif(dir==2)
        %
        F(1,:) = rhoV(:);
        F(2,:) = rhoU(:).*V(:);
        F(3,:) = rhoV(:).*V(:) + p(:);
        F(4,:) = rhoW(:).*V(:);
        F(5,:) = V(:).*(rhoE(:) + p(:));
        %
    else
        %
        F(1,:) = rhoW(:);
        F(2,:) = rhoU(:).*W(:);
        F(3,:) = rhoV(:).*W(:);
        F(4,:) = rhoW(:).*W(:) + p(:);
        F(5,:) = W(:).*(rhoE(:) + p(:));
        %
    end
