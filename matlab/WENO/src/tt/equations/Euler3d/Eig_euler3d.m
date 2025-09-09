function Eig = Eig_euler3d(Q,d,gam)
    % allocate memory
    Eig = zeros(size(Q));
    %
    rho  = Q(:,1);
    rhoU = Q(:,2);
    rhoV = Q(:,3);
    rhoW = Q(:,4);
    rhoE = Q(:,5);
    %
    U = rhoU./rho;
    V = rhoV./rho;
    W = rhoW./rho;
    %
    rhoe = rhoE - 0.5*rhoU.*U - 0.5*rhoV.*V - 0.5*rhoW.*W;
    p    = (gam-1)*rhoe;
    %
    if(min(p(:))<0) 
        error("Negative pressure")
    end
    %
    if(min(rho(:))<0) 
        error("Negative density")
    end
    % calculate speed of sound
    a = sqrt(gam*p./rho);
    if(d==1)
        temp = abs(U) + a;
    elseif(d==2)
        temp = abs(V) + a;
    elseif(d==3)
        temp = abs(W) + a;
    else
        temp = max(max(abs(U),abs(V)),abs(W)) + a;
    end
    %
    Eig(:,1) = temp(:);
    Eig(:,2) = temp(:);
    Eig(:,3) = temp(:);
    Eig(:,4) = temp(:);
    Eig(:,5) = temp(:);