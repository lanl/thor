function Q=BC(Q,x,y,z,t,gam)
    %
    th = pi/3;
    Ly = 1;
    % 
    [temp,idx16]=min(abs(x(:,1,1)-1/6));
    if(x(idx16,1,1)>1/6) 
        idx16 = idx16-1;
    end
    %
    xs  = 10*t/sin(th) + 1/6 + Ly*cot(th);
    [temp,idxs]=min(abs(x(:,end,1)-xs));
    if(x(idxs,end,1)>xs) 
        idxs = idxs-1;
    end
    %
    rhoL =  8;
    uL   =  8.25*sin(th);
    vL   = -8.25*cos(th);
    pL   =  116.5; 
    %
    rhoR = 1.4;
    uR   = 0;
    vR   = 0;
    pR   = 1; 
    % inflow (supersonic)
    Q(1,1:3,:,:)      = rhoL;
    Q(2,1:3,:,:)      = rhoL*uL;
    Q(3,1:3,:,:)      = rhoL*vL;
    Q(5,1:3,:,:)      = pL/(gam-1) + 0.5*rhoL*(uL^2+vL^2);
    % outflow (subsonic)
    rhoin = Q(1,[end-3 end-3 end-3],:,:);
    uin   = Q(2,[end-3 end-3 end-3],:,:)./rhoin;
    vin   = Q(3,[end-3 end-3 end-3],:,:)./rhoin;
    rhoein= Q(5,[end-3 end-3 end-3],:,:) - 0.5*rhoin.*(uin.^2+vin.^2);
    pin   = (gam-1)*rhoein;
    pex   = 2*pR-pin;
    %
    Q(1,end-2:end,:,:) = Q(1,[end-3 end-3 end-3],:,:);
    Q(2,end-2:end,:,:) = Q(2,[end-3 end-3 end-3],:,:);
    Q(3,end-2:end,:,:) = Q(3,[end-3 end-3 end-3],:,:);
    Q(5,end-2:end,:,:) = pex/(gam-1) + 0.5*rhoin.*(uin.^2+vin.^2);
    % y=0 & x<1/6
    Q(1,1:idx16,1:3,:) = rhoL;
    Q(2,1:idx16,1:3,:) = rhoL*uL;
    Q(3,1:idx16,1:3,:) = rhoL*vL;
    Q(5,1:idx16,1:3,:) = pL/(gam-1) + 0.5*rhoL*(uL^2+vL^2);
    % y=0 & x>1/6 (symmetry)
    Q(1,idx16+1:end,1:3,:) =  Q(1,idx16+1:end,[6 5 4],:);
    Q(2,idx16+1:end,1:3,:) =  Q(2,idx16+1:end,[6 5 4],:);
    Q(3,idx16+1:end,1:3,:) = -Q(3,idx16+1:end,[6 5 4],:);
    Q(5,idx16+1:end,1:3,:) =  Q(5,idx16+1:end,[6 5 4],:);
    % y=Ly & x<xs
    Q(1,1:idxs,end-2:end,:) = rhoL;
    Q(2,1:idxs,end-2:end,:) = rhoL*uL;
    Q(3,1:idxs,end-2:end,:) = rhoL*vL;
    Q(5,1:idxs,end-2:end,:) = pL/(gam-1) + 0.5*rhoL*(uL^2+vL^2);
    % y=Ly & x>xs
    Q(1,idxs+1:end,end-2:end,:) = rhoR;
    Q(2,idxs+1:end,end-2:end,:) = rhoR*uR;
    Q(3,idxs+1:end,end-2:end,:) = rhoR*vR;
    Q(5,idxs+1:end,end-2:end,:) = pR/(gam-1) + 0.5*rhoR*(uR^2+vR^2);
    % no velocity in z-direction
    Q(4,:,:,:)         = 0;
    % periodic in z
    Q(:,:,:,1:3)       = Q(:,:,:,end-5:end-3);
    Q(:,:,:,end-2:end) = Q(:,:,:,4:6);
    %
end