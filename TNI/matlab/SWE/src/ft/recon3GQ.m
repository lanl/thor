function Qq=recon5GQ(Qq,data,dir)
    %
    weps = data.weno.weps;

    % get size
    Nx_total = data.Nx_total;
    Ny_total = data.Ny_total;
    %
    if(dir==1)
        d = 2;
    elseif(dir==2)
        d = 1;
    end
    %
    i = [1+2*(1==dir)+2*(1==d):Nx_total-2*(1==dir)-2*(1==d)];
    j = [1+2*(2==dir)+2*(2==d):Ny_total-2*(2==dir)-2*(2==d)];
    % calculate shift indices
    im = i - (1==d); ip = i + (1==d);
    jm = j - (2==d); jp = j + (2==d);
    
    % conserved variable-based reconstruction   
    Qm  = Qq(:,im,jm,1);
    Q   = Qq(:,i,j,1);
    Qp  = Qq(:,ip,jp,1);
    %
    params = data.reconParams(Qm,Q,Qp,weps);
    %
    for ii=1:size(data.weno.gamma,2)
        %
        c0 = data.weno.gamma(1,ii); 
        c1 = data.weno.gamma(2,ii); 
        %
        [Omega0,Omega1] = data.calcOmega(params,c0,c1);
        %
        c = data.weno.C{ii};
        %
        Q0 = c(1,1)*Q   + c(1,2)*Qp;
        Q1 = c(2,1)*Qm  + c(2,2)*Q ;
        %
        if(ii<=data.quad.n)
            %
            Qq(:,i,j,ii) = Q0.*Omega0+ Q1.*Omega1;
            %
        else
            %
            imid = data.weno.imid;
            %
            Qq(:,i,j,imid) = data.weno.sigma_p*Qq(:,i,j,imid) ...
                           - data.weno.sigma_m*(Q0.*Omega0+ Q1.*Omega1);
            %
        end
        %
    end
    %
end