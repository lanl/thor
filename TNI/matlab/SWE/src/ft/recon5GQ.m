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
    imm = i - 2*(1==d); im = i - (1==d); ip = i + (1==d); ipp = i + 2*(1==d);
    jmm = j - 2*(2==d); jm = j - (2==d); jp = j + (2==d); jpp = j + 2*(2==d);
    
    % conserved variable-based reconstruction   
    Qmm = Qq(:,imm,jmm,1);
    Qm  = Qq(:,im,jm,1);
    Q   = Qq(:,i,j,1);
    Qp  = Qq(:,ip,jp,1);
    Qpp = Qq(:,ipp,jpp,1);
    %
    params = data.reconParams(Qmm,Qm,Q,Qp,Qpp,weps);
    %
    for ii=1:size(data.weno.gamma,2)
        %
        c0 = data.weno.gamma(1,ii); 
        c1 = data.weno.gamma(2,ii); 
        c2 = data.weno.gamma(3,ii);
        %
        [Omega0,Omega1,Omega2] = data.calcOmega(params,c0,c1,c2);
        %
        c = data.weno.C{ii};
        %
        Q0 = c(1,1)*Q   + c(1,2)*Qp +  c(1,3)*Qpp;
        Q1 = c(2,1)*Qm  + c(2,2)*Q  +  c(2,3)*Qp ;
        Q2 = c(3,1)*Qmm + c(3,2)*Qm +  c(3,3)*Q  ;
        %
        if(ii<=data.quad.n)
            %
            Qq(:,i,j,ii) = Q0.*Omega0+ Q1.*Omega1 + Q2.*Omega2;
            %
        else
            %
            imid = data.weno.imid;
            %
            Qq(:,i,j,imid) = data.weno.sigma_p*Qq(:,i,j,imid) ...
                           - data.weno.sigma_m*(Q0.*Omega0+ Q1.*Omega1 + Q2.*Omega2);
            %
        end
        %
    end
    %
end