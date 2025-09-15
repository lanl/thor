function data=recon(data,d)
    %
    weps = data.weno.weps;
    % get size
    Nx_total = data.Nx_total;
    Ny_total = data.Ny_total;
    Nz_total = data.Nz_total;
    %
    i = [1+2*(1==d):Nx_total-2*(1==d)];
    j = [1+2*(2==d):Ny_total-2*(2==d)];
    k = [1+2*(3==d):Nz_total-2*(3==d)];
    %
    % calculate shift indices
    imm = i - 2*(1==d); im = i - (1==d); ip = i + (1==d); ipp = i + 2*(1==d);
    jmm = j - 2*(2==d); jm = j - (2==d); jp = j + (2==d); jpp = j + 2*(2==d);
    kmm = k - 2*(3==d); km = k - (3==d); kp = k + (3==d); kpp = k + 2*(3==d);

    % reconstruct variables at the right side 
    Fmm = data.invFR(:,imm,jmm,kmm);
    Fm  = data.invFR(:,im,jm,km);
    F   = data.invFR(:,i,j,k);
    Fp  = data.invFR(:,ip,jp,kp);
    Fpp = data.invFR(:,ipp,jpp,kpp);
    %
    FR0 = ( 2*F   + 5*Fp -    Fpp)/6;
    FR1 = ( -Fm   + 5*F  +  2*Fp )/6;
    FR2 = ( 2*Fmm - 7*Fm + 11*F  )/6;
    %
    % compute smoothness parameters
    beta0 = (13/12*( F  - 2*Fp + Fpp).^2 + 1/4*(3*F - 4*Fp + Fpp).^2 + weps);
    %
    beta1 = (13/12*(Fm  -  2*F +  Fp).^2 + 1/4*(Fm - Fp).^2 + weps);
    %
    beta2 = (13/12*(Fmm - 2*Fm +   F).^2 + 1/4*(Fmm - 4*Fm + 3*F).^2 + weps);
    %
    params = paramsWeno5(beta0,beta1,beta2);
    %
    c0 = 0.3; c1 = 0.6; c2 = 0.1;
    %
    [Omega0R,Omega1R,Omega2R] = calcOmega(params,c0,c1,c2);
    %
    % reconstruct variables at the left side 
    Fmm = data.invFL(:,imm,jmm,kmm);
    Fm  = data.invFL(:,im,jm,km);
    F   = data.invFL(:,i,j,k);
    Fp  = data.invFL(:,ip,jp,kp);
    Fpp = data.invFL(:,ipp,jpp,kpp);
    %
    FL0 = ( 11*F   - 7*Fp + 2*Fpp)/6;
    FL1 = (  2*Fm  + 5*F  -   Fp )/6;
    FL2 = (   -Fmm + 5*Fm + 2*F  )/6;
    %
    % compute smoothness parameter
    beta0 = (13/12*( F  - 2*Fp + Fpp).^2 + 1/4*(3*F - 4*Fp + Fpp).^2 + weps);
    %
    beta1 = (13/12*(Fm  -  2*F +  Fp).^2 + 1/4*(Fm - Fp).^2 + weps);
    %
    beta2 = (13/12*(Fmm - 2*Fm +   F).^2 + 1/4*(Fmm - 4*Fm + 3*F).^2 + weps);
    %
    params = paramsWeno5(beta0,beta1,beta2);
    %
    c0 = 0.1; c1 = 0.6; c2 = 0.3;
    %
    [Omega0L,Omega1L,Omega2L] = calcOmega(params,c0,c1,c2);
    % 
    % calculate reconstructed variables
    %
    data.FR(:,i,j,k) = FR0.*Omega0R + FR1.*Omega1R + FR2.*Omega2R;
    %
    data.FL(:,i,j,k) = FL0.*Omega0L + FL1.*Omega1L + FL2.*Omega2L;  
    %
end