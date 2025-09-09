function data = recon5(data,d)
    %
    weps = data.weno.weps;
    % get size
    Nx_total = data.Nx_total;
    Ny_total = data.Ny_total;
    %
    i = [1+2*(1==d):Nx_total-2*(1==d)];
    j = [1+2*(2==d):Ny_total-2*(2==d)];
    %
    % calculate shift indices
    imm = i - 2*(1==d); im = i - (1==d); ip = i + (1==d); ipp = i + 2*(1==d);
    jmm = j - 2*(2==d); jm = j - (2==d); jp = j + (2==d); jpp = j + 2*(2==d);

    % conserved variable-based reconstruction
    Qmm = data.Q(:,imm,jmm);
    Qm  = data.Q(:,im,jm);
    Q   = data.Q(:,i,j);
    Qp  = data.Q(:,ip,jp);
    Qpp = data.Q(:,ipp,jpp);
    %
    params = data.reconParams(Qmm,Qm,Q,Qp,Qpp,weps);
    %
    % calculate weno weights on the right
    %
    c0 = 0.3; c1 = 0.6; c2 = 0.1;
    % reconstruct variables at the right side 
    QR0 = ( 2*Q   + 5*Qp -    Qpp)/6;
    QR1 = ( -Qm   + 5*Q  +  2*Qp )/6;
    QR2 = ( 2*Qmm - 7*Qm + 11*Q  )/6;
    %
    [Omega0R,Omega1R,Omega2R] = data.calcOmega(params,c0,c1,c2);
    %
    % calculate weno weights on the left
    %
    c0 = 0.1; c1 = 0.6; c2 = 0.3;
    % reconstruct variables at the left side 
    QL0 = ( 11*Q   - 7*Qp + 2*Qpp)/6;
    QL1 = (  2*Qm  + 5*Q  -   Qp )/6;
    QL2 = (   -Qmm + 5*Qm + 2*Q  )/6;
    %
    [Omega0L,Omega1L,Omega2L] = data.calcOmega(params,c0,c1,c2);
    %
    data.QR(:,i,j,1) = QR0.*Omega0R + QR1.*Omega1R + QR2.*Omega2R;
    %
    data.QL(:,i,j,1) = QL0.*Omega0L + QL1.*Omega1L + QL2.*Omega2L;  
    %
    data.QR(:,:,:,:) = recon5GQ(data.QR,data,d);
    data.QL(:,:,:,:) = recon5GQ(data.QL,data,d);
    %
end