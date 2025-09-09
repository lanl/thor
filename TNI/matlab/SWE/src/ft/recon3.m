function data = recon3(data,d)
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
    im = i - (1==d); ip = i + (1==d);
    jm = j - (2==d); jp = j + (2==d);

    % conserved variable-based reconstruction
    Qm  = data.Q(:,im,jm);
    Q   = data.Q(:,i,j);
    Qp  = data.Q(:,ip,jp);
    %
    params = data.reconParams(Qm,Q,Qp,weps);
    %
    % calculate weno weights on the right
    %
    c0 = 2/3; c1 = 1/3; 
    % reconstruct variables at the right side 
    QR0 = (   Q +   Qp)/2;
    QR1 = ( -Qm + 3*Q )/2;
    %
    [Omega0R,Omega1R] = data.calcOmega(params,c0,c1);
    %
    % calculate weno weights on the left
    %
    c0 = 1/3; c1 = 2/3; 
    % reconstruct variables at the left side 
    QL0 = ( 3*Q - Qp)/2;
    QL1 = (  Qm + Q )/2;
    %
    [Omega0L,Omega1L] = data.calcOmega(params,c0,c1);
    %
    % calculate reconstructed variables
    %
    data.QR(:,i,j,1) = QR0.*Omega0R + QR1.*Omega1R;
    %
    data.QL(:,i,j,1) = QL0.*Omega0L + QL1.*Omega1L; 
    %
    data.QR(:,:,:,:) = recon3GQ(data.QR,data,d);
    data.QL(:,:,:,:) = recon3GQ(data.QL,data,d);
    %
end