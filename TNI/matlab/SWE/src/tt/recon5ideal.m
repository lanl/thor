function data=recon5ideal(data,d)
    %
    % computes a third order reconstruction at the quadrature points over the cell interfaces in the TT format
    %
    % get size
    Neq  = data.tt.Neq;
    n    = [data.Nx_total;data.Ny_total];
    %
    % calculate shift indices
    %
    j = [3:n(d)-2];
    %
    jmm = j - 2;
    jm  = j - 1;
    jp  = j + 1;
    jpp = j + 2;   
    %
    % loop over equations
    %
    data.tt.QL0 = data.tt.Q;
    data.tt.QR0 = data.tt.Q;
    %
    for i=1:Neq
        %
        Qvar = data.tt.Q{i};
        %
        r    = Qvar.r;
        core = reshape(Qvar{d}(:),[r(d),n(d),r(d+1)]);
        %
        Qmm = core(:,jmm,:);
        Qm  = core(:,jm ,:);
        Q   = core(:,j  ,:);
        Qp  = core(:,jp ,:);
        Qpp = core(:,jpp,:);
        %
        QL = core;
        QR = core;
        %
        % calculate weno weights on the right
        %
        Omega0R = 0.3; Omega1R = 0.6; Omega2R = 0.1;
        %
        % reconstruct variables at the right side 
        QR0 = ( 2*Q   + 5*Qp -    Qpp)/6;
        QR1 = ( -Qm   + 5*Q  +  2*Qp )/6;
        QR2 = ( 2*Qmm - 7*Qm + 11*Q  )/6;
        %
        % calculate weno weights on the left
        %
        Omega0L = 0.1; Omega1L = 0.6; Omega2L = 0.3;
        %
        % reconstruct variables at the left side 
        QL0 = ( 11*Q   - 7*Qp + 2*Qpp)/6;
        QL1 = (  2*Qm  + 5*Q  -   Qp )/6;
        QL2 = (   -Qmm + 5*Qm + 2*Q  )/6;
        %
        % calculate the reconstructed variables
        %
        QR(:,j,:) = QR0.*Omega0R + QR1.*Omega1R + QR2.*Omega2R;
        QL(:,j,:) = QL0.*Omega0L + QL1.*Omega1L + QL2.*Omega2L;
        % 
        data.tt.QL0{i}{d} = QL;
        data.tt.QR0{i}{d} = QR;
        % 
        ttR = quad_tt(n,r,data.quad.n,d);
        ttL = quad_tt(n,r,data.quad.n,d);
        %
        for jj=1:data.tt.d
            %
            if(jj==d)
                ttR{jj} = reshape(QR,size(ttR{jj}));
                ttL{jj} = reshape(QL,size(ttL{jj}));
            else
                ttR{jj} = reshape(data.tt.Qq{i}{jj},size(ttR{jj}));
                ttL{jj} = reshape(data.tt.Qq{i}{jj},size(ttL{jj}));
            end
        end
        %
        data.tt.QR{i} = ttR;
        data.tt.QL{i} = ttL;
        %
    end
    %
end
%
function tt=quad_tt(N,r,nq,dir)
    %
    tt     = tt_tensor;
    %
    tt.n      = N*nq;
    tt.n(dir) = N(dir);
    tt.d      = 2;
    tt.r      = r;
    tt.ps     = cumsum([1;tt.n.*tt.r(1:tt.d).*tt.r(2:tt.d+1)]);
    tt.core   = zeros(tt.ps(tt.d+1)-1,1);
    %
end