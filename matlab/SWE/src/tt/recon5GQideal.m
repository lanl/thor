function data=recon5GQideal(data)
    %
    % computes a fifth order reconstruction at the quadrature points over the cell volume in the TT format
    %
    weps = data.weno.weps;
    nq   = data.quad.n;
    % get size
    Neq  = data.tt.Neq;
    n    = [data.Nx_total;data.Ny_total];
    %
    % loop over equations
    %
    for i=1:Neq
        %
        Qvar = data.tt.Q{i};
        %
        r = Qvar.r;
        %
        for d=1:data.tt.d
            %
            data.Qcore{d} = zeros(r(d),nq,n(d),r(d+1));
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
            core = reshape(Qvar{d}(:),[r(d),n(d),r(d+1)]);
            %
            Qmm = core(:,jmm,:);
            Qm  = core(:,jm ,:);
            Q   = core(:,j  ,:);
            Qp  = core(:,jp ,:);
            Qpp = core(:,jpp,:);
            %
            for ii=1:size(data.weno.gamma,2)
                %
                Omega0 = data.weno.gamma(1,ii); 
                Omega1 = data.weno.gamma(2,ii); 
                Omega2 = data.weno.gamma(3,ii);
                %
                c = data.weno.C{ii};
                %
                Q0 = c(1,1)*Q   + c(1,2)*Qp +  c(1,3)*Qpp;
                Q1 = c(2,1)*Qm  + c(2,2)*Q  +  c(2,3)*Qp ;
                Q2 = c(3,1)*Qmm + c(3,2)*Qm +  c(3,3)*Q  ;
                %
                Qrecon = Q0.*Omega0+ Q1.*Omega1 + Q2.*Omega2;
                %
                if(ii<=data.quad.n)
                    data.Qcore{d}(:,ii,j,:) = reshape(Qrecon,size(data.Qcore{d}(:,ii,j,:)));
                else
                    imid = data.weno.imid;
                    data.Qcore{d}(:,imid,j,:) = data.weno.sigma_p*data.Qcore{d}(:,imid,j,:) ...
                                              - data.weno.sigma_m*reshape(Qrecon,size(data.Qcore{d}(:,imid,j,:)));
                end
            end
            %
            % BC treatment to prevent amen_cross convergence issues when conserved variables are equal to 0
            %
            for ii=1:nq
                data.Qcore{d}(:,ii,[1:2,n(d)-1,n(d)],:) = core(:,[1:2,n(d)-1,n(d)],:);
            end
        end
        %
        data.tt.Qq{i}=quad_tt(nq*n,r);
        %
        for d=1:data.tt.d
            data.tt.Qq{i}{d}=reshape(data.Qcore{d},size(data.tt.Qq{i}{d}));
        end
        %
    end
    %
    return;
end
%
function tt=quad_tt(n,r)
    %
    tt      = tt_tensor;
    tt.n    = n;
    tt.d    = 2;
    tt.r    = r;
    tt.ps   = cumsum([1;tt.n.*tt.r(1:tt.d).*tt.r(2:tt.d+1)]);
    tt.core = zeros(tt.ps(tt.d+1)-1,1);
    %
end