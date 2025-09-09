function data=faceFlux(data,d)
    % compute + and - numerical fluxes for finite difference scheme
    Qstencil = getWeno5StencilNeq(data.tt.Q,d);
    %
    eps_cr = data.tt.eps_cr;
    %
    if(isempty(data.tt.F0L{d}) || isempty(data.tt.F0R{d}))
        rmax = 0;
        %
        for i=1:data.tt.Neq
            %
            r = max(data.tt.Q{i}.r);
            %
            if(r>rmax)
                idx  = i;
                rmax = r;
            end
            %
        end
        %
        Flux = cross_interpolation(Qstencil, @(x) computeFaceFlux(x,d,data), eps_cr, data.tt.Q{idx}, 0);
        %
    else
        Flux = cross_interpolation(Qstencil, @(x) computeFaceFlux(x,d,data), eps_cr, data.tt.F0L{d}, 0);
    end
    %
    Flux = q_unfold(Flux,data.tt.eps);
    %
    for i=1:data.tt.Neq
        %
        data.tt.FL{i} = Flux{i};
        data.tt.FR{i} = Flux{i+data.tt.Neq}; 
        %
    end
    %
    if(data.rki==1)
        rmaxL = 0;
        %
        for i=1:data.tt.Neq
            %
            rL = max(Flux{i}.r);
            %
            if(rL>rmaxL)
                idxL  = i;
                rmaxL = rL;
            end
            %
        end
        %
        data.tt.F0L{d} = Flux{idxL};
        %
    end
    %
end
%
function Fhat = computeFaceFlux(Q,d,data)
    %
    N   = size(Q,1);
    Neq = data.tt.Neq;
    %
    Q = reshape(Q,[N,Neq,5]);
    %
    % calculate cell centered fluxes and eigenvalues
    %
    Fhat  = zeros(N,2*Neq);
    FluxL = zeros(N,Neq,5);
    FluxR = zeros(N,Neq,5);
    EigQ  = zeros(N,Neq,5);
    %
    Flux(:,:,1) = data.F(Q(:,:,1),d,data.gam);
    Flux(:,:,2) = data.F(Q(:,:,2),d,data.gam);
    Flux(:,:,3) = data.F(Q(:,:,3),d,data.gam);
    Flux(:,:,4) = data.F(Q(:,:,4),d,data.gam);
    Flux(:,:,5) = data.F(Q(:,:,5),d,data.gam);
    %
    EigQ(:,:,1) = data.Eig(Q(:,:,1),d,data.gam);
    EigQ(:,:,2) = data.Eig(Q(:,:,2),d,data.gam);
    EigQ(:,:,3) = data.Eig(Q(:,:,3),d,data.gam);
    EigQ(:,:,4) = data.Eig(Q(:,:,4),d,data.gam);
    EigQ(:,:,5) = data.Eig(Q(:,:,5),d,data.gam);
    %
    EigQ = Q*max(EigQ(:));
    %
    FluxL(:)  = 0.5*(Flux(:) - EigQ(:));
    FluxR(:)  = 0.5*(Flux(:) + EigQ(:));
    %
    % Tranform to the characteristic field if needed
    %
    if(data.ReconChar)
        %
        FhatCL = zeros(N,Neq);
        FhatCR = zeros(N,Neq);
        %
        FluxCL = zeros(size(FluxL));
        FluxCR = zeros(size(FluxR));
        %
        temp  = zeros(N,1);
        %
        QhL = 0.5*(Q(:,:,3)+Q(:,:,2));
        QhR = 0.5*(Q(:,:,3)+Q(:,:,4));
        %
        [RhL,LhL]=eigVecs(QhL,d,data.gam);
        [RhR,LhR]=eigVecs(QhR,d,data.gam);
        %
        % Left
        %
        for k=1:5
            %
            for i=1:Neq
                %
                temp(:) = 0;
                %
                for j=1:Neq
                    %
                    temp(:) = temp(:) + LhL(:,i,j).*FluxL(:,j,k);
                    %
                end
                %
                FluxCL(:,i,k) = temp;
                %
            end
        end
        %
        % Right
        %
        for k=1:5
            %
            for i=1:Neq
                %
                temp(:) = 0;
                %
                for j=1:Neq
                    %
                    temp(:) = temp(:) + LhR(:,i,j).*FluxR(:,j,k);
                    %
                end
                %
                FluxCR(:,i,k) = temp;
                %
            end
        end
        %
    end
    %
    weps = data.weno.weps;
    %
    % Tranform back to the physical field
    %
    if(data.ReconChar)
        %
        FhatCL(:,:) = fun_teno5(FluxCL,-1,weps);
        FhatCR(:,:) = fun_teno5(FluxCR,+1,weps);
        % 
        % Left
        %
        for i=1:Neq
            %
            temp(:) = 0;
            %
            for j=1:Neq
                %
                temp(:) = temp(:) + RhL(:,i,j).*FhatCL(:,j);
                %
            end
            %
            Fhat(:,i) = temp;
            %
        end
        %
        % Right
        %
        for i=1:Neq
            %
            temp(:) = 0;
            %
            for j=1:Neq
                %
                temp(:) = temp(:) + RhR(:,i,j).*FhatCR(:,j);
                %
            end
            %
            Fhat(:,i+Neq) = temp;
            %
        end
        %
    else
        %
        Fhat(:,1:Neq)       = fun_teno5(FluxL,-1,weps);
        Fhat(:,Neq+1:2*Neq) = fun_teno5(FluxR,+1,weps);
        %
    end
    %
end