function data=localLF(data,d)
    %
    [~,idx] = getMaxRank(data.tt.Q);
    %
    Flux = cross_interpolation(data.tt.Q, @(x) fun_localLF(x,d,data.F,data.Eig,data.gam), data.tt.eps_cr, data.tt.Q{idx}, 0);
    %
    Flux = q_unfold(Flux,data.tt.eps);
    %
    for i=1:data.tt.Neq
        %
        j = i + data.tt.Neq;
        %
        data.tt.invFL{i} = Flux{i};
        data.tt.invFR{i} = Flux{j}; 
        %
    end
    %
end

function Flux = fun_localLF(Q,d,F,Eig,gam)
    %
    n = size(Q,2);
    %
    Flux = zeros(size(Q,1),2*n);
    %
    % Calculate Flux and Eigenvalues
    %
    F_inv  = F(Q,d,gam);
    LQ     = Eig(Q,d,gam);
    LQ(:)  = Q(:).*LQ(:);
    %
    Flux(:,1:n)     = 0.5*(F_inv - LQ); % left
    Flux(:,n+1:2*n) = 0.5*(F_inv + LQ); % right
    %
end