function [tt] = rounding_rand_orth(Y,R,rnew)
    % Algorithm 3.1 in Daas et al. (2023)
    % TT-rounding: randomize than orthogonolize
    %  - Y    : tensor to be rounded;
    %  - rnew : new rank (a scalar).
    % function
    %    Select a random Gaussian TT-tensor R with target TT-ranks rr(:) 
    %    W{:} = PartialContractionsRL(Y, R)
    %    X_1 = Y_1
    %    for n = 1 to N-1 do
    %      Z_n = V[X_n]           % V[X_n] is rr(n-1) x I(n) x r(n)
    %      P_n = Z_n * W{n}       % form the sketched matrix
    %      [V[X_n], ~] = QR(P_n)  % thin QR to compute an orthonormal basis
    %      M_n = V[X_n]^T @ Z_n   % form rr(n) x Rn matrix
    %      H(X_{n+1}) = M_n @ H[Y_{n+1}]
    %    end for
    % end function
    d = Y.d;
%    R = tt_rand(Y.n(1),d,rnew,1);
    rr= R.r;
    W = partial_contraction_rl(Y, R);
    X = cell(d, 1);
    X{1} = get_core(Y,1);
    for k = 1:d-1
        Zk = reshape(X{k},[size(X{k},1)*size(X{k},2),size(X{k},3)]);
        Pk = Zk * W{k};
        [VX, ~] = qr(Pk, "econ");
        X{k} = reshape(VX, [size(X{k},1),size(X{k},2),size(VX,2)]);
        Mk = VX' * Zk;
        HY = hstack(Y, k+1);
        HX = Mk  * HY;
        X{k+1} = permute(reshape(HX, [size(HX,1),Y.r(k+2),Y.n(k+1)]),[1,3,2]);
    end
    for i=2:d-1
        X{i}=permute(X{i},[2,1,3]);
    end
    X{1} = reshape(X{1}(:), [size(X{1},2),size(X{1},3)]);
    X{d} = permute(X{d},[2,1]);
    tt = tt_tensor(X);
end

