function [W] = partial_contraction_rl(X, Y)
    % Algorithm 2.3 in Daas et al. (2023):
    % Right-to-left contraction of tensors X and Y.
    % Require: Tensors X, Y with consistent dimensions in TT format and ranks
    %          rx(:) and ry(:), respectively.
    % Ensure:  Matrices {W(:)} satisfy W(n-1) = H[X_n] @ H[Y_n]^T for 1 < n <= N
    %   function [W] = PartialContractionsRL(X, Y)
    %     W_{N-1} = H[X_N] @ H[Y_N]^T 
    %     for n = N-1 down to 2 do
    %       V[Z_n] = V[X_n] @ W_n 
    %       W_{n-1} = H[Z_n] @ H[Y_n]^T
    %     end for
    %   end function
    d = X.d;
    W = cell(d,1);
    HX = hstack(X, d);
    HY = hstack(Y, d);
    W{d-1} = HX*HY'; % (X.r(d), X.n(d)*X.r(d+1)) * (Y.r(d), Y.n(d)*Y.r(d+1))'
    for k = d-1:-1:2
        VX = vstack(X, k);
        VZ = VX*W{k};
        Zk = reshape(VZ,[X.r(k),X.n(k),size(W{k},2)]);
        HZ = reshape(permute(Zk, [1,3,2]),[X.r(k),X.n(k)*size(W{k},2)]);
        HY = hstack(Y, k);
        W{k-1} = HZ*HY';
    end
    return;
end

