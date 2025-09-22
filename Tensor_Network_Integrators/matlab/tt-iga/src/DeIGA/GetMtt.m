function Mtt = GetMtt(knot1n, knot2n, ctpxn, ctpyn, tt_tol)
[Jt1, Jt2, R] = DecomposeJ(knot1n, knot2n, ctpxn, ctpyn, tt_tol);
for r = 1:R
    M1T = CalNN_1D(knot1n, Jt1(:,:,r));
    M2T = CalNN_1D(knot2n, Jt2(r,:,:));
    % M = M + kron(M2T, M1T); %full grid
    if r == 1
        Mtt = matrices_to_tt_matrix_fn({M1T,M2T});
    else
        Mtt = Mtt + matrices_to_tt_matrix_fn({M1T,M2T});
    end
end

end