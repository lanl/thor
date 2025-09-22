function Ftt = Getftt_3D(knot1n, knot2n, knot3n, ctpxn, ctpyn, ctpzn, ffunction, tt_tol)
ftt = Decomposef_3D(knot1n, knot2n, knot3n, ctpxn, ctpyn, ctpzn, ffunction, tt_tol);
R1 = ftt.r(2);
R2 = ftt.r(3);
fprintf('f: rank = %d %d\n', R1, R2);
for r1 = 1:R1
    for r2 = 1:R2
        %get vectors
        ft1 = ftt{1}(1,:,r1);
        ft2 = ftt{2}(r1,:,r2);
        ft3 = ftt{3}(r2,:,1);
        f1T = CalN_1D(knot1n, ft1);
        f2T = CalN_1D(knot2n, ft2);
        f3T = CalN_1D(knot3n, ft3);
        if r1*r2 == 1
            Ftt = matrices_to_tt_matrix_fn({f1T', f2T', f3T'});
        else
            Ftt = round(Ftt + matrices_to_tt_matrix_fn({f1T', f2T', f3T'}), tt_tol);
        end
    end
end

Ftt = round(Ftt, tt_tol);

end