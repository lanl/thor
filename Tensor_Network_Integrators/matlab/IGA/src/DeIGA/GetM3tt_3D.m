function MT = GetM3tt_3D(knot1n, knot2n, knot3n, ctpxn, ctpyn, ctpzn,tt_tol)
I3 = [1,0,0;0,1,0;0,0,1]; 
Jtt = DecomposeJ_3D(knot1n, knot2n, knot3n, ctpxn, ctpyn, ctpzn, tt_tol);
R1 = Jtt.r(2);
R2 = Jtt.r(3);
for r1 = 1:R1
    for r2 = 1:R2
        Jt1 = Jtt{1}(1,:,r1);
        Jt2 = Jtt{2}(r1,:,r2);
        Jt3 = Jtt{3}(r2,:,1);

        M1T = CalNN_1D(knot1n, Jt1);
        M2T = CalNN_1D(knot2n, Jt2);
        M3T = CalNN_1D(knot3n, Jt3);
        if r1 == 1 && r2 == 1
            MT = round(matrices_to_tt_matrix_fn({I3, M1T, M2T, M3T}), tt_tol);
        else
            MT = round(MT + matrices_to_tt_matrix_fn({I3, M1T, M2T, M3T}), tt_tol);
        end
    end
end

end