function KT = GetKtt(knot1n, knot2n, ctpxn, ctpyn, tt_tol)
Kf = GetKf_2D(knot1n, knot2n, ctpxn, ctpyn);

[n1, ~] = EvaluateKnot(knot1n);
[n2, ~] = EvaluateKnot(knot2n);
sz1 = [n1, n2];
sz2 = [n1, n2]; 

KT = tt_matrix(zeros(prod(sz1), prod(sz2)), 1e-12, sz1, sz2);
Max_Kf = max(Kf(:,:,:,:));

for i =1:2
    for j =1:2
        if max(abs(Kf(:,:,i,j))) > 1e-8*Max_Kf
            Kf_tt = tt_tensor(Kf(:,:,i,j), tt_tol);
            fprintf('K_%d%d: rank = %d \n', i, j, Kf_tt.r(2));
            for r1 = 1:Kf_tt.r(2)
                    Kx = Kf_tt{1}(1,:,r1);
                    Ky = Kf_tt{2}(r1,:,1);
                    if i==1 && j==1
                        KxT = CalDNDN_1D(knot1n, Kx);
                        KyT = CalNN_1D(knot2n, Ky);
                    elseif i==1 && j==2
                        KxT = CalNDN_1D(knot1n, Kx);
                        KyT = CalDNN_1D(knot2n, Ky);
                    elseif i==2 && j==1
                        KxT = CalDNN_1D(knot1n, Kx);
                        KyT = CalNDN_1D(knot2n, Ky);
                    elseif i==2 && j==2
                        KxT = CalNN_1D(knot1n, Kx);
                        KyT = CalDNDN_1D(knot2n, Ky);
                    end
                    KT = KT + matrices_to_tt_matrix_fn({KxT, KyT});
            end
        end
    end
end

KT = round(KT, tt_tol);
end
