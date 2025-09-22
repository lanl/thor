function [J1, J2, r] = DecomposeJ(knot1n, knot2n, ctpxn, ctpyn, tt_tol)
    %Gauss data
    [Gp1, ~, ~] = getIntegralData1D(knot1n);
    [Gp2, ~, ~] = getIntegralData1D(knot2n);

    nn1 = length(Gp1);
    nn2 = length(Gp2);

    J = zeros(nn1,nn2);

    for i = 1:nn1
        for j = 1:nn2
            [~, DBS_Dz1, DBS_Dz2, Inon1, Inon2] = GetBspline2Dn(knot1n, knot2n, Gp1(i), Gp2(j));

            ctpx = ctpxn(Inon1, Inon2);
            ctpy = ctpyn(Inon1, Inon2);

            DX1_Dzeta1 = ctpx(:)' * DBS_Dz1(:);
            DX2_Dzeta1 = ctpy(:)' * DBS_Dz1(:);
   
            DX1_Dzeta2 = ctpx(:)' * DBS_Dz2(:);
            DX2_Dzeta2 = ctpy(:)' * DBS_Dz2(:);
    
            % Compute Jacobian determinant |J|
            J(i,j) = abs(DX1_Dzeta1 * DX2_Dzeta2 - DX2_Dzeta1 * DX1_Dzeta2);
        end
    end
    
    %% decompose omega
    Jtt = tt_tensor(J, tt_tol);
    J1 = Jtt{1};
    J2 = Jtt{2};
    r = Jtt.r(2);
end
