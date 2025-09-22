function Kf = GetKf_2D(knot1n, knot2n, ctpxn, ctpyn)
    %Gauss data
    [Gp1, ~, ~] = getIntegralData1D(knot1n);
    [Gp2, ~, ~] = getIntegralData1D(knot2n);

    nn1 = length(Gp1);
    nn2 = length(Gp2);

    Kf = zeros(nn1, nn2, 2, 2);
    for i1 = 1:nn1
        for i2 = 1:nn2
            [~, DBS_Dz1, DBS_Dz2, Inon1, Inon2] = GetBspline2Dn(knot1n, knot2n, Gp1(i1), Gp2(i2));

            ctpx = ctpxn(Inon1, Inon2);
            ctpy = ctpyn(Inon1, Inon2);

            DX1_Dzeta1 = ctpx(:)' * DBS_Dz1(:);
            DX2_Dzeta1 = ctpy(:)' * DBS_Dz1(:);
   
            DX1_Dzeta2 = ctpx(:)' * DBS_Dz2(:);
            DX2_Dzeta2 = ctpy(:)' * DBS_Dz2(:);
    
            DX1Dz = [DX1_Dzeta1, DX1_Dzeta2];
            DX2Dz = [DX2_Dzeta1, DX2_Dzeta2];
            % Jacobian
            j = [DX1Dz; DX2Dz];
            % Inverse of Jacobian
            ij = inv(j);
            % Jacobian determinant |J|
            J = abs(det(j));
            G = J*(ij*ij');
            %J = abs(DX1_Dzeta1 * DX2_Dzeta2 - DX2_Dzeta1 * DX1_Dzeta2);
            % K factor
            %K11f(i1,i2) = (1/J) * (DX2_Dzeta2^2 + DX1_Dzeta2^2);
            %K12f(i1,i2) = -(1/J) * (DX2_Dzeta1 * DX2_Dzeta2 + DX1_Dzeta2 * DX1_Dzeta1);
            %K21f(i1,i2) = -(1/J) * (DX2_Dzeta1 * DX2_Dzeta2 + DX1_Dzeta2 * DX1_Dzeta1);
            %K22f(i1,i2) = (1/J) * (DX2_Dzeta1^2 + DX1_Dzeta1^2);

            Kf(i1, i2, :, :) = [G(1,1) G(1,2); G(2,1) G(2,2)];

        end
    end

end
