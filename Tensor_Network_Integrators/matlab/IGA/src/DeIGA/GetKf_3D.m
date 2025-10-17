function Kf = GetKf_3D(knot1n, knot2n, knot3n, ctpxn, ctpyn, ctpzn)
    %Gauss data
    [Gp1, ~, ~] = getIntegralData1D(knot1n);
    [Gp2, ~, ~] = getIntegralData1D(knot2n);
    [Gp3, ~, ~] = getIntegralData1D(knot3n);

    nn1 = length(Gp1);
    nn2 = length(Gp2);
    nn3 = length(Gp3);

   
    Kf = zeros(nn1, nn2, nn3, 3, 3);
    for i1 = 1:nn1
        for i2 = 1:nn2
            for i3 =1:nn3
                [~, DBS_Dz1, DBS_Dz2, DBS_Dz3, Inon1, Inon2, Inon3] = GetBspline3Dn(knot1n, knot2n, knot3n, ...
                                                                  Gp1(i1), Gp2(i2), Gp3(i3));
            
                ctpx = ctpxn(Inon1, Inon2, Inon3);
                ctpy = ctpyn(Inon1, Inon2, Inon3);
                ctpz = ctpzn(Inon1, Inon2, Inon3);
                %%
                DX1_Dzeta1 = ctpx(:)' * DBS_Dz1(:);
                DX2_Dzeta1 = ctpy(:)' * DBS_Dz1(:);
                DX3_Dzeta1 = ctpz(:)' * DBS_Dz1(:);

                DX1_Dzeta2 = ctpx(:)' * DBS_Dz2(:);
                DX2_Dzeta2 = ctpy(:)' * DBS_Dz2(:);
                DX3_Dzeta2 = ctpz(:)' * DBS_Dz2(:);

                DX1_Dzeta3 = ctpx(:)' * DBS_Dz3(:);
                DX2_Dzeta3 = ctpy(:)' * DBS_Dz3(:);
                DX3_Dzeta3 = ctpz(:)' * DBS_Dz3(:);
                %%
                DX1Dz = [DX1_Dzeta1, DX1_Dzeta2, DX1_Dzeta3];
                DX2Dz = [DX2_Dzeta1, DX2_Dzeta2, DX2_Dzeta3];
                DX3Dz = [DX3_Dzeta1, DX3_Dzeta2, DX3_Dzeta3];

                %% Compute Jacobian determinant |J|
                j = [DX1Dz; DX2Dz; DX3Dz];
                % Inverse of Jacobian
                ij = inv(j);
                % Jacobian determinant |J|
                J = abs(det(j));
                G = J*(ij*ij');
                % K factor
                Kf(i1, i2, i3, :, :) = [G(1,1) G(1,2) G(1,3);...
                                        G(2,1) G(2,2) G(2,3);... 
                                        G(3,1) G(3,2) G(3,3)];
            end
        end
    end

end
