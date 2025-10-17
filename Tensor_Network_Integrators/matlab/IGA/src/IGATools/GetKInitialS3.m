function KT = GetKInitialS3(knot1n, knot2n, ctpxn, ctpyn, NWn, m2DT)
    % Evaluate knots
    [n1, p1] = EvaluateKnot(knot1n);
    [n2, p2] = EvaluateKnot(knot2n);

    % Get Gaussian points
    [Gp1, Alpha1] = GaussPoint(p1);
    [Gp2, Alpha2] = GaussPoint(p2);
    
    % Compute Gpv1 and Gpv2
    [Gpv1, Gpv2] = GetGpv(Gp1, Gp2);
    
    % Compute AlphaV1 and AlphaV2
    numElements = (n1 - p1) * (n2 - p2);
    AlphaV1 = repmat([Alpha1(1); Alpha1(1); Alpha1(2); Alpha1(2)], numElements, 1);
    AlphaV2 = repmat([Alpha2(1); Alpha2(2); Alpha2(1); Alpha2(2)], numElements, 1);
    
    % Extract transformation matrices from m2DT
    Zeta11 = m2DT(:, :, 1);
    Zeta12 = m2DT(:, :, 2);
    Zeta21 = m2DT(:, :, 3);
    Zeta22 = m2DT(:, :, 4);
    
    numGP = p1 * p2;
    numPoints = numElements * numGP;
    
    ZetaG1 = zeros(numPoints, 1);
    ZetaG2 = zeros(numPoints, 1);
    J2 = zeros(numPoints, 1);

    count = 1;
    for ii = 1:(n1 - p1)
        for jj = 1:(n2 - p2)
            for kk = 1:numGP
                ZetaG1(count) = 0.5*(Zeta11(ii, jj) + Zeta12(ii, jj) + Gpv1(kk)*(Zeta12(ii, jj) - Zeta11(ii,jj)));
                ZetaG2(count) = 0.5*(Zeta21(ii, jj) + Zeta22(ii, jj) + Gpv2(kk)*(Zeta22(ii, jj) - Zeta21(ii,jj)));
                J2(count) = 0.25*(Zeta12(ii, jj) - Zeta11(ii, jj))*(Zeta22(ii, jj) - Zeta21(ii, jj));
                count = count + 1;
            end
        end
    end
       
    KT = 0;
    for ii = 1:numPoints
        KT = KT + GetKT3(knot1n, knot2n, ZetaG1(ii), ZetaG2(ii), J2(ii), ctpxn, ctpyn, NWn, AlphaV1(ii), AlphaV2(ii));
    end
end
