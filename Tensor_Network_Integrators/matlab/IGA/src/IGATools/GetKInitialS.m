function KT = GetKInitialS(knot1n, knot2n, ctpxn, ctpyn, NWn, m2DT)
    % Evaluate knot vectors
    [n1, p1] = EvaluateKnot(knot1n);
    [n2, p2] = EvaluateKnot(knot2n);

    % Get Gauss quadrature points and weights
    [Gp1, alpha1] = GaussPoint(p1);
    [Gp2, alpha2] = GaussPoint(p2);

    % Extract subdomain boundaries from m2DT
    zeta11 = m2DT(:, :, 1);
    zeta12 = m2DT(:, :, 2);
    zeta21 = m2DT(:, :, 3);
    zeta22 = m2DT(:, :, 4);

    % Compute J2 scaling factor
    J2 = ((zeta12 - zeta11) .* (zeta22 - zeta21)) / 4;

    % Initialize stiffness matrix
    KT = zeros(n1, n2, n1, n2); 

    % Loop over all Gauss quadrature points and elements
    for kk = 1:length(Gp1)
        for ll = 1:length(Gp2)
            for ii = 1:(n1 - p1)
                for jj = 1:(n2 - p2)
                    % Compute quadrature points in the parametric domain
                    zeta1 = (zeta11(ii, jj) + zeta12(ii, jj) + Gp1(kk) * (zeta12(ii, jj) - zeta11(ii, jj))) / 2;
                    zeta2 = (zeta21(ii, jj) + zeta22(ii, jj) + Gp2(ll) * (zeta22(ii, jj) - zeta21(ii, jj))) / 2;

                    % Compute element stiffness matrix contribution
                    KT = KT + GetKT(knot1n, knot2n, zeta1, zeta2, J2(ii, jj), ctpxn, ctpyn, NWn, alpha1(kk), alpha2(ll));
                end
            end
        end
    end
end
