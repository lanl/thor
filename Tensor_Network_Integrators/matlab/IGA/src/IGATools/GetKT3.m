function K = GetKT3(knot1n, knot2n, zeta1, zeta2, J2, ctpxn, ctpyn, NWn, alpha_x, alpha_y)

    [n1, ~] = EvaluateKnot(knot1n);
    [n2, ~] = EvaluateKnot(knot2n);
    
    Inon1 = ActiveBSplineIndices(knot1n, zeta1);
    Inon2 = ActiveBSplineIndices(knot2n, zeta2);
    
    [DND_zeta1, DND_zeta2] = GetNURBS2D1DfTT(knot1n, knot2n, zeta1, zeta2, Inon1, Inon2, NWn);

    DND_zeta1_non = DND_zeta1(Inon1, Inon2);
    DND_zeta2_non = DND_zeta2(Inon1, Inon2);

    ctpx = ctpxn(Inon1, Inon2);
    ctpy = ctpyn(Inon1, Inon2);

    DX1D_zeta1 = ctpx(:)'*DND_zeta1_non(:);
    DX2D_zeta1 = ctpy(:)'*DND_zeta1_non(:);
   
    DX1D_zeta2 = ctpx(:)'*DND_zeta2_non(:);
    DX2D_zeta2 = ctpy(:)'*DND_zeta2_non(:);

    % Compute Jacobian determinant
    J1 = abs(DX1D_zeta1 * DX2D_zeta2 - DX2D_zeta1 * DX1D_zeta2);
    % Initialize stiffness matrices
    K1 = zeros(n1, n2, n1, n2);
    K2 = zeros(n1, n2, n1, n2);
    K3 = zeros(n1, n2, n1, n2);
    K4 = zeros(n1, n2, n1, n2);

    for i = 1:n1
        for j = 1:n2
            for k = 1:n1
                for l = 1:n2
                    K1(i,j,k,l) = (1/J1) * (DX2D_zeta2^2 + DX1D_zeta2^2) ...
                        * DND_zeta1(i,j) * DND_zeta1(k,l);

                    K2(i,j,k,l) = -(1/J1) * (DX2D_zeta1 * DX2D_zeta2 + DX1D_zeta2 * DX1D_zeta1) ...
                        * DND_zeta2(i,j) * DND_zeta1(k,l);

                    K3(i,j,k,l) = -(1/J1) * (DX2D_zeta1 * DX2D_zeta2 + DX1D_zeta2 * DX1D_zeta1) ...
                        * DND_zeta1(i,j) * DND_zeta2(k,l);

                    K4(i,j,k,l) = (1/J1) * (DX2D_zeta1^2 + DX1D_zeta1^2) ...
                        * DND_zeta2(i,j) * DND_zeta2(k,l);
                end
            end
        end
    end
    % Compute final stiffness matrix
    K = alpha_x*alpha_y*J2*(K1 + K2 + K3 + K4);
    
end
