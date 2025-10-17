function M = GetM3(knot1n, knot2n, zeta1, zeta2, J2, ctpxn, ctpyn, alpha_x, alpha_y)

    [n1, ~] = EvaluateKnot(knot1n);
    [n2, ~] = EvaluateKnot(knot2n);
    
    [BS, DBS_Dz1, DBS_Dz2, Inon1, Inon2] = GetBspline2Df(knot1n, knot2n, zeta1, zeta2);

    %[DND_zeta1, DND_zeta2] = GetNURBS2D1DfTT(knot1n, knot2n, zeta1, zeta2, Inon1, Inon2, NWn);

    DND_zeta1_non = DBS_Dz1(Inon1, Inon2);
    DND_zeta2_non = DBS_Dz2(Inon1, Inon2);

    ctpx = ctpxn(Inon1, Inon2);
    ctpy = ctpyn(Inon1, Inon2);

    DX1D_zeta1 = ctpx(:)'*DND_zeta1_non(:);
    DX2D_zeta1 = ctpy(:)'*DND_zeta1_non(:);
   
    DX1D_zeta2 = ctpx(:)'*DND_zeta2_non(:);
    DX2D_zeta2 = ctpy(:)'*DND_zeta2_non(:);

    % Compute Jacobian determinant
    J1 = abs(DX1D_zeta1 * DX2D_zeta2 - DX2D_zeta1 * DX1D_zeta2);
    % Initialize stiffness matrices
    M = zeros(n1, n2, n1, n2);

    for i = 1:n1
        for j = 1:n2
            for k = 1:n1
                for l = 1:n2
                    M(i,j,k,l) = BS(i,j) * BS(k,l);
                end
            end
        end
    end
    % Compute final stiffness matrix
    M = alpha_x * alpha_y * J2 * J1 * (M);

end
