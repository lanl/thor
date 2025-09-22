function K_final = GetKT(knot1n, knot2n, zeta1, zeta2, J2, ctpxn, ctpyn, NWn, alpha_x, alpha_y)

    % Compute active B-spline indices
    Inon1 = ActiveBSplineIndices(knot1n, zeta1);
    Inon2 = ActiveBSplineIndices(knot2n, zeta2);

    % Compute 2D NURBS basis functions and their derivatives
    [~, DND_zeta1, DND_zeta2] = GetNURBS2D1DfT(knot1n, knot2n, zeta1, zeta2, Inon1, Inon2, NWn);

    % Compute physical coordinates X
    %X = sum(NN(Inon1, Inon2) .* ctpn(Inon1, Inon2), 'all');
    
    % Compute derivatives of physical coordinates w.r.t parametric space
    DXD_zeta1 = [sum(ctpxn(Inon1, Inon2).*DND_zeta1(Inon1, Inon2),'all'),sum(ctpyn(Inon1, Inon2).*DND_zeta1(Inon1, Inon2),'all')];
    DXD_zeta2  = [sum(ctpxn(Inon1, Inon2).*DND_zeta2(Inon1, Inon2),'all'),sum(ctpyn(Inon1, Inon2).*DND_zeta2(Inon1, Inon2),'all')];
    
    % Compute Jacobian determinant
    J1 = abs(DXD_zeta1(1) * DXD_zeta2(2) - DXD_zeta1(2) * DXD_zeta2(1));

    % Compute inverse Jacobian matrix
    D_zetaDX = (1 / J1) * [ DXD_zeta2(2), -DXD_zeta2(1); 
                           -DXD_zeta1(2),  DXD_zeta1(1)];

    % Compute stiffness matrix
    K = CalKT(DND_zeta1, DND_zeta2, D_zetaDX, knot1n, knot2n);

    % Final scaled stiffness matrix
    K_final = alpha_x * alpha_y * (J1 * J2) * K;
end
