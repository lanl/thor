function K = GetKT2(knot1n, knot2n, zeta1, zeta2, J2, ctpxn, ctpyn, NWn, alpha_x, alpha_y)

    [n1, p1] = EvaluateKnot(knot1n);
    [n2, p2] = EvaluateKnot(knot2n);
    
    Inon1 = ActiveBSplineIndices(knot1n, p1, zeta1);
    Inon2 = ActiveBSplineIndices(knot2n, p2, zeta2);
    
    [DND_zeta11, DND_zeta12, DND_zeta21, DND_zeta22] = ...
        GetNURBS2D1DfTT(knot1n, knot2n, zeta1, zeta2, Inon1, Inon2, NWn);
    
    % Compute derivatives of physical coordinates w.r.t parametric space
    DND_zeta1 = DND_zeta11(Inon1, Inon2)+DND_zeta12(Inon1, Inon2);
    DND_zeta2 = DND_zeta21(Inon1, Inon2)+DND_zeta22(Inon1, Inon2);

    DXD_zeta1 = [sum(ctpxn(Inon1, Inon2).*DND_zeta1,'all'),sum(ctpyn(Inon1, Inon2).*DND_zeta1,'all')];
    DXD_zeta2  = [sum(ctpxn(Inon1, Inon2).*DND_zeta2,'all'),sum(ctpyn(Inon1, Inon2).*DND_zeta2,'all')];
    
    % Compute Jacobian determinant
    J1 = abs(DXD_zeta1(1) * DXD_zeta2(2) - DXD_zeta1(2) * DXD_zeta2(1));
    
    % Initialize stiffness matrices
    K1 = zeros(n1, n2, n1, n2);
    K2 = zeros(n1, n2, n1, n2);
    K3 = zeros(n1, n2, n1, n2);
    K4 = zeros(n1, n2, n1, n2);
    
    for i = 1:n1
        for j = 1:n2
            for k = 1:n1
                for l = 1:n2
                    K1(i,j,k,l) = (1/J1) * (DXD_zeta2(2)^2 + DXD_zeta2(1)^2) ...
                        * (DND_zeta11(i,j) + DND_zeta12(i,j)) ...
                        * (DND_zeta11(k,l) + DND_zeta12(k,l));

                    K2(i,j,k,l) = -(1/J1) * (DXD_zeta1(2) * DXD_zeta2(2) + DXD_zeta2(1) * DXD_zeta1(1)) ...
                        * (DND_zeta21(i,j) + DND_zeta22(i,j)) ...
                        * (DND_zeta11(k,l) + DND_zeta12(k,l));

                    K3(i,j,k,l) = -(1/J1) * (DXD_zeta1(2) * DXD_zeta2(2) + DXD_zeta2(1) * DXD_zeta1(1)) ...
                        * (DND_zeta11(i,j) + DND_zeta12(i,j)) ...
                        * (DND_zeta21(k,l) + DND_zeta22(k,l));

                    K4(i,j,k,l) = (1/J1) * (DXD_zeta1(2)^2 + DXD_zeta1(1)^2) ...
                        * (DND_zeta21(i,j) + DND_zeta22(i,j)) ...
                        * (DND_zeta21(k,l) + DND_zeta22(k,l));
                end
            end
        end
    end
    
    % Compute final stiffness matrix
    K = alpha_x * alpha_y * J2 * (K1 + K2 + K3 + K4);
    
end
