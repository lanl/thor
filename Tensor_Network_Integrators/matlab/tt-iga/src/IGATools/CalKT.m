function K = CalKT(DND_zeta1, DND_zeta2, D_zetaDX, knot1n, knot2n)
    % Evaluate knot vectors
    [n1, ~] = EvaluateKnot(knot1n);
    [n2, ~] = EvaluateKnot(knot2n);

    % Initialize B tensor (2 x n1 x n2)
    B = zeros(2, n1, n2);

    % Assign values to B tensor
    for i = 1:n1
        for j = 1:n2
            B(1, i, j) = DND_zeta1(i, j);
            B(2, i, j) = DND_zeta2(i, j);
        end
    end

    % Initialize stiffness matrix K (n1 x n2 x n1 x n2)
    K = zeros(n1, n2, n1, n2);

    % Compute stiffness matrix
    for i = 1:n1
        for j = 1:n2
            for k = 1:n1
                for l = 1:n2
                    K(i, j, k, l) = ...
                        (D_zetaDX(1, 1)^2 + D_zetaDX(1, 2)^2) * B(1, i, j) * B(1, k, l) ...
                      + (D_zetaDX(1, 1) * D_zetaDX(2, 1) + D_zetaDX(1, 2) * D_zetaDX(2, 2)) * B(2, i, j) * B(1, k, l) ...
                      + (D_zetaDX(1, 1) * D_zetaDX(2, 1) + D_zetaDX(1, 2) * D_zetaDX(2, 2)) * B(1, i, j) * B(2, k, l) ...
                      + (D_zetaDX(2, 1)^2 + D_zetaDX(2, 2)^2) * B(2, i, j) * B(2, k, l);
                end
            end
        end
    end
end
