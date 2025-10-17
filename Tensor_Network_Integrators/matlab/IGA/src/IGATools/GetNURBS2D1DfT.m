function [N2D, DNURBS_zeta1, DNURBS_zeta2] = GetNURBS2D1DfT(knot1n, knot2n, zeta1, zeta2, Inon1, Inon2, NWn)
    % Evaluate knot vectors
    [n1, ~] = EvaluateKnot(knot1n);
    [n2, ~] = EvaluateKnot(knot2n);

    % Initialize arrays for basis functions and derivatives
    NonUBS1 = zeros(1, n1);
    DNonUBS1 = zeros(1, n1);
    NonUBS2 = zeros(1, n2);
    DNonUBS2 = zeros(1, n2);

    % Order of derivative
    m = 1;

    % Compute 1D basis functions and derivatives
    NN1= GetDNUBSf(knot1n,zeta1,m);
    NonUBS1(Inon1)  = NN1(:,1);
    DNonUBS1(Inon1) = NN1(:,2);
    
    NN2= GetDNUBSf(knot2n,zeta2,m);
    NonUBS2(Inon2)  = NN2(:,1);
    DNonUBS2(Inon2) = NN2(:,2);
        
    % Initialize 2D NURBS basis function matrices
    N2D = zeros(n1, n2);
    DNURBS_zeta1 = zeros(n1, n2);
    DNURBS_zeta2 = zeros(n1, n2);

    % Compute the weight function
    W = sum(NonUBS1(Inon1)' * NonUBS2(Inon2) .* NWn(Inon1, Inon2), 'all');
    
    % Compute weight derivatives
    DW_zeta1 = sum(DNonUBS1(Inon1)' * NonUBS2(Inon2) .* NWn(Inon1, Inon2), 'all');
    DW_zeta2 = sum(NonUBS1(Inon1)' * DNonUBS2(Inon2) .* NWn(Inon1, Inon2), 'all');

    % Compute the NURBS basis functions and their derivatives
    for i = 1:length(Inon1)
        for j = 1:length(Inon2)
            w_ij = NWn(Inon1(i), Inon2(j));
            phi_ij = NonUBS1(Inon1(i)) * NonUBS2(Inon2(j));

            N2D(Inon1(i), Inon2(j)) = (w_ij / W) * phi_ij;
            
            DNURBS_zeta1(Inon1(i), Inon2(j)) = (w_ij / W^2) * ...
                (DNonUBS1(Inon1(i)) * NonUBS2(Inon2(j)) * W - phi_ij * DW_zeta1);
            
            DNURBS_zeta2(Inon1(i), Inon2(j)) = (w_ij / W^2) * ...
                (NonUBS1(Inon1(i)) * DNonUBS2(Inon2(j)) * W - phi_ij * DW_zeta2);
        end
    end
end
