function [BS, DBS_Dz1, DBS_Dz2, Inon1, Inon2] = GetBspline2Df(knot1n, knot2n, zeta1, zeta2)
    
    Inon1 = ActiveBSplineIndices(knot1n, zeta1);
    Inon2 = ActiveBSplineIndices(knot2n, zeta2);

    % Evaluate knot vectors
    [n1, ~] = EvaluateKnot(knot1n);
    [n2, ~] = EvaluateKnot(knot2n);

    % Initialize global shape function arrays
    NonUBS1 = zeros(1, n1);
    DNonUBS1 = zeros(1, n1);
    NonUBS2 = zeros(1, n2);
    DNonUBS2 = zeros(1, n2);

    % Order of derivative
    m = 1;

    % Compute 1D basis functions and derivatives
    NN1 = GetDNUBSf(knot1n,zeta1,m);
    NonUBS1(Inon1)  = NN1(:,1);
    DNonUBS1(Inon1) = NN1(:,2);
    
    NN2 = GetDNUBSf(knot2n,zeta2,m);
    NonUBS2(Inon2)  = NN2(:,1);
    DNonUBS2(Inon2) = NN2(:,2);

    % Initialize 2D NURBS basis function derivatives
    BS = zeros(n1, n2);
    DBS_Dz1 = zeros(n1, n2);
    DBS_Dz2 = zeros(n1, n2);


    % Compute NURBS basis function derivatives
    for i = 1:length(Inon1)
        for j = 1:length(Inon2)
            BS(Inon1(i), Inon2(j)) = NonUBS1(Inon1(i)) * NonUBS2(Inon2(j));            
            DBS_Dz1(Inon1(i), Inon2(j)) = DNonUBS1(Inon1(i)) * NonUBS2(Inon2(j));
            DBS_Dz2(Inon1(i), Inon2(j)) = NonUBS1(Inon1(i)) * DNonUBS2(Inon2(j));
        end
    end
end

