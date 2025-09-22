function [BS1, DBS1, Inon1] = GetBspline1Df(knot1n, zeta1)
    Inon1 = ActiveBSplineIndices(knot1n, zeta1);
    [n1, ~] = EvaluateKnot(knot1n);
    BS1 = zeros(1, n1);
    DBS1 = zeros(1, n1);

    % Order of derivative
    m = 1;

    % Compute 1D basis functions and derivatives
    NN1 = GetDNUBSf(knot1n,zeta1,m);
    BS1(Inon1)  = NN1(:,1);
    DBS1(Inon1) = NN1(:,2);
 
end
