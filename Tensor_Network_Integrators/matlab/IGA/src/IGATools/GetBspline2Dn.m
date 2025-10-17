function [BS, DBS_Dz1, DBS_Dz2, Inon1, Inon2] = GetBspline2Dn(knot1n, knot2n, Gp1, Gp2)

[BS1, DBS1, Inon1] = GetBspline1Dn(knot1n, Gp1);
[BS2, DBS2, Inon2] = GetBspline1Dn(knot2n, Gp2);

n1 = length(BS1);
n2 = length(BS2);
% Initialize 2D Bspline and its derivatives
BS = zeros(n1, n2);
DBS_Dz1 = zeros(n1, n2);
DBS_Dz2 = zeros(n1, n2);


    % Compute NURBS basis function derivatives
    for i = 1:length(Inon1)
        for j = 1:length(Inon2)
            BS(i, j) = BS1(i) * BS2(j);            
            DBS_Dz1(i, j) = DBS1(i) * BS2(j);
            DBS_Dz2(i, j) = BS1(i) * DBS2(j);
        end
    end
end
