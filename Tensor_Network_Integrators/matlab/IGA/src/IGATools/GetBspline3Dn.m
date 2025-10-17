function [BS, DBS_Dz1, DBS_Dz2, DBS_Dz3, Inon1, Inon2, Inon3] = GetBspline3Dn(knot1n, knot2n, knot3n, ...
                                                                Gp1, Gp2, Gp3)

[BS1, DBS1, Inon1] = GetBspline1Dn(knot1n, Gp1);
[BS2, DBS2, Inon2] = GetBspline1Dn(knot2n, Gp2);
[BS3, DBS3, Inon3] = GetBspline1Dn(knot3n, Gp3);

n1 = length(BS1);
n2 = length(BS2);
n3 = length(BS3);

% Initialize 2D Bspline and its derivatives
BS = zeros(n1, n2, n3);
DBS_Dz1 = zeros(n1, n2, n3);
DBS_Dz2 = zeros(n1, n2, n3);
DBS_Dz3 = zeros(n1, n2, n3);

    % Compute basis function derivatives
    for i = 1:length(Inon1)
        for j = 1:length(Inon2)
            for k = 1:length(Inon3)
                BS(i, j, k) = BS1(i) * BS2(j) * BS3(k);            
                DBS_Dz1(i, j, k) = DBS1(i) * BS2(j) * BS3(k);
                DBS_Dz2(i, j, k) = BS1(i) * DBS2(j) * BS3(k);
                DBS_Dz3(i, j, k) = BS1(i) * BS2(j) * DBS3(k);
            end
        end
    end

end
