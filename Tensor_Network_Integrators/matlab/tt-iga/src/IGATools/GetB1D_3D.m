function [NN, B1D, J] = GetB1D_3D(k1, k2, k3, xi1, xi2, xi3, we, pp1, pp2, pp3, ctpxe, ctpye, ctpze)
[NN, DND1, DND2, DND3] = GetD1NURBS3Df(k1, k2, k3, xi1, xi2, xi3, pp1, pp2, pp3, we);
B1D = zeros(3,(pp1+1)*(pp2+1)*(pp3+1));

DX1Dz = [ctpxe*DND1, ctpxe*DND2, ctpxe*DND3];
DX2Dz = [ctpye*DND1, ctpye*DND2, ctpye*DND3];
DX3Dz = [ctpze*DND1, ctpze*DND2, ctpze*DND3];

% First Jacobian
j = [DX1Dz; DX2Dz; DX3Dz];
J = abs(det(j));


Dz1DX1 = (1/J)*(-j(2,3)*j(3,2) + j(2,2)*j(3,3));
Dz1DX2 = (1/J)*(j(1,3)*j(3,2) - j(1,2)*j(3,3));
Dz1DX3 = (1/J)*(-j(1,3)*j(2,2) + j(1,2)*j(2,3));

Dz2DX1 = (1/J)*(j(2,3)*j(3,1) - j(2,1)*j(3,3));
Dz2DX2 = (1/J)*(-j(1,3)*j(3,1) + j(1,1)*j(3,3));
Dz2DX3 = (1/J)*(j(1,3)*j(2,1) - j(1,1)*j(2,3));

Dz3DX1 = (1/J)*(-j(2,2)*j(3,1) + j(2,1)*j(3,2));
Dz3DX2 = (1/J)*(j(1,2)*j(3,1) - j(1,1)*j(3,2));
Dz3DX3 = (1/J)*(-j(1,2)*j(2,1) + j(1,1)*j(2,2));


DNDX1 = DND1*Dz1DX1 + DND2*Dz2DX1 + DND3*Dz3DX1;
DNDX2 = DND1*Dz1DX2 + DND2*Dz2DX2 + DND3*Dz3DX2;
DNDX3 = DND1*Dz1DX3 + DND2*Dz2DX3 + DND3*Dz3DX3;


B1D(1,1:length(B1D)) = DNDX1;
B1D(2,1:length(B1D)) = DNDX2;
B1D(3,1:length(B1D)) = DNDX3;

end