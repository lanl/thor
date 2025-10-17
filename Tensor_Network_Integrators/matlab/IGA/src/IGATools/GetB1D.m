function [NN, B1D, DetJ1] = GetB1D(k1,k2,xi1,xi2,we,pp1,pp2,ctpxe,ctpye)
[NN, DND1, DND2] = GetD1NURBS2Df(k1,k2,xi1,xi2,pp1,pp2,we);
B1D = zeros(2,(pp1+1)*(pp2+1));

DXD1 = [ctpxe*DND1,ctpye*DND1];
DXD2 = [ctpxe*DND2,ctpye*DND2];
% First Jacobian
DetJ1 = abs(DXD1(1)*DXD2(2) - DXD1(2)*DXD2(1));
% Derivative of Jacobian
Dxi1Dx1 = +DXD2(2)/DetJ1;
Dxi2Dx1 = -DXD1(2)/DetJ1;
Dxi1Dx2 = -DXD2(1)/DetJ1;
Dxi2Dx2 = +DXD1(1)/DetJ1;

xxterm1 =DND1(:)*Dxi1Dx1 + DND2(:)*Dxi2Dx1;
yyterm2 =DND1(:)*Dxi1Dx2 + DND2(:)*Dxi2Dx2;

B1D(1,1:length(B1D))=xxterm1;
B1D(2,1:length(B1D))=yyterm2;

end