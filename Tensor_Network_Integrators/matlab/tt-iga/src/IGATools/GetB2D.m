function [B2D,DetJ1] = GetB2D(k1,k2,xi1,xi2,we,pp1,pp2,ctpxe,ctpye)
[~, DND1, DND2] = GetD1NURBS2Df(k1,k2,xi1,xi2,pp1,pp2,we);
B2D = zeros(3,2*(pp1+1)*(pp2+1));
DXD1 = [ctpxe*DND1,ctpye*DND1];
DXD2 = [ctpxe*DND2,ctpye*DND2];
% First Jacobian
DetJ1 = abs(DXD1(1)*DXD2(2) - DXD1(2)*DXD2(1));
% Derivative of Jacobian
Dxi1Dx1 =+DXD2(2)/DetJ1;
Dxi2Dx1 =-DXD1(2)/DetJ1;
Dxi1Dx2 =-DXD2(1)/DetJ1;
Dxi2Dx2 =+DXD1(1)/DetJ1;

xxterm1 =DND1(:)*Dxi1Dx1+DND2(:)*Dxi2Dx1;
yyterm2 =DND1(:)*Dxi1Dx2+DND2(:)*Dxi2Dx2;
xyterm1 =DND1(:)*Dxi1Dx2+DND2(:)*Dxi2Dx2;
xyterm2 =DND1(:)*Dxi1Dx1+DND2(:)*Dxi2Dx1;

B2D(1,1:2:length(B2D)-1)=xxterm1;
B2D(2,2:2:length(B2D)-0)=yyterm2;

B2D(3,1:2:length(B2D)-1)=xyterm1;
B2D(3,2:2:length(B2D)-0)=xyterm2;

end