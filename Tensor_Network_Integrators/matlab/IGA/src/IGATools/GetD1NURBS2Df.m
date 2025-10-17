function [N2D, DND1, DND2] = GetD1NURBS2Df(k1,k2,xi1,xi2,pp1,pp2,we)
m=1;
NUB1     = GetDNUBSf(k1, xi1, m);
NonUBS1  = NUB1(:,1);
DNonUBS1 = NUB1(:,2);
NUB2     = GetDNUBSf(k2, xi2, m);
NonUBS2  = NUB2(:,1);
DNonUBS2 = NUB2(:,2);

N2D      = zeros((pp1 + 1)*(pp2 + 1),1);
DND1     = zeros((pp1 + 1)*(pp2 + 1),1);
DND2     = zeros((pp1 + 1)*(pp2 + 1),1);

N1N2 = Tproduct(NonUBS1, NonUBS2);
N1N2 = N1N2(:);

DN1N2 = Tproduct(DNonUBS1, NonUBS2);
DN1N2 = DN1N2(:);

N1DN2 = Tproduct(NonUBS1, DNonUBS2);
N1DN2 = N1DN2(:);

W   = we*(N1N2);
DW1 = we*(DN1N2);
DW2 = we*(N1DN2);

for i=1:(pp1+1)
    for j=1:(pp2+1)
        idx = (pp1+1)*(j-1) + i;
        N2D(idx)  = NonUBS1(i)*NonUBS2(j)*we(idx)/W;
        DND1(idx) = (we(idx)/W^2)*(DNonUBS1(i)*NonUBS2(j)*W - NonUBS1(i)*NonUBS2(j)*DW1);
        DND2(idx) = (we(idx)/W^2)*(NonUBS1(i)*DNonUBS2(j)*W - NonUBS1(i)*NonUBS2(j)*DW2);
    end
end

end