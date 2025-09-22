function [N3D, DND1, DND2, DND3] = GetD1NURBS3Df(k1, k2, k3, xi1, xi2, xi3, pp1, pp2, pp3, we)
m = 1;

NUB1     = GetDNUBSf(k1, xi1, m);
NonUBS1  = NUB1(:,1);
DNonUBS1 = NUB1(:,2);

NUB2     = GetDNUBSf(k2, xi2, m);
NonUBS2  = NUB2(:,1);
DNonUBS2 = NUB2(:,2);

NUB3     = GetDNUBSf(k3, xi3, m);
NonUBS3  = NUB3(:,1);
DNonUBS3 = NUB3(:,2);

N3D      = zeros((pp1 + 1)*(pp2 + 1)*(pp3 + 1), 1);
DND1     = zeros((pp1 + 1)*(pp2 + 1)*(pp3 + 1), 1);
DND2     = zeros((pp1 + 1)*(pp2 + 1)*(pp3 + 1) ,1);
DND3     = zeros((pp1 + 1)*(pp2 + 1)*(pp3 + 1) ,1);
%%
N1N2N3 = Tproduct3(NonUBS1, NonUBS2, NonUBS3);
N1N2N3 = N1N2N3(:); 
%%
DN1N2N3 = Tproduct3(DNonUBS1, NonUBS2, NonUBS3);
DN1N2N3 = DN1N2N3(:);

N1DN2N3 = Tproduct3(NonUBS1, DNonUBS2, NonUBS3);
N1DN2N3 = N1DN2N3(:);

N1N2DN3 = Tproduct3(NonUBS1, NonUBS2, DNonUBS3);
N1N2DN3 = N1N2DN3(:); 
%%

W   = we*(N1N2N3);
DW1 = we*(DN1N2N3);
DW2 = we*(N1DN2N3);
DW3 = we*(N1N2DN3);

m1 = pp1+1;
m2 = pp2+1;
m3 = pp3+1;

for i1=1:m1
    for i2=1:m2
        for i3=1:m3
            idx = m1*(i2-1) + i1 + m1*m2*(i3-1);
            N3D(idx)  = NonUBS1(i1)*NonUBS2(i2)*NonUBS3(i3)*we(idx)/W;
            DND1(idx) = (we(idx)/W^2)*...
                (DNonUBS1(i1)*NonUBS2(i2)*NonUBS3(i3)*W - NonUBS1(i1)*NonUBS2(i2)*NonUBS3(i3)*DW1);
            DND2(idx) = (we(idx)/W^2)*...
                (NonUBS1(i1)*DNonUBS2(i2)*NonUBS3(i3)*W - NonUBS1(i1)*NonUBS2(i2)*NonUBS3(i3)*DW2);
            DND3(idx) = (we(idx)/W^2)*...
                (NonUBS1(i1)*NonUBS2(i2)*DNonUBS3(i3)*W - NonUBS1(i1)*NonUBS2(i2)*NonUBS3(i3)*DW3);
        end
    end
end


end