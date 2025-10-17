function [Ke, Me] = GetKLaplace(ctpxe, ctpye, k1, k2, pp1, pp2, we, m2D)
[Gpx, wx] = GaussPoint(pp1+1);
[Gpy, wy] = GaussPoint(pp2+1);
% The last index is NOE
Ke = zeros((pp1+1)*(pp2+1),(pp1+1)*(pp2+1),length(ctpxe));
Me = zeros((pp1+1)*(pp2+1),(pp1+1)*(pp2+1),length(ctpxe));

for k=1:length(ctpxe)
      Zetax1 = m2D(k, 1);
      Zetax2 = m2D(k, 2);
      Zetay1 = m2D(k, 3);
      Zetay2 = m2D(k, 4); 
  for ii=1:length(Gpx)
    for jj=1:length(Gpy)
      [NN, B1D, J1] = GetB1D(k1,k2,...
          (Zetax1 + Zetax2 + Gpx(ii)*(Zetax2 - Zetax1))/2,...
          (Zetay1 + Zetay2 + Gpy(jj)*(Zetay2 - Zetay1))/2,...
          we(k,:),pp1,pp2,ctpxe(k,:),ctpye(k,:));
      J2 = ((Zetax2 - Zetax1)*(Zetay2 - Zetay1))/4;
      Ke(:,:,k) = Ke(:,:,k) + wx(ii)*wy(jj)*J1*J2*(B1D'*B1D);
      Me(:,:,k) = Me(:,:,k) + wx(ii)*wy(jj)*J1*J2*Tproduct(NN,NN);      
    end
  end
end
end