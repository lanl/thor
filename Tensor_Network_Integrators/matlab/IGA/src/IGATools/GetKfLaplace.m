function [Ke, fe] = GetKfLaplace(ctpxe,ctpye,k1,k2,pp1,pp2,we,m2D, ffunction)
[Gpx,wx] = GaussPoint(pp1);
[Gpy,wy] = GaussPoint(pp2);

Ke = zeros((pp1+1)*(pp2+1),(pp1+1)*(pp2+1),length(ctpxe));
fe = zeros((pp1+1)*(pp2+1),length(ctpxe));

for k = 1:length(ctpxe)
      Zetax1 = m2D(k, 1);
      Zetax2 = m2D(k, 2);
      Zetay1 = m2D(k, 3);
      Zetay2 = m2D(k, 4); 
  for ii = 1:length(Gpx)
    for jj = 1:length(Gpy)
      [NN, B1D, J1] = GetB1D(k1,k2,...
          (Zetax1 + Zetax2 + Gpx(ii)*(Zetax2 - Zetax1))/2,...
          (Zetay1 + Zetay2 + Gpy(jj)*(Zetay2 - Zetay1))/2,...
          we(k,:),pp1,pp2,ctpxe(k,:),ctpye(k,:));
      J2 = ((Zetax2 - Zetax1)*(Zetay2 - Zetay1))/4;
      Ke(:,:,k) = Ke(:,:,k) + wx(ii)*wy(jj)*J1*J2*(B1D'*B1D);
      %fe(:,k) = fe(:,k) + wx(ii)*wy(jj)*J1*J2*NN*ffunction;
      x = ctpxe(k,:)*NN;
      y = ctpye(k,:)*NN;
      fe(:,k) = fe(:,k) + wx(ii)*wy(jj)*J1*J2*NN*(cos(x)*sin(y));
    end
  end
end
end