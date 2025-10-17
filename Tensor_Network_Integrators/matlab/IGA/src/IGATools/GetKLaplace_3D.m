function [Ke, Me] = GetKLaplace_3D(ctpxe, ctpye, ctpze, k1, k2, k3, pp1, pp2, pp3, we, m3D)
[Gpx, wx] = GaussPoint(pp1+1);
[Gpy, wy] = GaussPoint(pp2+1);
[Gpz, wz] = GaussPoint(pp3+1);

[nn, ~] = size(ctpxe);
% The last index is NOE
Ke = zeros((pp1+1)*(pp2+1)*(pp3+1), (pp1+1)*(pp2+1)*(pp3+1), nn);
Me = zeros((pp1+1)*(pp2+1)*(pp3+1), (pp1+1)*(pp2+1)*(pp3+1), nn);

for k=1:nn
    Zetax1 = m3D(k, 1);
    Zetax2 = m3D(k, 2);
    Zetay1 = m3D(k, 3);
    Zetay2 = m3D(k, 4);
    Zetaz1 = m3D(k, 5);
    Zetaz2 = m3D(k, 6);
    for ii=1:length(Gpx)
        for jj=1:length(Gpy)
            for kk=1:length(Gpz)
                [NN, B1D, J1] = GetB1D_3D(k1, k2, k3,...
                                         (Zetax1 + Zetax2 + Gpx(ii)*(Zetax2 - Zetax1))/2,...
                                         (Zetay1 + Zetay2 + Gpy(jj)*(Zetay2 - Zetay1))/2,...
                                         (Zetaz1 + Zetaz2 + Gpz(kk)*(Zetaz2 - Zetaz1))/2,...
                                         we(k,:), pp1, pp2, pp3, ctpxe(k,:), ctpye(k,:), ctpze(k,:));
                 J2 = ((Zetax2 - Zetax1)*(Zetay2 - Zetay1)*(Zetaz2 - Zetaz1))/8;
                 Ke(:,:,k) = Ke(:,:,k) + wx(ii)*wy(jj)*wz(kk)*J1*J2*(B1D'*B1D);
                 Me(:,:,k) = Me(:,:,k) + wx(ii)*wy(jj)*wz(kk)*J1*J2*Tproduct(NN,NN);      
            end
        end
    end

end
end