function [Ke, Me] = GetKM_Elas_3D(ctpxe, ctpye, ctpze, k1, k2, k3, pp1, pp2, pp3, we, m2D)
[Gpx, wx] = GaussPoint(pp1+1);
[Gpy, wy] = GaussPoint(pp2+1);
[Gpz, wz] = GaussPoint(pp3+1);
CC = [1.0, 0.7, 0.7, 0.0, 0.0, 0.0;...
      0.7, 1.0, 0.7, 0.0, 0.0, 0.0;...
      0.7, 0.7, 1.0, 0.0, 0.0, 0.0;...
      0.0, 0.0, 0.0, 0.5, 0.0, 0.0;...
      0.0, 0.0, 0.0, 0.0, 0.5, 0.0;...
      0.0, 0.0, 0.0, 0.0, 0.0, 0.5];

% The last index is NOE
Ke = zeros(3*(pp1+1)*(pp2+1)*(pp3+1),3*(pp1+1)*(pp2+1)*(pp3+1),length(ctpxe));
Me = zeros(3*(pp1+1)*(pp2+1)*(pp3+1),3*(pp1+1)*(pp2+1)*(pp3+1),length(ctpxe));

for k=1:length(ctpxe)
    Zetax1 = m2D(k, 1);
    Zetax2 = m2D(k, 2);
    Zetay1 = m2D(k, 3);
    Zetay2 = m2D(k, 4);
    Zetaz1 = m2D(k, 5);
    Zetaz2 = m2D(k, 6);
    for ii=1:length(Gpx)
        for jj=1:length(Gpy)
            for kk=1:length(Gpz)
                [NN, B3D, J1] = GetB3D_3D(k1, k2, k3,...
                                         (Zetax1 + Zetax2 + Gpx(ii)*(Zetax2 - Zetax1))/2,...
                                         (Zetay1 + Zetay2 + Gpy(jj)*(Zetay2 - Zetay1))/2,...
                                         (Zetaz1 + Zetaz2 + Gpz(kk)*(Zetaz2 - Zetaz1))/2,...
                                         we(k,:), pp1, pp2, pp3, ctpxe(k,:), ctpye(k,:), ctpze(k,:));
                 J2 = ((Zetax2 - Zetax1)*(Zetay2 - Zetay1)*(Zetaz2 - Zetaz1))/8;
                 Ke(:,:,k) = Ke(:,:,k) + wx(ii)*wy(jj)*wz(kk)*J1*J2*(B3D'*CC*B3D);
                 Me(:,:,k) = Me(:,:,k) + wx(ii)*wy(jj)*wz(kk)*J1*J2*(NN'*NN);      
            end
        end
    end

end
end