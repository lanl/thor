function [L2, u_error] = CalH1L2_Ring(u_in, u_out, r_int, r_out,... 
                      ue, ctpxe, ctpye, ctpze, k1, k2, k3, pp1, pp2, pp3, we, m2D)
[Gpx, wx] = GaussPoint(pp1+1);
[Gpy, wy] = GaussPoint(pp2+1);
[Gpz, wz] = GaussPoint(pp3+1);
L2n = 0; L2d = 0;
[nn, n1] = size(ctpxe);
u_error = zeros(nn, n1);
for k=1:nn
    cnt = 1;
    Zetax1 = m2D(k, 1);
    Zetax2 = m2D(k, 2);
    Zetay1 = m2D(k, 3);
    Zetay2 = m2D(k, 4);
    Zetaz1 = m2D(k, 5);
    Zetaz2 = m2D(k, 6);
    for ii=1:length(Gpx)
        for jj=1:length(Gpy)
            for kk=1:length(Gpz)
                [NN, ~, J1] = GetB1D_3D(k1, k2, k3,...
                                         (Zetax1 + Zetax2 + Gpx(ii)*(Zetax2 - Zetax1))/2,...
                                         (Zetay1 + Zetay2 + Gpy(jj)*(Zetay2 - Zetay1))/2,...
                                         (Zetaz1 + Zetaz2 + Gpz(kk)*(Zetaz2 - Zetaz1))/2,...
                                         we(k,:), pp1, pp2, pp3, ctpxe(k,:), ctpye(k,:), ctpze(k,:));
                 J2 = ((Zetax2 - Zetax1)*(Zetay2 - Zetay1)*(Zetaz2 - Zetaz1))/8;
                 xe = ctpxe(k,:)*NN;
                 ye = ctpye(k,:)*NN;
                 ze = ctpze(k,:)*NN;
                 ug = ue(k,:)*NN;

                 %Analytical solution
                 r = sqrt(xe^2 + ye^2);
                 ua = (u_in*log(r_out/r)+ u_out*log(r/r_int))/log(r_out/r_int);
                 du = abs(ug - ua);
                 u_error(k, cnt) = du;
                 cnt = cnt +1;
                 %dug = B1D.ue(kk);
                 L2n = L2n + wx(ii)*wy(jj)*wz(kk)*J1*J2*sqrt(du*du);
                 L2d = L2d + wx(ii)*wy(jj)*wz(kk)*J1*J2*sqrt(abs(ua)); 
            end
        end
    end

end
L2 = L2n/L2d;
end