function [L2, u_error] = CalH1L2_Ring_par(u_in, u_out, r_int, r_out,... 
                      ue, ctpxe, ctpye, ctpze, k1, k2, k3, pp1, pp2, pp3, we, m2D)

% Gauss points and weights
[Gpx, wx] = GaussPoint(pp1+1);
[Gpy, wy] = GaussPoint(pp2+1);
[Gpz, wz] = GaussPoint(pp3+1);

% Element count and size
[nn, n1] = size(ctpxe);

% Output initialization
u_error = zeros(nn, n1);
L2n_vec = zeros(nn, 1);
L2d_vec = zeros(nn, 1);

% Parallel loop
parfor k = 1:nn
    L2n_local = 0;
    L2d_local = 0;
    u_error_local = zeros(1, n1);
    cnt = 1;

    % Elemental bounds
    Zetax1 = m2D(k, 1); Zetax2 = m2D(k, 2);
    Zetay1 = m2D(k, 3); Zetay2 = m2D(k, 4);
    Zetaz1 = m2D(k, 5); Zetaz2 = m2D(k, 6);

    for ii = 1:length(Gpx)
        for jj = 1:length(Gpy)
            for kk = 1:length(Gpz)
                % Reference to physical space mapping
                xref = (Zetax1 + Zetax2 + Gpx(ii)*(Zetax2 - Zetax1))/2;
                yref = (Zetay1 + Zetay2 + Gpy(jj)*(Zetay2 - Zetay1))/2;
                zref = (Zetaz1 + Zetaz2 + Gpz(kk)*(Zetaz2 - Zetaz1))/2;

                [NN, ~, J1] = GetB1D_3D(k1, k2, k3,...
                                         xref, yref, zref,...
                                         we(k,:), pp1, pp2, pp3, ...
                                         ctpxe(k,:), ctpye(k,:), ctpze(k,:));
                J2 = ((Zetax2 - Zetax1)*(Zetay2 - Zetay1)*(Zetaz2 - Zetaz1))/8;

                % Physical coordinates
                xe = ctpxe(k,:) * NN;
                ye = ctpye(k,:) * NN;
                ze = ctpze(k,:) * NN;

                % Numerical solution
                ug = ue(k,:) * NN;

                % Analytical solution
                r = sqrt(xe^2 + ye^2);
                ua = (u_in*log(r_out/r)+ u_out*log(r/r_int))/log(r_out/r_int);

                % Error computation
                du = abs(ug - ua);
                u_error_local(cnt) = du;

                wprod = wx(ii)*wy(jj)*wz(kk);
                L2n_local = L2n_local + wprod*J1*J2*sqrt(du*du);
                L2d_local = L2d_local + wprod*J1*J2*sqrt(abs(ua));
                cnt = cnt + 1;
            end
        end
    end

    L2n_vec(k) = L2n_local;
    L2d_vec(k) = L2d_local;
    u_error(k,:) = u_error_local;
end

% Final L2 norm
L2 = sum(L2n_vec)/sum(L2d_vec);
end
