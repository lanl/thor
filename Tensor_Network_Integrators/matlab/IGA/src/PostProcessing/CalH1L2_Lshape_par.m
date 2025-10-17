function L2 = CalH1L2_Lshape_par(ue, ctpxe, ctpye, ctpze, k1, k2, k3, pp1, pp2, pp3, we, m2D)
    % --- Configure the parallel pool ---
    maxAllowed = parcluster('local').NumWorkers;
    desired = min(10, maxAllowed);  % Cap at 15 or available limit    

    p = gcp('nocreate');
    if isempty(p) || p.NumWorkers ~= desired
        if ~isempty(p)
            delete(p);
        end
        pool = parpool(desired);
        pool.IdleTimeout = 120;  % in minutes
    end

    % --- Gauss points and weights ---
    [Gpx, wx] = GaussPoint(pp1 + 1);
    [Gpy, wy] = GaussPoint(pp2 + 1);
    [Gpz, wz] = GaussPoint(pp3 + 1);

    wx = wx(:); wy = wy(:); wz = wz(:);
    [nn, ~] = size(ctpxe);

    L2n_k = zeros(nn, 1);  % parallel accumulators
    L2d_k = zeros(nn, 1);

    pi2 = pi^2;

    % --- Parallel loop over elements ---
    parfor k = 1:nn
        Zetax1 = m2D(k, 1); Zetax2 = m2D(k, 2);
        Zetay1 = m2D(k, 3); Zetay2 = m2D(k, 4);
        Zetaz1 = m2D(k, 5); Zetaz2 = m2D(k, 6);

        J2 = ((Zetax2 - Zetax1) * (Zetay2 - Zetay1) * (Zetaz2 - Zetaz1)) / 8;

        L2n_local = 0;
        L2d_local = 0;

        for ii = 1:length(Gpx)
            zeta1 = (Zetax1 + Zetax2 + Gpx(ii) * (Zetax2 - Zetax1)) / 2;
            wxi = wx(ii);
            for jj = 1:length(Gpy)
                zeta2 = (Zetay1 + Zetay2 + Gpy(jj) * (Zetay2 - Zetay1)) / 2;
                wyj = wy(jj);
                for kk = 1:length(Gpz)
                    zeta3 = (Zetaz1 + Zetaz2 + Gpz(kk) * (Zetaz2 - Zetaz1)) / 2;
                    wzk = wz(kk);

                    % Evaluate basis function and Jacobian
                    [NN, ~, J1] = GetB1D_3D(k1, k2, k3, zeta1, zeta2, zeta3, ...
                                            we(k,:), pp1, pp2, pp3, ...
                                            ctpxe(k,:), ctpye(k,:), ctpze(k,:));

                    xe = ctpxe(k,:) * NN;
                    ye = ctpye(k,:) * NN;
                    ug = ue(k,:) * NN;

                    ua = (1 / (2 * pi2)) * sin(pi * xe) * sin(pi * ye);
                    du = ug - ua;

                    w = wxi * wyj * wzk;
                    L2n_local = L2n_local + w * J1 * J2 * abs(du);
                    L2d_local = L2d_local + w * J1 * J2 * sqrt(abs(ua));
                end
            end
        end

        L2n_k(k) = L2n_local;
        L2d_k(k) = L2d_local;
    end

    L2 = sum(L2n_k) / sum(L2d_k);
end
