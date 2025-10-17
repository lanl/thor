function L2 = CalH1L2t_Lshape(utt, knot1n, knot2n, knot3n, ctpxn, ctpyn, ctpzn)
    %Gauss data
    [Gp1, J2d1, alphad1] = getIntegralData1D(knot1n);
    [Gp2, J2d2, alphad2] = getIntegralData1D(knot2n);
    [Gp3, J2d3, alphad3] = getIntegralData1D(knot3n);

    nn1 = length(Gp1);
    nn2 = length(Gp2);
    nn3 = length(Gp3);

    L2n = 0;
    L2d = 0;
    for i1 = 1:nn1
        for i2 = 1:nn2
            for i3 =1:nn3
                [NN, DBS_Dz1, DBS_Dz2, DBS_Dz3, Inon1, Inon2, Inon3] = GetBspline3Dn(knot1n, knot2n, knot3n, ...
                                                                                    Gp1(i1), Gp2(i2), Gp3(i3));
                ctpx = ctpxn(Inon1, Inon2, Inon3);
                ctpy = ctpyn(Inon1, Inon2, Inon3);
                ctpz = ctpzn(Inon1, Inon2, Inon3);
                ue = full(utt(Inon1, Inon2, Inon3));
                
                %%
                xg = dot(ctpx(:), NN(:));
                yg = dot(ctpy(:), NN(:));
                zg = dot(ctpz(:), NN(:));
                ug = dot(ue(:), NN(:));
                %Compute Jacobian determinant |J|
                J1 = ComputeJacobian(ctpx, ctpy, ctpz, DBS_Dz1, DBS_Dz2, DBS_Dz3);
                J1 = abs(det(J1));

                % L2 norm
                J2 = J2d1(i1)*J2d2(i2)*J2d3(i3);
                w = alphad1(i1)*alphad2(i2)*alphad3(i3);
                
                % Analytical solution
                ua = (1/(2*pi^2))*sin(pi*xg)*sin(pi*yg);
                % Pointwise error
                du = abs(ug - ua);

                L2n = L2n + w*J1*J2*sqrt(du*du);
                L2d = L2d + w*J1*J2*sqrt(abs(ua)); 
            end
        end
    end
    L2 = L2n/L2d;
end
