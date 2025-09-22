%This function is used for 2D variable like displacement with u_x, u_y
function NN = CalNN_2D(knot1n, Jt)
[Gd1, J2d1, alphad1] = getIntegralData1D(knot1n);
NN = 0;

for ii = 1:length(Gd1)
    [N, ~] = GetBspline1Df(knot1n, Gd1(ii));
    N2 = NtoN2(N);
    % Multiply with factor of numerical integral
    NN = NN + alphad1(ii) * J2d1(ii) * Jt(ii) * (N2'*N2);
end

end