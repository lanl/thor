function DNN = CalDNN_1D(knot1n, Jt)
[Gd1, J2d1, alphad1] = getIntegralData1D(knot1n);
DNN = 0;

for ii = 1:length(Gd1)
    [N, DN] = GetBspline1Df(knot1n, Gd1(ii));
    % Multiply with factor of numerical integral
    DNN = DNN + alphad1(ii) * J2d1(ii) * Jt(ii) * Tproduct(DN, N);
end

end