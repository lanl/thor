function NDN = CalNDN_1D(knot1n, Jt)
[Gd1, J2d1, alphad1] = getIntegralData1D(knot1n);
NDN = 0;

for ii = 1:length(Gd1)
    [N, DN] = GetBspline1Df(knot1n, Gd1(ii));
    % Multiply with factor of numerical integral
    NDN = NDN + alphad1(ii) * J2d1(ii) * Jt(ii) * Tproduct(N, DN);
end

end