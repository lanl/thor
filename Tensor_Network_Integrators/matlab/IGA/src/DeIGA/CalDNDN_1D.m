function DNDN = CalDNDN_1D(knot1n, Jt)
[Gd1, J2d1, alphad1] = getIntegralData1D(knot1n);
DNDN = 0;

for ii = 1:length(Gd1)
    [~, DN] = GetBspline1Df(knot1n, Gd1(ii));
    % Multiply with factor of numerical integral
    DNDN = DNDN + alphad1(ii) * J2d1(ii) * Jt(ii) * Tproduct(DN, DN);
end

end