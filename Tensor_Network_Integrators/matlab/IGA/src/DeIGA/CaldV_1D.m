function dV = CaldV_1D(knot1n, Jt)
[Gd1, J2d1, alphad1] = getIntegralData1D(knot1n);
dV = 0;

for ii = 1:length(Gd1)
    % Multiply with factor of numerical integral
    dV = dV + alphad1(ii) * J2d1(ii) * Jt(ii);
end

end