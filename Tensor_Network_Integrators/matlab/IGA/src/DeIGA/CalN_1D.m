%This function is used for scalar variable like PF, temperature,
%etc
function NN = CalN_1D(knot1n, Jt)
[Gd1, J2d1, alphad1] = getIntegralData1D(knot1n);
NN = 0;

for ii = 1:length(Gd1)
    [N, ~] = GetBspline1Df(knot1n, Gd1(ii));
    % Multiply with factor of numerical integral
    NN = NN + alphad1(ii) * J2d1(ii) * Jt(ii) * N;
end

end