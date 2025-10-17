function activeIndices = ActiveBSplineIndices(knot, zeta)
    [~, p] = EvaluateKnot(knot);
    i = sum(knot <= zeta);
    activeIndices = (i - p - 1 : i - 1) + 1;
end
