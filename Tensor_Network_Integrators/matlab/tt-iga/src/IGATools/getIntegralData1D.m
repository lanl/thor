function [G1D, J1D, Alpha1D] = getIntegralData1D(knot)
    m1D = macro1D(knot);
    [~, p] = EvaluateKnot(knot);

    % Get Gaussian points
    [Gp, Alpha] = GaussPoint(p+1);
    G1D = zeros(length(m1D), length(Gp));
    J1D = zeros(length(m1D), length(Gp));
    Alpha1D = zeros(length(m1D), length(Gp)); 
    for i = 1:length(m1D)
        for j = 1:length(Gp)
            G1D(i, j) = 0.5*(m1D(i, 2) + m1D(i, 1) + Gp(j)*(m1D(i, 2) - m1D(i, 1)));
            J1D(i, j) = 0.5*(m1D(i, 2) - m1D(i, 1));
            Alpha1D(i, j) = Alpha(j);
        end    
    end
    % Vectorize
    G1D = G1D(:);
    J1D = J1D(:);
    Alpha1D =Alpha1D(:);