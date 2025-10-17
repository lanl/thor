function S = GetNURBS2DDomain(knot1,knot2)
[n, q] = EvaluateKnot(knot1);
[m, p] = EvaluateKnot(knot2);
S      = zeros(n*m,4);
    for l=1:m
        for k=1:n    
        S(n*(l-1)+k,:) = [knot1(k), knot1(k + q + 1), knot2(l), knot2(l + p + 1)];
        end
    end
end

