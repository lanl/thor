function S = GetNURBS3DDomain(knot1, knot2, knot3)
[n1, q1] = EvaluateKnot(knot1);
[n2, q2] = EvaluateKnot(knot2);
[n3, q3] = EvaluateKnot(knot3);
S        = zeros(n1*n2*n3,6);

% The order of extraction is 1st-2nd-3rd direction
for i1=1:n1
    for i2=1:n2
        for i3=1:n3
            idx = n1*(i2-1) + i1 + n1*n2*(i3-1);
            S(idx, :) = [knot1(i1), knot1(i1 + q1 + 1),...
                         knot2(i2), knot2(i2 + q2 + 1),...
                         knot3(i3), knot3(i3 + q3 + 1)];
        end
    end
end

end
