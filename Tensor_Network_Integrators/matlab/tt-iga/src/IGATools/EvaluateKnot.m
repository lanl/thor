function [n,p] = EvaluateKnot(knot)
p = 0;
for i=1:length(knot)
 if knot(i) == 0
     p = p + 1;
 end
end
p = p - 1;
n = length(knot) - p - 1;

end

