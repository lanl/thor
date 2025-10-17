function [Bezier, sp] = GetBezierExtraction1D(Ck, N, ith, kold, knew)
[s1, s2] = GetNRelate(kold, knew);
Ce = zeros(length(s1),length(s2));
kR = ReduceKnotFull(knew);
sp = [kR(ith) kR(ith+1)];
for i=1:length(s1(ith))
    for j=1:length(s2(ith))
    Ce(i,j) = Ck(s1(ith,i),s2(ith,j));
    end
end

Ne = N(s1(ith,:));
Bezier = inverse(Ce)*Ne;

end

