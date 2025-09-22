function Eknot = GetEknot(knot)
Rk    = ReduceKnotFull(knot);
Eknot = zeros(length(Rk)-1,2);
for i=1:length(Eknot)
Eknot(i,:) = [Rk(i),Rk(i+1)];
end
end

