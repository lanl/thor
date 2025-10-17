function [s1M, s2M] = GetNRelate(k1,k2)

[n1, p1] = EvaluateKnot(k1);
[n2, p2] = EvaluateKnot(k2);

N1span = [k1(1:n1);k1((p1+2):(p1+n1+1))]';
N2span = [k2(1:n2);k2((p2+2):(p2+n2+1))]';
Ek = GetEknot(k2);

s1 = zeros(length(Ek),n1);
s2 = zeros(length(Ek),n2);

s1M = zeros(length(Ek),p1+1);
s2M = zeros(length(Ek),p2+1);

for i=1:length(Ek)
  for j=1:n1
   if CheckIn(Ek(i,:), N1span(j,:)) == 1
     s1(i, j) = j;
   end
  end
  for k=1:n2
   if CheckIn(Ek(i,:), N2span(k,:)) == 1
     s2(i, k) = k;
   end
  end
end
 
for i=1:length(Ek)
    s1M(i,:) = ZeroOff(s1(i,:));
    s2M(i,:) = ZeroOff(s2(i,:));
end

end

