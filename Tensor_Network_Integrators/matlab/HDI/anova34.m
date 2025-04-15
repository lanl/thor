function fx = anova34(A)

persistent Count
if isempty(Count)
  Count = 0;
end

if nargin == 0  % Reply counter and reset it
  fx    = Count;
  Count = 0;
else
  ak = 1;
  fx =  prod((abs(4*A - 3) - 0.25),2);
  Count = Count + 1;
end

return;
end