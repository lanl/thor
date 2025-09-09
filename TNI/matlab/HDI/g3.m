function fx = g3(A)
%This function takes counts the number of function evaluations f(x) counts
%= counts +1, until we call f then it displayes the number of evaluations
%and 

persistent Count
if isempty(Count)
  Count = 0;
end

if nargin == 0  % Reply counter and reset it
  fx    = Count;
  Count = 0;
  return;
end

d  = size(A,2);
c  = 1;
fx = (1+ sum(c*A,2)).^(-d-1);

Count = Count + 1;

end