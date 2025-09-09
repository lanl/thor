function fx = g5(A)

persistent Count
if isempty(Count)
  Count = 0;
end

if nargin == 0  % Reply counter and reset it
  fx    = Count;
  Count = 0;
else
  c  = 1;
  fx = exp( - sum(c*A,2) );
  Count = Count + 1;
end

return;
end