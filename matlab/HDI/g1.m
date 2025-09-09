function fx = g1(A)

persistent Count
if isempty(Count)
  Count = 0;
end

if nargin == 0  % Reply counter and reset it
  fx    = Count;
  Count = 0;
else
  c  = 0.5;
  w  = 0;
  fx = cos(2*pi*w + 2*pi*sum(c*A,2)) ;
  fx = fx.*fx;
  Count = Count + 1;
end

return;
end