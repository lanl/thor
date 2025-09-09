function fx = g2( A )
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

wl = 0.;
cl = 1.;
s  = 4/pi;

fx = prod( s*1./( cl.^(2) + (A - wl).^2), 2 );
Count = Count + 1;
end