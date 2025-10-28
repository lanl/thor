function [cr] = get_core(tt, n)
% return n-th core of a tensor
%   - tt: tt-tensor
%   - n : core number
  d = tt.d;
  cr = reshape(tt.core(tt.ps(n):tt.ps(n+1)-1),[tt.r(n),tt.n(n),tt.r(n+1)]);
end

