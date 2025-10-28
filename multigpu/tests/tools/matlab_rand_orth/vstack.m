function [V] = vstack(tt, n)
% hstack: return n-th core of a tensor in vstack form
%   - tt: tt-tensor
%   - n : core number
  V = reshape(tt.core(tt.ps(n):tt.ps(n+1)-1),[tt.r(n)*tt.n(n),tt.r(n+1)]);
end

