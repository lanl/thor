function [H] = hstack(tt, n)
% hstack: return n-th core of a tensor in hstack form
%   - tt: tt-tensor
%   - n : core number
  d = tt.d; 
  cr = reshape(tt.core(tt.ps(n):tt.ps(n+1)-1),[tt.r(n),tt.n(n),tt.r(n+1)]);
  H = reshape(permute(cr,[1,3,2]),[tt.r(n),tt.n(n)*tt.r(n+1)]);
end

