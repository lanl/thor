function [ratio, tt_numel, full_numel] = compress_ratio_tt(tt, r,n)
%this function will compute the storage ratio of a tt compared to full grid

r = tt.r;
if isa(tt,'tt_matrix')
  n = tt.n.*tt.m;
elseif isa(tt,'tt_tensor')
  n = tt.n;
end
tt_numel = 0;
for i = 1:numel(n)
  tt_numel = tt_numel + r(i)*n(i)*r(i+1); 
end
full_numel = prod(n);
ratio = tt_numel/full_numel;
