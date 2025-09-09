function X = tt_set_zero_boundaries(X,I)
% set the boundary ot X tensor to zeros
% I is a cell array specifying which boundary to set to zeroes
N = X.n;
G = core2cell(X);

if nargin<2
  I = cell(1,numel(N));
  for i = 1:numel(N)
    I{i} = [1,N(i)];
  end
end

for i = 1:numel(N)
  G{i}(:,I{i},:) = 0;
end
X = cell2core(tt_tensor,G);