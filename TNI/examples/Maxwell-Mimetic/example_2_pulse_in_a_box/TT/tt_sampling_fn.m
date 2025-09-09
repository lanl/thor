function X = tt_sampling_fn(X,C)
% X is the tensor
% C is a cell array containing indices in each dimension

G = core2cell(X);
for i = 1:numel(C)
  G{i} = G{i}(:,C{i},:);
end
X = cell2core(tt_tensor,G);