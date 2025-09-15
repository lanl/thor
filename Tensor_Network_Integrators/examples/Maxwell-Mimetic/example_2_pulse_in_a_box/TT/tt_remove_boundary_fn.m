function X = tt_remove_boundary_fn(X, dims)

G = core2cell(X);
for i = 1: numel(dims)
  d = dims(i);
  G{d} = G{d}(:,2:end-1,:);
end
X = cell2core(tt_tensor, G);