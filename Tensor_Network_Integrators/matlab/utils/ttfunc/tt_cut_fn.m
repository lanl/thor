function X = tt_cut_fn(X, dims, idxs)
% Cut tt X on the dimension in dims using idx in idxs
% idxs has dimension of len(dims)x2. Each row is the first and last index
% of the corresponding dimension

for i = 1:numel(dims)
  d = dims(i);
  if d== ndims(X)
    X{d} = X{d}(:,idxs(i,1):idxs(i,2));
  else
    X{d} = X{d}(:,idxs(i,1):idxs(i,2),:);
  end
end

end