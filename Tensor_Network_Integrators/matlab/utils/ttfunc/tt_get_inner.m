function X = tt_get_inner(X,Cidx)

N = X.n;
G = core2cell(X);
G1 = cell(1,numel(N));

if nargin<2
  for i = 1:numel(N)
    G1{i} = G{i}(:,[2:end-1],:);
  end
else
  for i = 1:numel(N)
    G1{i} = G{i}(:,Cidx{i},:);
  end
end

X = cell2core(tt_tensor,G1);