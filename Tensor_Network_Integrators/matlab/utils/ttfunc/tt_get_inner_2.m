function X = tt_get_inner_2(X,dims)

N = X.n;
G = core2cell(X);
G1 = cell(1,numel(N));

if nargin<2
  for i = 1:numel(N)
    G1{i} = G{i}(:,[2:end-1],:);
  end
else
  for i = 1:numel(N)
    if sum(i==dims)>0
      G1{i} = G{i}(:,[2:end-1],:);
    else
      G1{i} = G{i};
    end
  end
end


X = cell2core(tt_tensor,G1);