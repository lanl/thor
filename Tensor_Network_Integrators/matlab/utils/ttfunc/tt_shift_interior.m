function X = tt_shift_interior(X, dim, s, zero_bnd)
%s is the shift index

if nargin<4
  zero_bnd = 0;
end

N = X.n;
G = core2cell(X);
if numel(dim)==1
  G = tt_shift_one_dim(G,dim,s, N);
else
  for i = 1:numel(dim)
    G = tt_shift_one_dim(G,dim(i),s(i), N);
  end
end
X = cell2core(tt_tensor,G);

if zero_bnd; X = tt_set_zero_boundaries(X); end

end
%%
function G = tt_shift_one_dim(G, dim, s, N)

if s>=0
  G{dim}(:,[1:N(dim)-s],:) = G{dim}(:,[1+s:N(dim)],:);
else
  G{dim}(:,[1-s:N(dim)],:) = G{dim}(:,[1:N(dim)+s],:);
end
end