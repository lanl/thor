function V_tt = tt_expand_vector(v, N, idim)
% place the vector v into the 'idim' dimension of the rank-1 tensor train with
% dimensions N
ndim = numel(N);
G = cell(1, ndim);
for i = 1:ndim
  G{i} = ones(1,N(i));
end
G{idim} = reshape(v,1,[]);
V_tt = cell2core(tt_tensor,G);