function ttmat = matrices_to_tt_matrix_fn(C)
% C is a cell array of matrices
d = numel(C);
G = cell(1,d);
for i = 1:d
  G{i} = reshape(C{i},[1,size(C{i}),1]);
end
ttmat = cell2core(tt_matrix,G);