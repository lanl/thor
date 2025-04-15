function ttmat = matrices_to_qtt_matrix_fn(C,tol)
% C is a cell array of matrices
d = numel(C);
G = {};
for i = 1:d
  %   G{i} = reshape(C{i},[1,size(C{i}),1]);
  A = double(C{i});
  n = size(A,1);
  l = log2(n);
  Aqtt = tt_matrix(A,tol,2*ones(1,l),2*ones(1,l));
  G = [G;core2cell(Aqtt)];
end
%% energy dimension
% G = [G;reshape(double(C{d}),[1,size(C{d}),1])];

%%
ttmat = cell2core(tt_matrix,G);