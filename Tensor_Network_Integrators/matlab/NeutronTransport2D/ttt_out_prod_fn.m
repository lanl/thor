function T = ttt_out_prod_fn(C)
% C is an array of matrix/tensor

T = tensor(C{1});
for i = 2:numel(C)
 T = ttt(T,tensor(C{i}));
end