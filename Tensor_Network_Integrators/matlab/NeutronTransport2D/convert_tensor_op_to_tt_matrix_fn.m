function Htt = convert_tensor_op_to_tt_matrix_fn(H,d)
permuted_dim = [1:2:d,2:2:d];
H_for_tt = permute(double(H),permuted_dim);
Htt = tt_matrix(H_for_tt);
end