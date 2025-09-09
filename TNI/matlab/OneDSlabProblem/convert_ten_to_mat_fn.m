function A = convert_ten_to_mat_fn(A,nx,nmu,nE)
A = double(reshape(permute(A,[1,3,5,2,4,6]),[nx*nmu*nE,nx*nmu*nE]));
end