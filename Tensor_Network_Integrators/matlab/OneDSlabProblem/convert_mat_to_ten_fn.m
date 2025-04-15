function A = convert_mat_to_ten_fn(A,nx,nmu,nE)
A = reshape(permute(A,[1,4,2,5,3,6]),[nx,nx,nmu,nmu,nE,nE]);
end