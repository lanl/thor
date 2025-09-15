function tt_sig_ext = extend_tt_sigmaSF_fn(sigsmat,tol)

Ns = size(sigsmat);
d = numel(Ns) - 2;
S = reshape(sigsmat,[Ns(1:end-2), Ns(end)*Ns(end-1)]);
tt_sigs = tt_tensor(S,tol);

n = tt_sigs.n;
r = tt_sigs.r;
G = cell(1,3);
for i = 1:d
  crc = tt_sigs{i};
  crc = reshape(crc,[r(i)*n(i),r(i+1)]);
  I = diag_ten(r(i+1),n(i));
  I = reshape(I,r(i+1),[]);
  crc = tensorprod(crc,I,2,1);
  crc = reshape(crc,[r(i),n(i),n(i),r(i+1)]);
  G{i} = crc;
end
G{d+1} = reshape(tt_sigs{d+1},r(d+1),Ns(end-1),Ns(end),r(d+2));

tt_sig_ext = cell2core(tt_matrix,G);