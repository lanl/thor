function sig_op = make_sigma_operator_fn(tt_sig)

G = core2cell(tt_sig);
r = tt_sig.r;
n = tt_sig.n;
d = tt_sig.d;

G_op = cell(1,d);
for k = 1:d
  tempten = zeros(r(k),n(k),n(k),r(k+1));
  for i = 1:r(k)
    for j = 1:r(k+1)
      tempten(i,:,:,j) = diag(G{k}(i,:,j));
    end
  end
  G_op{k} = tempten;
end

sig_op = cell2core(tt_matrix,G_op);

end