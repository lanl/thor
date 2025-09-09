function tt_op = make_tt_to_operator(Xtt)
% only for 2D rightnow, 1 energy group right now
G = core2cell(Xtt);
r = Xtt.r;
n = Xtt.n;
d = Xtt.d;

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

tt_op = cell2core(tt_matrix,G_op);

end