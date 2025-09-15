function Hqtt = convert_H_to_qtt_fn(Htt, tol)
% d is a list of dimension that need to convert to qtt
G = Htt;
Gnew = {};
for i = 1:2
  A = squeeze(G{i});
  n = size(A,1);
  l = log2(n);
  Aqtt = tt_matrix(A,tol,2*ones(1,l),2*ones(1,l));
  Gnew = [Gnew; core2cell(Aqtt)];
end
%add the third dimension
Gnew = [Gnew;G{3}];

Hqtt = cell2core(tt_matrix,Gnew);
