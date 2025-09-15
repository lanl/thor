function Xnew = tt_weighted_sum_fn(X,dim,w)
% X is a TT
% weighted sum the 'dim' dimension with weight vector w
G = core2cell(X);
Gnew = cell(ndims(X)-1,1); 
cr = G{dim};
cr2 = tensorprod(cr,w,2,1);
for i = 1:dim-1
  Gnew{i} = G{i};
end
Gnew{dim} = tensorprod(cr2, G{dim+1},2,1);
for i = (dim+1):(ndims(X)-1)
  Gnew{i} = G{i+1};
end

Xnew = cell2core(tt_tensor,Gnew);