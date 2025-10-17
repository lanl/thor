function x = tt_orthogonolize(x)
  %
  G = core2cell(x);
  %
  r = x.r;
  n = x.n;
  d = numel(G);
  %
  for i=1:d-1
    %
    cr = reshape(G{i},[r(i)*n(i),r(i+1)]);
    %
    [tempq,tempr] = qr(cr,"econ");
    %
    G{i} = reshape(tempq,r(i),n(i),r(i+1));
    %
    if(i<d)
      G{i+1} = tensorprod(tempr,G{i+1},2,1);
    end
    %
  end
  %
  x = cell2core(tt_tensor,G);
  %
end