function y=tt_contract(A,x,eps_tt)
  %
  GA = core2cell(A);
  Gx = core2cell(x);
  %
  dA = length(GA);
  dx = length(Gx);
  %
  temp      = cell2core(tt_tensor,GA(dA-dx+1:dA));
  temp      = dot(temp,x);
  GA{dA-dx} = tensorprod(GA{dA-dx},temp,3,1);
  y         = cell2core(tt_tensor,[GA(1:dA-dx)]);
  y         = round(y, eps_tt);
  %
end