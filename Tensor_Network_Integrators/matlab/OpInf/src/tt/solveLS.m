function O = solveLS(D,R,gamma,eps_tt)
  %
  D   = tt_pinv(D,8,eps_tt,gamma);
  %
  G  = cell(9,1);
  GD = core2cell(D);
  GR = core2cell(R);
  %
  G(1:4) = GR(1:4);
  %
  if(ndims(GD{1})==2)
    temp = squeeze(GR{end})*GD{1}';
  else
    temp = squeeze(GR{end})*squeeze(GD{1});
  end
  %
  G{5} = tensorprod(temp,GD{2},2,1);
  %
  G(6:12) = GD(3:9);
  %
  O = cell2core(tt_tensor,G);
  %
end
%