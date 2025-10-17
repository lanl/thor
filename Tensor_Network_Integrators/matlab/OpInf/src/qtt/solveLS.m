function O = solveLS(D,R,sp_dim,t_dim,gamma,eps_tt)
  %
  D   = tt_pinv(D,D.d-t_dim,eps_tt,gamma);
  %
  G  = cell(3*(sp_dim+1),1);
  GD = core2cell(D);
  GR = core2cell(R);
  %
  G(1:sp_dim+1) = GR(1:sp_dim+1);
  %
  tempR = cell2core(tt_tensor,GR(end-t_dim+1:end));
  tempD = cell2core(tt_tensor,GD(1:1+t_dim-1));
  temp  = squeeze(dot(tempR,tempD));
  %
  G{sp_dim+2} = tensorprod(temp,GD{t_dim+1},2,1);
  %
  G(sp_dim+3:end) = GD(t_dim+2:end);
  %
  O = cell2core(tt_tensor,G);
  %
  O = round(O,eps_tt);
  %
end
%