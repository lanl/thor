function Xttpinv = tt_pinv(Xtt,coreidx,tt_tol,gamma)
  %
  % This function computes the TT-format of the pseudo-inverse
  % of the matricization of Xtt(Nu x Nv), where
  % Nu = prod(Xtt.n(1:coreidx))
  % Nv = prod(Xtt.n(coreidx+1:end))
  %
  % Step 0: Check regularization parameter
  %
  if(~exist("gamma","var"))
    gamma = 0;
  end
  %
  % Step 1: Get cores
  %
  G = core2cell(Xtt);
  n = Xtt.n;
  r = Xtt.r;
  %
  % Step 2: Left-orthogonolize G1 & G2
  %
  [Gl, Rl] = tt_orthogonolize(G(1:coreidx),n(1:coreidx),r(1:coreidx+1),1);
  %
  % Step 3: Right-orthogonolize G3
  %
  [Gr, Rr] = tt_orthogonolize(G(coreidx+1:end),n(coreidx+1:end),r(coreidx+1:end),-1);
  %
  % Step 4: Compute SVD of Rl*Rr
  %
  [u,s,v] = svd(Rl*Rr,"econ");
  %
  if(gamma>0)
    %
    s = diag(s);
    %
    denom = s.^2 + gamma;
    %
    s = diag(s./denom);
    %
  else
    %
    s = pinv(s);
    %
  end
  %
  % Step 5: Multiply the last core of Gl by u
  %
  Gl{end} = tensorprod(Gl{end},u,3,1);
  %
  % Step 6: Multiply v by the first core of Gr 
  %
  Gr{1} = tensorprod(v,Gr{1},1,1);
  %
  % Step 7: Transpose Gl
  %
  Utt = cell2core(tt_tensor,Gl);
  %
  UTtt = permute(tt_reshape(Utt,[Utt.n',Utt.r(end)],tt_tol),[Utt.d+1,1:Utt.d],tt_tol);
  %
  UTtt = tt_reshape(UTtt,Utt.n',tt_tol,Utt.r(end));
  %
  % Step 8: Transpose Gr
  %
  Vtt = cell2core(tt_tensor,Gr);
  %
  VTtt = permute(tt_reshape(Vtt,[Vtt.r(1),Vtt.n'],tt_tol),[1+(1:Vtt.d),1],tt_tol);
  %
  VTtt = tt_reshape(VTtt,Vtt.n',tt_tol,1,Vtt.r(1));
  %
  % Step 9: Assemble pseudo-inverse of Xtt
  %
  Gu = core2cell(UTtt);
  Gv = core2cell(VTtt);
  %
  Gv{end} = tensorprod(Gv{end},s,3,1);
  %
  Xttpinv = cell2core(tt_tensor,[Gv;Gu]); 
  %
end