function O = solveLS(D,R)
  %
  % Solve the LS problem |Ax-b| using SVD: x = V(S^-1)(U^T)b 
  %
  [u,s,v] = svd(D,"econ");
  %
  O = (v*((u'*R)./diag(s)))';
  %
end