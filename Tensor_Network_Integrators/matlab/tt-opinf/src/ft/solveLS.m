function O = solveLS(D,R,gamma)
  %
  [Nx1,Ny1,Nz1,Neq1,Nx2,Ny2,Nz2,Neq2,K] = size(D);
  %
  N1 = Nx1*Ny1*Nz1*Neq1;
  N2 = Nx2*Ny2*Nz2*Neq2;
  %
  D = reshape(D,[N1*N2,K]);
  %
  % Solve the LS problem |Ax-b| using SVD: x = V(S^-1)(U^T)b 
  %
  [u,s,v] = svd(D,"econ");
  %
  u = u';
  %
  u = reshape(u,[K,Nx1,Ny1,Nz1,Neq1,Nx2,Ny2,Nz2,Neq2]);
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
  O = tensorprod(R,v,5,1);
  %
  O = tensorprod(O,s,5,1);
  %
  O = tensorprod(O,u,5,1);
  %
end