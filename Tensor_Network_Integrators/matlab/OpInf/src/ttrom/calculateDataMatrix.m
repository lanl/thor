function D = calculateDataMatrix(Xh,U,K,n,p,s,j1,j2,isA,isF,isRegularization,g1,g2)
  %
  % Calculate column size of the data matrix for linear, source and quadratic terms
  %
  Dsz = 0;
  %
  if(isA)
    Dsz = Dsz + n;
  end
  %
  Dsz = Dsz + p; 
  %
  if(isF) 
    Dsz = Dsz + s;
  end  
  %
  % Calculate additional row size due to regularization
  %
  sz = 0;
  %
  if(g1>0)
    if(isA)
      sz = sz + n;
    end
    sz = sz + p;
  end
  %
  if(g2>0 && isF) 
    sz = sz + s;
  end
  %
  % Allocate memory
  %
  D = zeros(K+sz,Dsz);
  %
  % Get the diagonal indices corresponding to regularization terms
  %
  if(isRegularization)
    %
    idx = K+1;
    %
    if(g1<=0)
      if(isA)
        idx = idx + (K+sz)*(n+p);
      else
        idx = idx + (K+sz)*(p);
      end
    end 
    %
    idx = idx:K+sz+1:(K+sz)*Dsz;
    %
  end
  %
  i1 = 1;
  r1 = 1;
  %
  % Contribution of the diffusion term
  %
  if(isA)
    i2 = i1 + n - 1; 
    D(1:K,i1:i2) = Xh(j1:j2,:);
    %
    if(g1>0)
      %
      r2            = r1 + n - 1;
      D(idx(r1:r2)) = g1;
      r1            = r2 + 1;
      %
    end
    %
    i1 = i2 + 1;
    %
  end
  %
  % Contribution of the source term
  %
  if(p>0)
    i2           = i1 + p - 1;
    D(1:K,i1:i2) = U;
    i1           = i2 + 1;
  end
  %
  if(g1>0)
    r2            = r1 + p - 1;
    D(idx(r1:r2)) = g1;
    r1            = r2 + 1;
  end
  %
  %
  % Contribution of the convection term
  %
  if(isF)
    %
    X = Xh(j1:j2,:)'; % n x K
    %
    for i=1:n
      %
      i2 = i1 + i - 1;
      %
      Xi = zeros(i,K);
      %
      for k=1:K
        Xi(:,k) = X(i,k)*X(1:i,k);
      end
      %
      D(1:K,i1:i2) = Xi';
      %
      if(g2>0)
        r2            = r1 + i -1;
        D(idx(r1:r2)) = g2;
        r1            = r2 + 1;
      end
      %
      %
      i1 = i2 + 1;
      %
    end
    %
  end
  %
end