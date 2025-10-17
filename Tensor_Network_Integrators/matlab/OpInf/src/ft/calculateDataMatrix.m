function D = calculateDataMatrix(x,U,K,j1,j2,isA,isF)
  %
  % Calculate column size of the data matrix for linear, source and quadratic terms
  %
  [Nx,Ny,Nz,Neq,~] = size(x);
  %
  p = 0;
  %
  Dsz = zeros(1,9);
  %
  Dsz(1:4) = [Nx Ny Nz Neq];
  Dsz(end) = K;
  %
  if(isA)
    Dsz(5:end-1) = Dsz(5:end-1) + 1;
  end
  %
  if(~isempty(U))
    error("u(t) is not supported yet")
  end
  %
  if(isF) 
    Dsz(5:end-1) =  Dsz(5:end-1) + [Nx Ny Nz Neq];
  end  
  %
  % Allocate memory
  %
  D = zeros(Dsz);
  %
  i1 = 1;
  %
  % Contribution of the diffusion term
  %
  if(isA)
    %
    D(:,:,:,:,1,1,1,1,1:K) = reshape(x(:,:,:,:,j1:j2),[Nx Ny Nz Neq 1 1 1 1 K]);
    %
    i1 = i1 + 1;
    %
  end
  %
  % Contribution of the source term
  %
  if(p>0)
    error("u(t) is not yet implemented")
  end
  %
  % Contribution of the convection term
  %
  if(isF)
    %
    for i=1:K
      %temp1 = reshape(x(:,:,:,:,j1:j2), [Nx,Ny,Nz,Neq,1,1,1,1,K]);
      %temp2 = reshape(x(:,:,:,:,j1:j2), [1,1,1,1,Nx,Ny,Nz,Neq,K]);
      %D(:,:,:,:,i1:end,i1:end,i1:end,i1:end,:) = temp1.*temp2;
      %
      D(:,:,:,:,i1:end,i1:end,i1:end,i1:end,i) = tensorprod(squeeze(x(:,:,:,:,j1+i-1)),squeeze(x(:,:,:,:,j1+i-1)));
      %
    end
    %
  end
  %
end