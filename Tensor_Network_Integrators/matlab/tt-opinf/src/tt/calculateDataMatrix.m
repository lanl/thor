function D = calculateDataMatrix(x,U,K,j1,j2,isA,isF,eps_tt)
  %
  G      = core2cell(x);
  G{end} = G{end}(:,j1:j2,:);
  x      = cell2core(tt_tensor,G);
  %
  % Calculate column size of the data matrix for linear, source and quadratic terms
  %
  Nx  = size(x,1);
  Ny  = size(x,2);
  Nz  = size(x,3);
  Neq = size(x,4);
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
  D = tt_zeros(Dsz);
  %
  % Contribution of the convection term first
  %
  if(isF)
    %
    xtemp1 = tt_reshape(x,[Nx,Ny,Nz,Neq,1,1,1,1,K]);
    xtemp2 = tt_reshape(x,[1,1,1,1,Nx,Ny,Nz,Neq,K]);
    %
    xtemp1 = tt_broadcast(xtemp1, 5:8,[Nx Ny Nz Neq]);
    xtemp2 = tt_broadcast(xtemp2, 1:4,[Nx Ny Nz Neq]);
    %
    x2 = cross_interpolation({xtemp1 xtemp2},@(xxx) myfun(xxx),eps_tt,[],0);
    %x2 = xtemp1.*xtemp2;
    %x2 = round(x2,eps_tt);
    %
  end
  %
  i1 = 1;
  %
  % Contribution of the diffusion term
  %
  if(isA)
    %
    x = tt_reshape(x,[Nx,Ny,Nz,Neq,1,1,1,1,K]);
    %
    G = core2cell(x); 
    %
    for idx=1:4
      Gtemp        = zeros(x.r(4+idx),Dsz(4+idx),x.r(5+idx)); 
      Gtemp(:,1,:) = G{4+idx};
      G{4+idx}     = Gtemp;
    end
    %
    x = cell2core(tt_tensor,G);
    %
    D = round(D + x, eps_tt);
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
    G = core2cell(x2); 
    %
    for idx=1:4
      Gtemp             = zeros(x2.r(4+idx),Dsz(4+idx),x2.r(5+idx)); 
      Gtemp(:,i1:end,:) = G{4+idx};
      G{4+idx}          = Gtemp;
    end
    %
    x2 = cell2core(tt_tensor,G);
    %
    D = round(D + x2, eps_tt);
    %
  end
  %
  D = round(D, eps_tt);
  %
end
%
function y = myfun(x)
  y = x(:,1).*x(:,2);
end
%