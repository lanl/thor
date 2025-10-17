function [D,sp_dim,t_dim] = calculateDataMatrix(x,U,K,j1,j2,isA,isF,eps_tt,isQtt)
  %
  G      = core2cell(x);
  G{end} = G{end}(:,j1:j2,:);
  x      = cell2core(tt_tensor,G);
  %
  sp_dim = numel(x.n(1:end-2)); 
  %
  if(isQtt)
    %
    t_dim  = numel(factor(x.n(end)));
    %
    n = [x.n(1:end-2)' x.n(end-1) factor(x.n(end))];
    %
    x = tt_reshape(x,n);
    %
    x = round(x,eps_tt);
    %
  else
    %
    t_dim = 1;
    %
  end
  %
  % Calculate column size of the data matrix for linear, source and quadratic terms
  %
  p = 0;
  %
  Dsz = zeros(1,sp_dim+1+numel(x.n));
  %
  Dsz(1:sp_dim+1)      = x.n(1:sp_dim+1)';
  Dsz(end-t_dim+1:end) = x.n(end-t_dim+1:end)';
  %
  if(isA)
    Dsz(sp_dim+2:2*sp_dim+2) = Dsz(sp_dim+2:2*sp_dim+2) + 1;
  end
  %
  if(~isempty(U))
    error("u(t) is not supported yet")
  end
  %
  if(isF) 
    Dsz(sp_dim+2:2*sp_dim+2) =  Dsz(sp_dim+2:2*sp_dim+2) + x.n(1:sp_dim+1)';
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
    xtemp1 = tt_reshape(x,[x.n(1:sp_dim+1)',ones(1,sp_dim+1),x.n(end-t_dim+1:end)']);
    xtemp2 = tt_reshape(x,[ones(1,sp_dim+1),x.n(1:sp_dim+1)',x.n(end-t_dim+1:end)']);
    %
    xtemp1 = tt_broadcast(xtemp1, sp_dim+2:2*sp_dim+2,x.n(1:sp_dim+1)');
    xtemp2 = tt_broadcast(xtemp2, 1:sp_dim+1,         x.n(1:sp_dim+1)');
    %
    x2 = cross_interpolation({xtemp1 xtemp2},@(xxx) myfun(xxx),eps_tt,[],0);
    %
  end
  %
  i1 = 1;
  %
  % Contribution of the diffusion term
  %
  if(isA)
    %
    x = tt_reshape(x,[x.n(1:sp_dim+1)',ones(1,sp_dim+1),x.n(end-t_dim+1:end)']);
    %
    G = core2cell(x); 
    %
    for idx=1:sp_dim+1
      Gtemp           = zeros(x.r(sp_dim+1+idx),Dsz(sp_dim+1+idx),x.r(sp_dim+2+idx)); 
      Gtemp(:,1,:)    = G{sp_dim+1+idx};
      G{sp_dim+1+idx} = Gtemp;
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
    for idx=1:sp_dim+1
      Gtemp             = zeros(x2.r(sp_dim+1+idx),Dsz(sp_dim+1+idx),x2.r(sp_dim+2+idx)); 
      Gtemp(:,i1:end,:) = G{sp_dim+1+idx};
      G{sp_dim+1+idx}   = Gtemp;
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