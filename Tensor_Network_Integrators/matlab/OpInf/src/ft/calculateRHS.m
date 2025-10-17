function RHS=calculateRHS(t,x,u,A,F,B,isA,isF,Nx,Ny,Nz,Neq)
  %
  %fprintf("t=%e\n",t);
  x = squeeze(reshape(x,[Nx,Ny,Nz,Neq]));
  %
  d = ndims(x);
  %
  RHS = zeros(size(x));
  %
  if(~isempty(B))
    error("u(t) is not yet implemented")
  end 
  %
  if(isA)
    %
    Abeg = d+1;
    Aend = ndims(A);
    %
    RHS(:,:,:,:) = RHS + tensorprod(A,x,Abeg:Aend,1:d); 
    %
  end
  %
  if(isF)
    %
    x2 = tensorprod(x,x);
    %
    Fbeg  = d+1;
    Fend  = ndims(F);
    %
    RHS(:,:,:,:) = RHS + tensorprod(F,x2,Fbeg:Fend,1:2*d); 
    %
  end
  %
  RHS = RHS(:);
  %
end
%