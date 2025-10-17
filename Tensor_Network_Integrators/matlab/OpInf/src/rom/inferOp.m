function [A,F,B,D,R]=inferOp(x,j1,j2,U,V,tint_order,dt,isA,isF,isRegularization,g1,g2)
  %
  K = length(j1:j2);
  n = size(V,2);
  s = n*(n+1)/2;
  %
  p = size(U,2);
  %
  % Project data to the POD space
  %
  R = calculateTimeDerivative(x,K,n,p,s,j1,j2,isA,isF,tint_order,dt,g1,g2);
  %
  D = calculateDataMatrix(x,U,K,n,p,s,j1,j2,isA,isF,isRegularization,g1,g2);
  %
  O = solveLS(D,R);
  %
  % Set outputs
  %
  i1 = 1;
  %
  % the diffusion term
  %
  if(isA)
    i2 = i1 + n - 1; 
    A  = O(:,i1:i2);
    i1 = i2 + 1;
  else
    A = [];
  end
  %
  % the source term
  %
  if(p>0)
    i2 = i1 + p - 1;
    B  = O(:,i1:i2);
    i1 = i2 + 1;
  else
    B = [];
  end
  %
  % the convection term
  %
  i2 = i1 + s - 1;
  if(isF)
    F = O(:,i1:i2);
  else
    F = [];
  end
  % 
end
%