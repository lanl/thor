function [A,F,B,D,R]=inferOp(x,j1,j2,U,tint_order,dt,isA,isF,gamma)
  %
  K = length(j1:j2);
  %
  p = size(U,2);
  %
  R = calculateTimeDerivative(x,K,j1,j2,tint_order,dt);
  %
  D = calculateDataMatrix(x,U,K,j1,j2,isA,isF);
  %
  O = solveLS(D,R,gamma);
  %
  % Set outputs
  %
  i1 = 1;
  %
  % the diffusion term
  %
  if(isA)
    A  = O(:,:,:,:,:,:,:,:,1,1,1,1);
    i1 = i1 + 1;
  else
    A = [];
  end
  %
  % the source term
  %
  if(p>0)
    %
    error("u(t) is not yet implemented");
    %
  else
    B = [];
  end
  %
  % the convection term
  %
  if(isF)
    %
    F = O(:,:,:,:,:,:,:,:,i1:end,i1:end,i1:end,i1:end);
    %
  else
    F = [];
  end
  % 
end
%