function [A,F,B,D,R]=inferOp(x,j1,j2,U,tint_order,dt,isA,isF,gamma,eps_tt)
  %
  K = length(j1:j2);
  %
  p = size(U,2);
  %
  %fprintf("Calculating time derivatives\n")
  %
  R = calculateTimeDerivative(x,K,j1,j2,tint_order,dt);
  %
  D = calculateDataMatrix(x,U,K,j1,j2,isA,isF,eps_tt);
  %
  O = solveLS(D,R,gamma,eps_tt);
  %
  % Set outputs
  %
  G = core2cell(O);
  %
  Nx  = size(D,1);
  Ny  = size(D,2);
  Nz  = size(D,3);
  Neq = size(D,4);
  %
  i1 = 1;
  %
  % the diffusion term
  %
  if(isA)
    A  = [G(1:8);G{9}(:,1,:);G{10}(:,1,:);G{11}(:,1,:);G{12}(:,1,:)];
    A  = cell2core(tt_tensor,A);
    A  = tt_reshape(A,[Nx,Ny,Nz,Neq,Nx,Ny,Nz,Neq]);
    A  = round(A,eps_tt);
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
    F  = [G(1:8);G{9}(:,i1:end,:);G{10}(:,i1:end,:);G{11}(:,i1:end,:);G{12}(:,i1:end,:)];
    F  = cell2core(tt_tensor,F);
    %
    F = tt_reshape(F,[Nx,Ny,Nz,Neq,Nx,Ny,Nz,Neq,Nx,Ny,Nz,Neq]);
    F = round(F,eps_tt);
    %
  else
    F = [];
  end
end
%