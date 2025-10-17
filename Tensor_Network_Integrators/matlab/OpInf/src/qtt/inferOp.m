function [A,F,B,D,R]=inferOp(x,j1,j2,U,tint_order,dt,isA,isF,gamma,eps_tt,isQtt)
  %
  K = length(j1:j2);
  %
  p = size(U,2);
  %
  %fprintf("Calculating time derivatives\n")
  %
  R = calculateTimeDerivative(x,K,j1,j2,tint_order,dt,eps_tt,isQtt);
  %
  %fprintf("Calculating data matrix\n")
  %
  [D,sp_dim,t_dim]  = calculateDataMatrix(x,U,K,j1,j2,isA,isF,eps_tt,isQtt);
  %
  %fprintf("Solving LS\n")
  %
  O = solveLS(D,R,sp_dim,t_dim,gamma,eps_tt);
  %
  %fprintf("Assembling outputs A,F\n")
  %
  % Set outputs
  %
  G = core2cell(O);
  %
  i1 = 1;
  %
  % the diffusion term
  %
  if(isA)
    A = cell(3*(sp_dim+1),1);
    %
    A(1:sp_dim+1)            = G(1:sp_dim+1);
    A(sp_dim+2:2*(sp_dim+1)) = G(sp_dim+2:2*(sp_dim+1));
    %
    for k=1:sp_dim+1
      A{2*sp_dim+2+k} = G{2*sp_dim+2+k}(:,1,:);
    end
    %
    A  = cell2core(tt_tensor,A);
    %
    A  = tt_reshape(A,A.n(1:2*(sp_dim+1)));
    %
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
    F = cell(3*(sp_dim+1),1);
    %
    F(1:sp_dim+1)            = G(1:sp_dim+1);
    F(sp_dim+2:2*(sp_dim+1)) = G(sp_dim+2:2*(sp_dim+1));
    %
    for k=1:sp_dim+1
      F{2*sp_dim+2+k} = G{2*sp_dim+2+k}(:,i1:end,:);
    end

    F  = cell2core(tt_tensor,F);
    %
    F = round(F,eps_tt);
    %
  else
    F = [];
  end
  %
end
%