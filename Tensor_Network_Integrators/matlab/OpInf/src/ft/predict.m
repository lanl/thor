function x=predict(x,u,A,F,B,Tfinal,dt,tint_order,isA,isF)
  %
  [Nx,Ny,Nz,Neq] = size(x);
  %
  if(isA)
    A = squeeze(A);
  end
  %
  if(isF)
    F = squeeze(F);
  end
  %
  if(tint_order<3)
    [~, x] = ode23(@(tt,xx)calculateRHS(tt,xx,u,A,F,B,isA,isF,Nx,Ny,Nz,Neq), 0:dt:Tfinal, x);
  else
    [~, x] = ode45(@(tt,xx)calculateRHS(tt,xx,u,A,F,B,isA,isF,Nx,Ny,Nz,Neq), 0:dt:Tfinal, x);
  end
  %
  x = reshape(x(end,:),[Nx,Ny,Nz,Neq]);
  %
end