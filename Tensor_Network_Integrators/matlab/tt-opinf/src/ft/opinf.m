% ------------------------------------------------------------------------%
% Inputs arguments                                                        %
% ------------------------------------------------------------------------%
% x                : snapshots                                            %
% u                : function handle for time dependent source terms      %
% t                : time levels for the snapshots                        %
% K                : number of snapshots                                  %
% isA              : flag to determine if the diffusion operator is used  %
% isF              : flag to determine if the convection operator is used %
% tint_order       : time derivative and implicit time integration order  %
% dt               : time step for snapshots                              %
% dtp              : timestep for prediction                              %
% tp               : final time for prediction                            %
% gamma            : regularization parameter                             %
% ------------------------------------------------------------------------%
function ft = opinf(x,u,t,K,isA,isF,tint_order,dt,dtp,tp,gamma)
  %
  % determine training timestep interval
  %
  j1 = 1 + tint_order;
  j2 = j1 + K - 1;
  %
  % prepare the source term
  %
  if(~isempty(u))
    %
    U = zeros(K,length(u(1)));
    %
    for j=1:K
      U(j,:) = u(t(j-1+j1));
    end
    %
    error("Non-empty u is not implemented")
    %
  else
    %
    u = @(tt) 0;
    U = [];
    %
  end  
  %
  % learn operators
  %
  tic;
  t1 = toc;
  [A,F,B,D,R] = inferOp(x(:,:,:,:,1:j2),j1,j2,U,tint_order,dt,isA,isF,gamma);
  t2 = toc;
  %
  % perform prediction
  %
  [Nx,Ny,Nz,Neq,~] = size(D);
  t3 = toc;
  xp = predict(reshape(x(:,:,:,:,j1-1),[Nx,Ny,Nz,Neq]),u,A,F,B,tp,dtp,tint_order,isA,isF);
  t4 = toc;
  %
  % prepare output
  %
  ft.j1       = j1;
  ft.j2       = j2;
  ft.xp       = xp;
  ft.A        = A;
  ft.F        = F;
  ft.B        = B;
  ft.D        = D;
  ft.R        = R;
  ft.Tlearn   = t2-t1;
  ft.Tpredict = t4-t3;
  %
end