% ------------------------------------------------------------------------%
% Inputs arguments                                                        %
% ------------------------------------------------------------------------%
% x                : snapshots                                            %
% u                : function handle for time dependent source terms      %
% t                : time levels for the snapshots                        %
% K                : number of snapshots                                  %
% param            : energy retained or the number of the POD basis       %
% isA              : flag to determine if the diffusion operator is used  %
% isF              : flag to determine if the convection operator is used %
% tint_order       : time derivative and implicit time integration order  %
% dt               : time step for snapshots                              %
% dtp              : timestep for prediction                              %
% tp               : final time for prediction                            %
% g1               : regularization parameter for A and B                 %
% g2               : regularization parameter for F                       %
% ------------------------------------------------------------------------%
function rom = opinf(x,u,t,K,param,isA,isF,...
                     tint_order,dt,dtp,tp,g1,g2)
  %
  if(g1>0 || g2>0)
    isRegularization = true;
  else
    isRegularization = false;
  end
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
  else
    %
    u = @(tt) 0;
    U = [];
    %
  end  
  %
  % calculate POD basis
  %
  tic;
  t01=toc;
  [V,n,retained_energy] = calculatePOD(x(:,1:j2),param); % setting j from 1 to j2 is to ensure consistency with ttrom
  Xh = x(:,1:j2)'*V; % to consistenly compare to ttrom
  t02=toc;
  rom.Tpod = t02-t01;
  %
  % learn operators
  %
  t1=toc;
  [A,F,B,D,R] = inferOp(Xh,j1,j2,U,V,tint_order,dt,isA,isF,isRegularization,g1,g2);
  t2=toc;
  %
  % prepare initial conditions for prediction time stepping
  %
  xphat = V'*x(:,1:j1-1);
  %
  % perform prediction
  %
  m = n;
  %
  t3 = toc;
  [xphat,t_ode] = predict(xphat,u,A,F,B,m,tp,dtp,tint_order,isA,isF);
  t4 = toc;
  %
  xp = V(:,1:m)*xphat(:,end);
  %
  % prepare output
  %
  rom.j1              = j1;
  rom.j2              = j2;
  rom.xp              = xp;
  rom.xphat           = xphat';
  rom.t_ode           = t_ode;
  rom.Xh              = Xh;
  rom.n               = n;
  rom.V               = V;
  rom.A               = A;
  rom.F               = F;
  rom.B               = B;
  rom.D               = D;
  rom.R               = R;
  rom.retained_energy = retained_energy;
  rom.Tlearn          = t2-t1;
  rom.Tpredict        = t4-t3;
  %
  t5=toc;
  %
  rom.Telapsed = t5-t01;
  %
end