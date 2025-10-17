% ------------------------------------------------------------------------%
% Inputs arguments                                                        %
% ------------------------------------------------------------------------%
% x                : either snapshots or a struct to define input handling%
% u                : function handle for time dependent source terms      %
% t                : time levels for the snapshots                        %
% K                : number of snapshots                                  %
% isA              : flag to determine if the diffusion operator is used  %
% isF              : flag to determine if the convection operator is used %
% tint_order       : time derivative and implicit time integration order  %
% dt               : time step for snapshots                              %
% dtp              : timestep for prediction                              %
% tp               : final time for prediction                            %
% g1               : regularization parameter for A and B                 %
% g2               : regularization parameter for F                       %
% eps_tt           : TT truncation error                                  %
% ------------------------------------------------------------------------%
function rom = opinf(x,u,t,K,isA,isF,tint_order,dt,dtp,tp,g1,g2,eps_tt)
  %
  tic;
  t01=toc;
  %
  % determine training timestep interval
  %
  j1 = 1 + tint_order;
  j2 = j1 + K - 1;
  if(isstruct(x))
    cfg = x;
    [xtt,Tpod] = get_snapshots_from_file_with_cross(cfg,eps_tt,j2);
  else
    cfg = [];
  end
  %
  if(g1>0 || g2>0)
    isRegularization = true;
  else
    isRegularization = false;
  end
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
  if(~isstruct(cfg))
    temp = x(:,1:j2);
    xtt = tt_tensor(temp,eps_tt);
    %xtt = cross_interpolation(size(temp),@(xx)myfun(xx,temp),eps_tt,0,true);
    %xtt = tt_reshape(xtt,[factor(xtt.n(1:end-1)) xtt.n(end)]);
    %xtt = tt_orthogonolize(xtt);
    t02=toc;
    rom.Tpod = t02-t01;
  else
    rom.Tpod = Tpod;
  end
  %
  % learn operators
  %
  t1=toc;
  [A,F,B,D,R] = inferOp(xtt{numel(xtt.n)},j1,j2,U,tint_order,dt,isA,isF,isRegularization,g1,g2);
  t2=toc;
  %
  % prepare initial conditions for prediction time stepping
  %
  xphat = xtt{numel(xtt.n)}(:,1:j1-1);
  %
  % perform prediction
  %
  n = xtt.r(end-1);
  m = n;
  %
  t3 = toc;
  [xphat, t_ode] = predict(xphat,u,A,F,B,m,tp,dtp,tint_order,isA,isF);
  t4 = toc;
  %
  G               = core2cell(xtt);
  rom.G           = G;
  G{numel(xtt.n)} = xphat(:,end);
  xptt            = cell2core(tt_tensor,G);
  xp              = full(xptt); 
  %
  % prepare output
  %
  rom.j1              = j1;
  rom.j2              = j2;
  rom.xp              = xp;
  rom.xptt            = xptt;
  rom.xphat           = xphat';
  rom.t_ode           = t_ode;
  rom.n               = n;
  rom.retained_energy = 0;
  rom.A               = A;
  rom.F               = F;
  rom.B               = B;
  rom.D               = D;
  rom.R               = R;
  rom.Tlearn          = t2-t1;
  rom.Tpredict        = t4-t3;
  %
  t5=toc;
  %
  rom.Telapsed = t5-t01;
  %
end
%
function y=myfun(n,Q)
  dims = size(Q);
  j    = ones(size(n,1),1);
  mult = 1;
  for k = 1:size(n,2)
    j = j + (n(:,k)-1)*mult;
    mult = mult * dims(k);
  end
  y = Q(j);
end