% ------------------------------------------------------------------------%
% Inputs arguments                                                        %
% ------------------------------------------------------------------------%
% x                : snapshots as a 3D standard tensor or TT              %
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
% eps_tt           : TT truncation                                        %
% ------------------------------------------------------------------------%
function tt = opinf(x,u,t,K,isA,isF,tint_order,dt,dtp,tp,gamma,eps_tt)
  %
  tic;
  %
  if(isa(x, 'tt_tensor'))
    tt.Tpod = 0; % snapshots already calculated in tt format
  else
    %
    if(ndims(x)~=5)
      error("Snapshots must be a tensor X(i,j,k,q,n)\ni: x-index\nj: y-index\nk: z-index\nq: equation index\nn: time-index");
    end
    %
    t1=toc;
    x = tt_tensor(x,eps_tt);
    t2=toc;
    tt.Tpod = t2-t1;
    %
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
    error("Non-empty u is not implemented")
    %
  else
    %
    u = @(tt) 0;
    U = [];
    %
  end  
  %
  %% adjust max time step value for x
  %
  G      = core2cell(x);
  G{end} = G{end}(:,1:j2,:);
  x      = cell2core(tt_tensor,G);
  %
  %% learn operators
  %
  t1=toc;
  [A,F,B,D,R] = inferOp(x,j1,j2,U,tint_order,dt,isA,isF,gamma,eps_tt);
  t2=toc;
  %
  %% perform prediction
  %
  t3 = toc;
  xp = predict(x(:,:,:,:,j1-1),u,A,F,B,tp,dtp,tint_order,isA,isF,eps_tt);
  t4 = toc;
  %
  %% prepare output
  %
  tt.j1       = j1;
  tt.j2       = j2;
  tt.xp       = xp;
  tt.A        = A;
  tt.F        = F;
  tt.B        = B;
  tt.D        = D;
  tt.R        = R;
  tt.Tlearn   = t2-t1;
  tt.Tpredict = t4-t3;
  %
end