function [U1,k,Ctt] =  FG_3D_full_Solver(testcase,SN,tol,eps,a,b)

Ku = testcase.Ku;
dKudu = testcase.dKudu;
bfun = testcase.bfun;
dbfun = testcase.dbfun;

rhsfn = testcase.rhsfn;
Fu = testcase.Fu;
dFudu = testcase.dFudu;
gfn = testcase.gfn; %boundary function

%%
N1 = SN+1;
type = 'Chebyshev';

%%
% space matrix
[tempDMA,temppts] = compute_derivative_mat_fn(SN,type);

XPt = (b-a)./2.*temppts + (a+b)/2;
DMA = (2/(b-a))*tempDMA;

Lapmat = DMA*DMA; %Laplace matrix

%time matrix
At = 2*tempDMA;
Tvec = (temppts+1)*0.5;
I = eye(N1,N1);

%% create grid in tt
Itt = {ones(1,N1),ones(1,N1),ones(1,N1),ones(1,N1)};
Ctt{1} = cell2core(tt_tensor,{Tvec',ones(1,N1),ones(1,N1),ones(1,N1)});
for ic = 2:4
  temp = Itt;
  temp{ic} = XPt';
  Ctt{ic} = cell2core(tt_tensor,temp);
end

%% build du/dt operatorsin TT format
It = 2:N1;
Ix = 2:N1-1;
% du/dt operator
Btt = matrices_to_tt_matrix_fn({At(It,It),I(Ix,Ix),I(Ix,Ix),I(Ix,Ix)});
Bmaptt = matrices_to_tt_matrix_fn({At(It,:),I(Ix,:),I(Ix,:),I(Ix,:)});
B = full(Btt);
Bmap = full(Bmaptt);

%% Laplace Operator
Laptt = matrices_to_tt_matrix_fn({I(It,It),Lapmat(Ix,Ix),I(Ix,Ix),I(Ix,Ix)}) ...
  + matrices_to_tt_matrix_fn({I(It,It),I(Ix,Ix),Lapmat(Ix,Ix),I(Ix,Ix)})...
  + matrices_to_tt_matrix_fn({I(It,It),I(Ix,Ix),I(Ix,Ix),Lapmat(Ix,Ix)});
Lapmaptt = matrices_to_tt_matrix_fn({I(It,:),Lapmat(Ix,:),I(Ix,:),I(Ix,:)}) ...
  + matrices_to_tt_matrix_fn({I(It,:),I(Ix,:),Lapmat(Ix,:),I(Ix,:)})...
  + matrices_to_tt_matrix_fn({I(It,:),I(Ix,:),I(Ix,:),Lapmat(Ix,:)});

Lap = full(Laptt);
Lapmap = full(Lapmaptt);

%% Convection operator

if isa(bfun{1},'function_handle')

  Convecmat = DMA;
  Convtt1 = matrices_to_tt_matrix_fn({I(It,It),Convecmat(Ix,Ix),I(Ix,Ix),I(Ix,Ix)});
  Convtt2 = matrices_to_tt_matrix_fn({I(It,It),I(Ix,Ix),Convecmat(Ix,Ix),I(Ix,Ix)});
  Convtt3 = matrices_to_tt_matrix_fn({I(It,It),I(Ix,Ix),I(Ix,Ix), Convecmat(Ix,Ix)});

  %
  Convttmap1 = matrices_to_tt_matrix_fn({I(It,:),Convecmat(Ix,:),I(Ix,:),I(Ix,:)});
  Convttmap2 = matrices_to_tt_matrix_fn({I(It,:),I(Ix,:),Convecmat(Ix,:),I(Ix,:)});
  Convttmap3 = matrices_to_tt_matrix_fn({I(It,:),I(Ix,:),I(Ix,:), Convecmat(Ix,:)});

end

%% convert to full matrix
Conv1 = full(Convtt1);
Conv2 = full(Convtt2);
Conv3 = full(Convtt3);
Convmap1 = full(Convttmap1);
Convmap2 = full(Convttmap2);
Convmap3 = full(Convttmap3);

%% compute F interior
g1tt = amen_cross(Ctt, @(x) cross_fun_nD(x,rhsfn),tol,'verb',0);
%get the inner of ftt
g1tt = tt_get_inner(g1tt,{2:N1,2:N1-1,2:N1-1,2:N1-1});
g1 = full(g1tt);

%% compute gBC
% compute BC over the whole domain
G = amen_cross(Ctt, @(x) cross_fun_nD(x,gfn),tol,'verb',0);
Gint = tt_set_zero_boundaries(G,{1,[1,N1],[1,N1],[1,N1]});
Gbc = full(round(G- Gint,tol));

%% Solution of the nonlinear system....Newton method

Fbc = @(y) (Bmap-diag(Ku(y))*Lapmap + ...
  diag(bfun{1}(y))*Convmap1 +diag(bfun{2}(y))*Convmap2 ...
  + diag(bfun{3}(y))*Convmap3 ...
  )*Gbc; %boundary function


F = @(y) B*y - (Ku(y)).*(Lap*y) + diag(bfun{1}(y))*(Conv1*y) ...
  +diag(bfun{2}(y))*(Conv2*y) + diag(bfun{3}(y))*(Conv3*y)...
  - Fu(y) - g1 + Fbc(y);

%%%% Jacobian
JFbc = @(y) (-1)*diag(dKudu(y).*(Lapmap*Gbc)) ...
             + diag(bfun{1}(y).*(Convmap1*Gbc)) ...
             + diag(bfun{2}(y).*(Convmap2*Gbc))...
             + diag(bfun{3}(y).*(Convmap3*Gbc));

Jfun = @(y) B - (diag(Ku(y))*Lap + diag((dKudu(y)).*(Lap*y)))...
  + (diag(bfun{1}(y))*Conv1 + diag((dbfun{1}(y)).*(Conv1*y)))...
  + (diag(bfun{2}(y))*Conv2 + diag((dbfun{2}(y)).*(Conv2*y)))...
  + (diag(bfun{3}(y))*Conv3 + diag((dbfun{3}(y)).*(Conv3*y)))...
 - diag(dFudu(y)) + JFbc(y);


%%%%%%
U0 = zeros(numel(g1),1); % Initial guess
% U0 = exacttt + 5e-1*rand(size(g1));
nrmFU0 = norm(F(U0));
maxiter = 1000;

for k = 1:maxiter
  tsiter = datetime;
  fprintf('\n******** Newton iter = %d ************\n\n', k)
  
  %compute Jacobian with new u
  Jmatk = Jfun(U0);

  % use GMRES to solve for du
  b = -F(U0);
  % du = gmres(@(v) Jvfun(F,U0,1e-6,v), b, [], gmres_eps, 100,[],[]);
  du = Jmatk\b;
  % du = gmres(@(v) Jvfun(v), b, [], eps, 300,[],[]);
  %update U1
  for alpha = [1,1/2,1/4,1/8,1/16]
    U1 = U0 + alpha*du;
    if norm(F(U1))<norm(F(U0))
      break;
    end
  end

  % U1tt = tt_tensor(reshape(U1,N1-1,N1-2,N1-2),eps)
  %compute local error
  local_err=norm(U1-U0)/norm(U0);
  res_norm = norm(F(U1))/nrmFU0;
  fprintf('apha = %.2e \n', alpha);
  fprintf('u_err = %.5e,  Fu_ratio = %.5e \n', local_err, res_norm)
  fprintf('iter time = %.2f \n', seconds(datetime-tsiter));
  if (local_err<eps) || (res_norm<eps)
    break;
  end
  U0 = U1;

end
