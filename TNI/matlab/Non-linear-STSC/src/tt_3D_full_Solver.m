function [U1,k,Ctt] =  tt_3D_full_Solver(testcase,SN,tol,eps,a,b)

Ku = testcase.Ku;
dKudu = testcase.dKudu;
bfun = testcase.bfun;
dbfun = testcase.dbfun;

rhsfn = testcase.rhsfn;
Fu = testcase.Fu;
dFudu = testcase.dFudu;
gfn = testcase.gfn;

%%
epsk = 1e-2;
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

%% create grid in tt for t,x,y,z
Itt = {ones(1,N1),ones(1,N1),ones(1,N1),ones(1,N1)};
Ctt{1} = cell2core(tt_tensor,{Tvec',ones(1,N1),ones(1,N1),ones(1,N1)});
for ic = 2:4
  temp = Itt;
  temp{ic} = XPt';
  Ctt{ic} = cell2core(tt_tensor,temp);
end

%% compute gBC
% compute BC over the whole domain
G = amen_cross(Ctt, @(x) cross_fun_nD(x,gfn),tol,'verb',0);
Gint = tt_set_zero_boundaries(G,{1,[1,N1],[1,N1],[1,N1]});
Gbc = round(G- Gint,tol);

%% build du/dt operatorsin TT format
It = 2:N1;
Ix = 2:N1-1;
% du/dt operator
Btt = matrices_to_tt_matrix_fn({At(It,It),I(Ix,Ix),I(Ix,Ix),I(Ix,Ix)});

Bmaptt = matrices_to_tt_matrix_fn({At(It,:),I(Ix,:),I(Ix,:),I(Ix,:)});

%% Laplace Operator
Laptt = matrices_to_tt_matrix_fn({I(It,It),Lapmat(Ix,Ix),I(Ix,Ix),I(Ix,Ix)}) ...
  + matrices_to_tt_matrix_fn({I(It,It),I(Ix,Ix),Lapmat(Ix,Ix),I(Ix,Ix)})...
  + matrices_to_tt_matrix_fn({I(It,It),I(Ix,Ix),I(Ix,Ix),Lapmat(Ix,Ix)});

Lapmaptt = matrices_to_tt_matrix_fn({I(It,:),Lapmat(Ix,:),I(Ix,:),I(Ix,:)}) ...
  + matrices_to_tt_matrix_fn({I(It,:),I(Ix,:),Lapmat(Ix,:),I(Ix,:)})...
  + matrices_to_tt_matrix_fn({I(It,:),I(Ix,:),I(Ix,:),Lapmat(Ix,:)});


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

%% compute F interior
g1tt = amen_cross(Ctt, @(x) cross_fun_nD(x,rhsfn),tol,'verb',0);
%get the inner of ftt
g1tt = tt_get_inner(g1tt,{2:N1,2:N1-1,2:N1-1,2:N1-1});

%% Check for constant function
for i = 1:3
  ucheck = tt_zeros(g1tt.n,numel(g1tt.n));
  if ~isa(bfun{i}(ucheck),'tt_tensor')
    bfun{i} = @(u,t,x,y,z) bfun{i}(u)*tt_ones(ucheck.n);
  end
  if ~isa(dbfun{i}(ucheck),'tt_tensor')
    dbfun{i} = @(u,t,x,y,z) dbfun{i}(u)*tt_ones(ucheck.n);
  end
end

%% Solution of the nonlinear system....Newton method
F = @(y) round(amen_mv(Btt,y,tol,'verb',0) - ...
  round(round(Ku(y),tol).*amen_mv(Laptt,y,tol,'verb',0),tol) ...
  + round(round(bfun{1}(y),tol).*amen_mv(Convtt1,y,tol,'verb',0),tol) ...
  + round(round(bfun{2}(y),tol).*amen_mv(Convtt2,y,tol,'verb',0),tol) ...
  + round(round(bfun{3}(y),tol).*amen_mv(Convtt3,y,tol,'verb',0),tol) ...
  - round(Fu(y),tol) - g1tt,tol) ...
  + Fbc(y, Bmaptt, Ku, Lapmaptt, bfun, Convttmap1, Convttmap2, Convttmap3, Gbc, tol);

% Initial guess
U0 = tt_zeros(g1tt.n,numel(g1tt.n));
% U0 = exacttt + 1e-2*tt_rand(g1tt.n,numel(g1tt.n),1); % Initial guess

%%
nrmFU0 = norm(F(U0));
maxiter = 200;

for k = 1:maxiter
  tsiter = datetime;
  % epsk = epsk*epsfactor;
  fprintf('\n******** Newton iter = %d ************\n\n', k);
  fprintf('epsk = %.5e \n', epsk);

  % compute Jacobian with JFbc
  Jmatk = Jfuntt(Btt,Laptt,Ku, dKudu, bfun, dbfun, Convtt1, Convtt2, Convtt3,dFudu,U0,tol) +...
    JFbc(U0, dKudu, Lapmaptt, bfun, Convttmap1, Convttmap2, Convttmap3, Gbc, tol);

  % Jmatk = Jfuntt(Btt,Laptt,Ku,dKudu,dFudu,U0,tol);

  % use GMRES to solve for du
  b = -F(U0);

  du = amen_solve2(Jmatk,b,epsk,'nswp',30,'verb',1);

  tegmres = datetime;
  fprintf('Linear Solve Time = %.2f\n', seconds(tegmres-tsiter));

  %update U1
  for alpha = [1,1/2,1/4,1/8,1/16]
    U1 = round(U0 + alpha*du, max(epsk,tol));
    if norm(F(U1))<norm(F(U0))
      break;
    end
  end

  %compute local error
  local_err=norm(U1-U0)/norm(U0);
  res_norm = norm(F(U1))/nrmFU0;

  % #* set epsk
  epsk = min(min([local_err,res_norm]),epsk);

  fprintf('apha = %.2e \n', alpha);
  fprintf('u_err = %.5e,  Fu_ratio = %.5e \n', local_err, res_norm)
  fprintf('iter time = %.2f \n', seconds(datetime-tsiter));
  if (local_err<eps) || (res_norm<eps)
    break;
  end
  U0 = U1;
end
end

%% subroutines
function fval = Fbc(y, Bmap, Ku, Lapmap, bfun, Convmap1, Convmap2, Convmap3, Gbc, tol)
diagKu = make_tt_to_operator(round(Ku(y),tol));

diagbf1 = make_tt_to_operator(round(bfun{1}(y),tol));
diagbf2 = make_tt_to_operator(round(bfun{2}(y),tol));
diagbf3 = make_tt_to_operator(round(bfun{3}(y),tol));

Amap = round(Bmap - amen_mm(diagKu,Lapmap,tol,'verb',0) ...
  + amen_mm(diagbf1,Convmap1,tol,'verb',0)...
  + amen_mm(diagbf2,Convmap2,tol,'verb',0)...
  + amen_mm(diagbf3,Convmap3,tol,'verb',0),tol);
fval = amen_mv(Amap,Gbc,tol,'verb',0);
end

function fval = JFbc(y, dKudu, Lapmap, bfun, Convmap1, Convmap2, Convmap3, Gbc, tol)
fval = ...
make_tt_to_operator(round((-1)*dKudu(y).*amen_mv(Lapmap,Gbc,tol,'verb',0),tol))...
+ make_tt_to_operator(round(bfun{1}(y).*amen_mv(Convmap1,Gbc,tol,'verb',0),tol))...
+ make_tt_to_operator(round(bfun{2}(y).*amen_mv(Convmap2,Gbc,tol,'verb',0),tol))...
+ make_tt_to_operator(round(bfun{3}(y).*amen_mv(Convmap3,Gbc,tol,'verb',0),tol));
end

function Jmat = Jfuntt(B,Lap,Ku,dKudu,bfun,dbfun,Conv1,Conv2,Conv3,dFudu,y,tol)
%FG-Jfun = @(y) B - (diag(y.^2)*Lap + 2*diag(y.*(Lap*y)) - 3*diag(y.^2));
% Jfun = @(y) B - (diag(Ku(y))*Lap + diag((dKudu(y)).*(Lap*y)) - diag(dFudu(y)));
diagKu = make_tt_to_operator(round(Ku(y),tol));
diagdKuduLy = make_tt_to_operator(round(dKudu(y).*amen_mv(Lap,y,tol,'verb',0),tol));

diagbf1 = make_tt_to_operator(round(bfun{1}(y),tol));
diagbf2 = make_tt_to_operator(round(bfun{2}(y),tol));
diagbf3 = make_tt_to_operator(round(bfun{3}(y),tol));

diagdbf1 = make_tt_to_operator(round(dbfun{1}(y).*amen_mv(Conv1,y,tol,'verb',0),tol));
diagdbf2 = make_tt_to_operator(round(dbfun{2}(y).*amen_mv(Conv2,y,tol,'verb',0),tol));
diagdbf3 = make_tt_to_operator(round(dbfun{3}(y).*amen_mv(Conv3,y,tol,'verb',0),tol));

diagdFudu = make_tt_to_operator(round(dFudu(y),tol));
Jmat = round(B ...
  - (amen_mm(diagKu,Lap,tol,'verb',0) + diagdKuduLy) ...
  - (amen_mm(diagbf1,Conv1,tol,'verb',0) + diagdbf1) ...
  - (amen_mm(diagbf2,Conv2,tol,'verb',0) + diagdbf2) ...
  - (amen_mm(diagbf3,Conv3,tol,'verb',0) + diagdbf3) ...
  - diagdFudu,tol);

end

