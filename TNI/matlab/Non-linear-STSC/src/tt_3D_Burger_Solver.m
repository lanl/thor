function [U1,k,Ctt] =  tt_3D_Burger_Solver(testcase,SN,tol,eps,a,b)

bfun = testcase.bfun;
dbfun = testcase.dbfun;
gfn = testcase.gfn; %boundary function
%%
N1 = SN+1;
type = 'Chebyshev';

epsk = 1e-1;
%%
% space matrix
[tempDMA,temppts] = compute_derivative_mat_fn(SN,type);

XPt = (b-a)./2.*temppts + (a+b)/2;
DMA = (2/(b-a))*tempDMA;

Lapmat = DMA*DMA; %Laplace matrix

%time matrix
t0 = 0;
t1 = 1.0;
Tvec = (t1-t0)./2.*temppts + (t1+t0)/2;
At = (2/(t1-t0))*tempDMA;
% At = 2*tempDMA;
% Tvec = (temppts+1)*0.5;
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

%% compute gBC
% compute BC over the whole domain
G = amen_cross(Ctt, @(x) cross_fun_nD(x,gfn),tol,'verb',0);
Gint = tt_set_zero_boundaries(G,{1,[1,N1],[1,N1],[1,N1]});
Gbc = round(G- Gint,tol);

%% Solution of the nonlinear system....Newton method

% Fbc = @(y) (Bmap - Lapmap + ...
%   diag(bfun{1}(y))*Convmap1 +diag(bfun{2}(y))*Convmap2 ...
%   + diag(bfun{3}(y))*Convmap3 ...
%   )*Gbc; %boundary function

% F = @(y) B*y - (Lap*y) + diag(bfun{1}(y))*(Conv1*y) ...
%   + diag(bfun{2}(y))*(Conv2*y) + diag(bfun{3}(y))*(Conv3*y)...
%   + Fbc(y, Bmaptt, Lapmaptt, bfun, Convttmap1, Convttmap2, Convttmap3, Gbc, tol);

F = @(y) round(amen_mv(Btt,y,tol,'verb',0) ...
  - amen_mv(Laptt,y,tol,'verb',0)...
  + round(round(bfun{1}(y),tol).*amen_mv(Convtt1,y,tol,'verb',0),tol) ...
  + round(round(bfun{2}(y),tol).*amen_mv(Convtt2,y,tol,'verb',0),tol) ...
  + round(round(bfun{3}(y),tol).*amen_mv(Convtt3,y,tol,'verb',0),tol),tol) ...
  + Fbc(y, Bmaptt, Lapmaptt, bfun, Convttmap1, Convttmap2, Convttmap3, Gbc, tol);


% %%%% Jacobian
% JFbc = @(y)  diag(bfun{1}(y).*(Convmap1*Gbc)) ...
%   + diag(bfun{2}(y).*(Convmap2*Gbc))...
%   + diag(bfun{3}(y).*(Convmap3*Gbc));
% 
% Jfun = @(y) B - Lap ...
%   + (diag(bfun{1}(y))*Conv1 + diag((dbfun{1}(y)).*(Conv1*y)))...
%   + (diag(bfun{2}(y))*Conv2 + diag((dbfun{2}(y)).*(Conv2*y)))...
%   + (diag(bfun{3}(y))*Conv3 + diag((dbfun{3}(y)).*(Conv3*y)))...
%   + JFbc(y);

%%%%%%
% exacttt = amen_cross(Ctt, @(x) cross_fun_nD(x,testcase.exactfn),tol);
% exacttt = tt_get_inner(exacttt,{2:N1,2:N1-1,2:N1-1,2:N1-1});

U0 = tt_zeros([N1-1,N1-2,N1-2,N1-2],4);
% U0 = exacttt + 1e-5*rand(size(exacttt));

%%
nrmFU0 = norm(F(U0));
maxiter = 100;
du = U0;
for k = 1:maxiter
  tsiter = datetime;
  % epsk = epsk*epsfactor;
  fprintf('\n******** Newton iter = %d ************\n\n', k);
  fprintf('epsk = %.5e \n', epsk);

  % compute Jacobian with JFbc
  Jmatk = Jfuntt(Btt,Laptt,bfun, Convtt1, Convtt2, Convtt3,U0,tol) +...
    JFbc(U0, bfun, Convttmap1, Convttmap2, Convttmap3, Gbc, tol);

  % Jmatk = Jfuntt(Btt,Laptt,Ku,dKudu,dFudu,U0,tol);

  % use GMRES to solve for du
  b = -F(U0);
  
  % du0 = round(du,epsk);

  du = amen_solve2(Jmatk,b,epsk,'nswp',30,'verb',1,'x0');

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
  % epsk = min([local_err,res_norm]);

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
function fval = Fbc(y, Bmap, Lapmap, bfun, Convmap1, Convmap2, Convmap3, Gbc, tol)

diagbf1 = make_tt_to_operator(round(bfun{1}(y),tol));
diagbf2 = make_tt_to_operator(round(bfun{2}(y),tol));
diagbf3 = make_tt_to_operator(round(bfun{3}(y),tol));

Amap = round(Bmap - Lapmap ...
  + amen_mm(diagbf1,Convmap1,tol,'verb',0)...
  + amen_mm(diagbf2,Convmap2,tol,'verb',0)...
  + amen_mm(diagbf3,Convmap3,tol,'verb',0),tol);
fval = amen_mv(Amap,Gbc,tol,'verb',0);
end

function fval = JFbc(y, bfun, Convmap1, Convmap2, Convmap3, Gbc, tol)
fval = ...
  make_tt_to_operator(round(bfun{1}(y).*amen_mv(Convmap1,Gbc,tol,'verb',0),tol))...
  + make_tt_to_operator(round(bfun{2}(y).*amen_mv(Convmap2,Gbc,tol,'verb',0),tol))...
  + make_tt_to_operator(round(bfun{3}(y).*amen_mv(Convmap3,Gbc,tol,'verb',0),tol));
end

function Jmat = Jfuntt(B,Lap, bfun, Conv1,Conv2,Conv3, y,tol)
% Jfun = @(y) B - Lap ...
%   + (diag(bfun{1}(y))*Conv1 + diag((dbfun{1}(y)).*(Conv1*y)))...
%   + (diag(bfun{2}(y))*Conv2 + diag((dbfun{2}(y)).*(Conv2*y)))...
%   + (diag(bfun{3}(y))*Conv3 + diag((dbfun{3}(y)).*(Conv3*y)))...
%   + JFbc(y);
diagbf1 = make_tt_to_operator(round(bfun{1}(y),tol));
diagbf2 = make_tt_to_operator(round(bfun{2}(y),tol));
diagbf3 = make_tt_to_operator(round(bfun{3}(y),tol));

diagdbf1 = make_tt_to_operator(amen_mv(Conv1,y,tol,'verb',0));
diagdbf2 = make_tt_to_operator(amen_mv(Conv2,y,tol,'verb',0));
diagdbf3 = make_tt_to_operator(amen_mv(Conv3,y,tol,'verb',0));

Jmat = round(B - Lap ...
  - (amen_mm(diagbf1,Conv1,tol,'verb',0) + diagdbf1) ...
  - (amen_mm(diagbf2,Conv2,tol,'verb',0) + diagdbf2) ...
  - (amen_mm(diagbf3,Conv3,tol,'verb',0) + diagdbf3), tol);

end