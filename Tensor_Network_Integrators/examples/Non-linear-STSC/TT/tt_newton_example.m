addpath(genpath('../../../matlab/Linear-STSC/src/'))
addpath(genpath('../../../matlab/Non-linear-STSC/src/'))
addpath(genpath('../../../matlab/utils/chebfun/'))
addpath(genpath('../../../matlab/utils/tt-toolbox/'))
addpath(genpath('../../../matlab/utils/ttfunc/'))

close all; clear; clc;

%% create test case
n = 12;
d = 4;
r = 2;
epstt = 2e-8;
utrue = tt_rand(n,d,2);
ftemp = @(u,eps) round(amen_cross({u}, @(x) exp(-x),eps,'verb',0) - ...
  amen_cross({u}, @(x) x.^2,eps,'verb',0), eps);
ftemptrue = ftemp(utrue,eps);
F = @(u,eps) round(ftemp(u,eps) - ftemptrue,eps);
Jfun = @(u,eps) round(...
  -make_tt_to_operator(amen_cross({u}, @(x) exp(-x),eps,'verb',0)) ...
  - 2*make_tt_to_operator(u)...
  ,eps);

%%
Itertime = [];
R = [];
ResNorm = [];
LocErr = [];
%%
U0 = 15*tt_rand(n,d,1);
epsk = 1e-1;
eps = 1e-6;
nrmFU0 = norm(F(U0,epsk));
maxiter = 1000;

for k = 1:maxiter
  tsiter = datetime;
  % epsk = epsk*epsfactor;
  fprintf('\n******** Newton iter = %d ************\n\n', k);
  fprintf('epsk = %.5e \n', epsk);

  % compute Jacobian with JFbc
  Jmatk = Jfun(U0,epstt);

  % Jmatk = Jfuntt(Btt,Laptt,Ku,dKudu,dFudu,U0,tol);

  % use GMRES to solve for du
  b = -F(U0,epstt);

  du = amen_solve2(Jmatk,b,epsk,'nswp',30,'verb',1);

  tegmres = datetime;
  fprintf('Linear Solve Time = %.2f\n', seconds(tegmres-tsiter));

  %update U1
  for alpha = [1,1/2,1/4,1/8,1/16]
    fprintf('testing alpha = %.5f \n',alpha);
    U1 = round(U0 + alpha*du, max(epsk,epstt));
    if norm(F(U1,epstt))<norm(F(U0,epstt))
      break;
    end
  end
  
  % save the rank
  R = [R, U1.r];

  %compute local error
  local_err=norm(U1-U0)/norm(U0);
  res_norm = norm(F(U1,epstt))/nrmFU0;
  ResNorm = [ResNorm, res_norm];
  LocErr = [LocErr, local_err];
  % #* set epsk
  epsk = min(min([local_err,res_norm]),epsk);
  % epsk = 0.5*min([local_err,res_norm]);

  fprintf('apha = %.2e \n', alpha);
  fprintf('u_err = %.5e,  Fu_ratio = %.5e \n', local_err, res_norm)
  Itertime = [Itertime, seconds(datetime-tsiter)];
  fprintf('iter time = %.2f \n', seconds(datetime-tsiter));
  if (local_err<eps) && (res_norm<eps)
    break;
  end
  U0 = U1;
end


%%
uest = U0;
error = norm(utrue - uest)/norm(utrue);
fprintf('Err = %.2e \n', error)

save('../plot_data/tt_newton_example.mat','Itertime','R','error','ResNorm','LocErr');



