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
utrue = full(tt_rand(n,d,2));
ftemp = @(u) exp(-u) - u.^3;
F = @(u) ftemp(u) - ftemp(utrue);
Jfun = @(u) -diag(exp(-u)) - 3*diag(u.^2);

%%
Itertime = [];
ResNorm = [];
LocErr = [];
%%
U0 = full(10*tt_rand(n,d,1));
eps = 1e-6;

nrmFU0 = norm(F(U0));
maxiter = 1000;

for k = 1:maxiter
  tsiter = datetime;
  fprintf('\n******** Newton iter = %d ************\n\n', k)

  %compute Jacobian with new u
  Jmatk = Jfun(U0);

  b = -F(U0);
  du = Jmatk\b;

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
  ResNorm = [ResNorm, res_norm];
  LocErr = [LocErr, local_err];

  fprintf('apha = %.2e \n', alpha);
  fprintf('u_err = %.5e,  Fu_ratio = %.5e \n', local_err, res_norm)
  Itertime = [Itertime,seconds(datetime-tsiter)];
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

save('../plot_data/fg_newton_example.mat','Itertime','error','ResNorm','LocErr');
