function [alpha, ttPsi1, Itertime, convrate, lambda, k_num_iter] = tt_SolveEig_Alpha(ttH, ttS, ttF, ttiV, params,...
  mv_opts, amen_solve_opts)

%{ 
Use k-effective to iteratively solve alph-eigenvalue problem
%}

%% get parameters
n = ttH.n; 
d = numel(n);
niter = params.niter;
epsi = params.epsi;
tt_tol = params.tt_tol;

%% initialization
ttPsi0 = tt_ones(n,d);
evm = 0.01;
convrate = nan(1,niter);
lambda = nan(1,niter);
Itertime = zeros(1,niter);
k_num_iter = zeros(1,niter);
for iter = 0:niter
  ts = datetime;
  fprintf('************************************\n')
  fprintf('Iteration = %d/%d \n convrate = %.5e \n |1-lambda| = %.5e \n elapsed time = %.2f s \n', ...
    iter,niter, convrate(max(iter-1,1)),abs(1-lambda(max(iter-1,1))), sum(Itertime) );
  fprintf('************************************\n')
  %% update alpha
  if ( iter == 0 )
    alpha(iter+1) = 0.0;
  elseif ( iter == 1 )
    alpha(iter+1) = alpha + evm;
  else
    alpha(iter+1) = alpha(iter-1) + ( 1 - k(iter-1) )/(k(iter)-k(iter-1))*(alpha(iter)-alpha(iter-1));
  end
  
  %% solve keffective problem
  LHS = round(ttH + alpha(iter+1)*ttiV,tt_tol);
  [k(iter+1), ttPsi1, ktime] = ...
  tt_fixed_point_eig_solve(LHS, ttS, ttF, params, mv_opts, amen_solve_opts);
  
  k_num_iter(iter+1) = sum(ktime>0);
  %% check convergence
  if iter>0
    lambda(iter) = sum(amen_mv(ttF,ttPsi1,tt_tol, mv_opts))...
      ./sum(amen_mv(ttF,ttPsi0,tt_tol, mv_opts));
    convrate(iter) = norm(ttPsi1-ttPsi0)/norm(ttPsi0);
    conv_test2 = 0;
    num_iter_check = 5;
    if iter> num_iter_check + 2
      conv_test2 = convrate(iter)>convrate(iter-num_iter_check);
    end

    if abs(lambda(iter)-1)<epsi && (convrate(iter) < epsi)
      fprintf('The fixed point scheme converged at iter = %d \n', iter)
      break;
    elseif conv_test2
      fprintf('The convergence rate increases at iter = %d \n', iter)
      break;
    end
  end
  ttPsi0 = ttPsi1;
  Itertime(iter+1) = seconds(datetime-ts);


end