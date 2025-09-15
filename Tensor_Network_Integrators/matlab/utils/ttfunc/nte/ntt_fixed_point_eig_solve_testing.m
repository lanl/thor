function [k, ttPsi1, Itertime, convrate, lambda, angnegnum] = ...
  ntt_fixed_point_eig_solve_testing(ttH, ttS, ttF, ...
  params, mv_opts, amen_solve_opts,round_type)

%% get parameters
n = ttH.n; 
d = numel(n);
% fixed point schemes
% ttPsi0 = tt_rand_pos(n,d,1);
niter = params.niter;
epsi = params.epsi;
tt_tol = params.tt_tol;

% mv_nswp = params.mv_nswp;
% amens_solve_nswp = params.amens_solve_nswp;
% kickrank = params.kickrank;
amen_cross = @amen_cross_zero;
%% initialization
ttPsi0 = tt_ones(n,d);
k = 1.0;

convrate = nan(1,niter);
lambda = nan(1,niter);
Itertime = zeros(1,niter);
% mv_nswp = 50;
% amens_solve_nswp = 50;
% kickrank = 2;
try
  epsk = params.epsk;
catch
  epsk = 1e-2;
end

% get Hinv
try 
  ttHinv = params.ttHinv;
  useHinv = 1;
catch
  useHinv = 0;
end


%% Fixed-point opertoin
for iter = 1:niter
  tic
  fprintf('Iteration = %d/%d \nconvrate = %.5e \n|1-lambda| = %.5e \n k = %.5e \nelapsed time = %.2f s \n', ...
    iter,niter, convrate(max(iter-1,1)),abs(1-lambda(max(iter-1,1))),k, sum(Itertime) );
  RHS = amen_mv(ttS,ttPsi0,tt_tol, mv_opts) ...
    + 1/k*amen_mv(ttF,ttPsi0,tt_tol, mv_opts);

  %  Solve a linear system for Psi1 = Htt\RHS;
  if useHinv
    ttPsi1 = amen_mv(ttHinv,RHS,epsk);
  else
    amen_solve_opts{end} = round(ttPsi0, tt_tol);
    ttPsi1 = amen_solve2(ttH, RHS, epsk, amen_solve_opts);
  end
  
  %ttPsi1 = dmrg_solve3(ttH,RHS,tt_tol);
  % angttPsi1 = compute_angflux_fn(ttPsi1);
  % angnegnum(iter) = sum(full(angttPsi1)<0);
  
  fprintf('Number of negative entry in Psi = %d \n',sum(full(ttPsi1)<0));
 
  % if sum(full(angttPsi1)<0)
  %   error('There is negative entry in angular flux');
  % end

  %% round ttPsi1
  if strcmp(round_type,'cross')
    ttPsi1 = amen_cross({ttPsi1},@(x) max(x,0),tt_tol*0.1,'verb',0);
  elseif strcmp(round_type,'full')
    %full projection
    Psi1 = max(tt_full(ttPsi1),0);
    % fprintf('Number of negative entry in Psi = %d \n', sum(Psi1(:)<0));
    ttPsi1 = ntt_decomp_fn(Psi1,int8(ttPsi1.r(2:end-1)*1.1));
  end
 
  % angttPsi1 = compute_angflux_fn(ttPsi1);

  angnegnum(iter) = sum(full(ttPsi1)<0);
  fprintf('After amen cross projection \n')
  fprintf('Number of negative entry in ttPsi = %d \n', sum(full(ttPsi1)<0));
  %% update eigenvalue
  lambda(iter) = sum(amen_mv(ttF,ttPsi1,tt_tol, mv_opts))...
    ./sum(amen_mv(ttF,ttPsi0,tt_tol, mv_opts));
  k = k*lambda(iter);
  %   fprintf('iter = %d - ||Psi1 - Psi0|| = %.5e \n',iter, norm(ttPsi1-ttPsi0));
  %relative error
  
  %%
  convrate(iter) = norm(ttPsi1-ttPsi0)/norm(ttPsi0);

  %% update epsk
  epsk = min(max([abs(lambda(iter)-1),convrate(iter),tt_tol]),epsk);

  % tt_tol = epsk;
  
  conv_test2 = 0;
  num_iter_check = 100;
  % if iter> num_iter_check
  %   conv_test2 = convrate(iter)>convrate(iter-num_iter_check);
  % end
  
  Itertime(iter) = toc;

  if abs(lambda(iter)-1)<epsi || (convrate(iter) < epsi)
    fprintf('The fixed point scheme converged at iter = %d \n', iter)
    fprintf('Iteration = %d/%d \nconvrate = %.5e \n|1-lambda| = %.5e \n k = %.5e \nelapsed time = %.2f s \n', ...
      iter,niter, convrate(max(iter,1)),abs(1-lambda(max(iter,1))),k, sum(Itertime) );
    break;
  elseif conv_test2
    fprintf('The convergence rate increases at iter = %d \n', iter)
    fprintf('Iteration = %d/%d \nconvrate = %.5e \n|1-lambda| = %.5e \n k = %.5e \nelapsed time = %.2f s \n', ...
    iter,niter, convrate(max(iter,1)),abs(1-lambda(max(iter,1))),k, sum(Itertime) );
    break;
  end
  ttPsi0 = ttPsi1;
  

end
%normalize ttPsi1
% ttPsi1 = ttPsi1/norm(ttPsi1);
end