function [k, ttPsi1, Itertime, convrate, lambda] = ...
  tt_fixed_point_eig_solve(ttH, ttS, ttF, params, mv_opts, amen_solve_opts)

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
%% initialization
ttPsi0 = tt_ones(n,d);
k = 1.0;

convrate = nan(1,niter);
lambda = nan(1,niter);
Itertime = zeros(1,niter);
% mv_nswp = 50;
% amens_solve_nswp = 50;
% kickrank = 2;

for iter = 1:niter
  tic
  fprintf('Iteration = %d/%d \nconvrate = %.5e \n|1-lambda| = %.5e \n k = %.5e \nelapsed time = %.2f s \n', ...
    iter,niter, convrate(max(iter-1,1)),abs(1-lambda(max(iter-1,1))),k, sum(Itertime) );
  RHS = amen_mv(ttS,ttPsi0,tt_tol, mv_opts) ...
    + 1/k*amen_mv(ttF,ttPsi0,tt_tol, mv_opts);

  %  Solve a linear system for Psi1 = Htt\RHS;
  % opts = {'local_prec','cjacobi','kicktype','svd','resid_damp',10,'nswp',10};
  % amen_solve_opts = {'nswp', amens_solve_nswp,'kickrank', kickrank};
  amen_solve_opts{end} = round(ttPsi0, tt_tol);
  ttPsi1 = amen_solve2(ttH, RHS, epsi, amen_solve_opts);
  
  %ttPsi1 = dmrg_solve3(ttH,RHS,tt_tol);
  %% update eigenvalue
  lambda(iter) = sum(amen_mv(ttF,ttPsi1,tt_tol, mv_opts))...
    ./sum(amen_mv(ttF,ttPsi0,tt_tol, mv_opts));
  k = k*lambda(iter);
  %   fprintf('iter = %d - ||Psi1 - Psi0|| = %.5e \n',iter, norm(ttPsi1-ttPsi0));
  %relative error
  
  %%
  convrate(iter) = norm(ttPsi1-ttPsi0)/norm(ttPsi0);

  conv_test2 = 0;
  num_iter_check = 5;
  if iter> num_iter_check
    % conv_test2 = convrate(iter)>convrate(iter-num_iter_check);
  end
  
  if abs(lambda(iter)-1)<epsi && (convrate(iter) < epsi)
    fprintf('The fixed point scheme converged at iter = %d \n', iter)
    break;
  elseif conv_test2
    fprintf('The convergence rate increases at iter = %d \n', iter)
    break;
  end
  ttPsi0 = ttPsi1;
  Itertime(iter) = toc;

end
%normalize ttPsi1
% ttPsi1 = ttPsi1/norm(ttPsi1);
end