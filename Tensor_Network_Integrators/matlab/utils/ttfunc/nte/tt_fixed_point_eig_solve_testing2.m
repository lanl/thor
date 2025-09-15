function [k, ttPsi1, Itertime, convrate, lambda, fluxnegnum, ttPsi1rank] = ...
  tt_fixed_point_eig_solve_testing2(ttH, ttS, ttF, params, mv_opts, amen_solve_opts)

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
ttPsi1rank = [];
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
  fprintf('Iteration = %d/%d --convrate = %.5e |1-lambda| = %.5e --k = %.5e\n', ...
    iter,niter, convrate(max(iter-1,1)),abs(1-lambda(max(iter-1,1))),k);
    fprintf('epsk = %.5e\n',epsk)

  RHS = amen_mv(ttS,ttPsi0,tt_tol, mv_opts) ...
    + 1/k*amen_mv(ttF,ttPsi0,tt_tol, mv_opts);

  %  Solve a linear system for Psi1 = Htt\RHS;
  amen_solve_opts{end} = round(ttPsi0, tt_tol);
  ttPsi1 = amen_solve2(ttH, RHS, tt_tol, amen_solve_opts);
  
  ttPsi1rank = [ttPsi1rank; ttPsi1.r];

  %% compute angular flux
  % G = core2cell(ttPsi1);
  % G{4} = sum(G{4},2);
  % G{3} = tensorprod(G{3},G{4},3,1);
  % display(ttH)
  % display(ttPsi1)
  % display(G)
  % angttPsi1 = cell2core(tt_tensor,G(1:3));
  angttPsi1 = compute_angflux_fn(ttPsi1);
  fluxnegnum(iter) = sum(full(angttPsi1)<0); % keep track of negative number in flux
  

  %% update eigenvalue
  lambda(iter) = sum(amen_mv(ttF,ttPsi1,tt_tol, mv_opts))...
    ./sum(amen_mv(ttF,ttPsi0,tt_tol, mv_opts));

  k = k*lambda(iter);

  %%
  convrate(iter) = norm(ttPsi1-ttPsi0)/norm(ttPsi0);

  %% update epsk
  % epsk = min(max([abs(lambda(iter)-1),convrate(iter),tt_tol]),epsk);
  % epsk = min([abs(lambda(iter)-1),convrate(iter),tt_tol]);

  epsk = exp(0.5*(log(abs(lambda(iter)-1)) + log(convrate(iter))))
  epsk = min(epsk,1e-2);

  Itertime(iter) = toc;

  if abs(lambda(iter)-1)<epsi && (convrate(iter) < epsi)
    fprintf('The fixed point scheme converged at iter = %d \n', iter)
    fprintf('Iteration = %d/%d \nconvrate = %.5e \n|1-lambda| = %.5e \n k = %.5e \nelapsed time = %.2f s \n', ...
      iter,niter, convrate(max(iter,1)),abs(1-lambda(max(iter,1))),k, sum(Itertime) );
    break;
  end
  ttPsi0 = ttPsi1;
end
end