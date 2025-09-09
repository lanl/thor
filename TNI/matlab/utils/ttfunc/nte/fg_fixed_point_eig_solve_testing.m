function [k, Psi1, Itertime, convrate, lambda] = ...
  fg_fixed_point_eig_solve_testing(H, S, F, params)

%% get parameters
% n = H.n; 
% d = numel(n);
% fixed point schemes
% ttPsi0 = tt_rand_pos(n,d,1);
niter = params.niter;
epsi = params.epsi;

%% initialization
% Psi0 = tt_full(tt_ones(n,d));
Psi0 = ones(size(H,1),1);
k = 1.0;

convrate = nan(1,niter);
lambda = nan(1,niter);
Itertime = zeros(1,niter);

%% Fixed-point opertoin
for iter = 1:niter
  tic
  % fprintf('Iteration = %d/%d \nconvrate = %.5e \n|1-lambda| = %.5e \n k = %.5e \nelapsed time = %.2f s \n', ...
  %   iter,niter, convrate(max(iter-1,1)),abs(1-lambda(max(iter-1,1))),k, sum(Itertime) );

  fprintf('Iteration = %d/%d --convrate = %.5e |1-lambda| = %.5e --k = %.5e\n', ...
    iter,niter, convrate(max(iter-1,1)),abs(1-lambda(max(iter-1,1))),k);

  % RHS = amen_mv(S,ttPsi0,tt_tol, mv_opts) ...
  %   + 1/k*amen_mv(F,ttPsi0,tt_tol, mv_opts);

  RHS = S*Psi0 + 1/k*F*Psi0;
  %  Solve a linear system for Psi1 = Htt\RHS;
  Psi1 = H\RHS;
  % if useHinv
  %   ttPsi1 = amen_mv(ttHinv,RHS,epsk);
  % else
  %   amen_solve_opts{end} = round(ttPsi0, tt_tol);
  %   ttPsi1 = amen_solve2(H, RHS, epsk, amen_solve_opts);
  % end

  % Psi1 = max(Psi1,0);  
  fprintf('Number of negative entry in Psi = %d \n',sum(Psi1(:)<0));
  

  %% update eigenvalue
  % lambda(iter) = sum(amen_mv(F,ttPsi1,tt_tol, mv_opts))...
  %   ./sum(amen_mv(F,ttPsi0,tt_tol, mv_opts));
  % 
  lambda(iter) = sum(F*Psi1,'all')./sum(F*Psi0,'all');
  k = k*lambda(iter);
  
  %%
  convrate(iter) = norm(Psi1-Psi0)/norm(Psi0);
  % convrate(iter) = max(abs((Psi1-Psi0)./Psi0));
  Itertime(iter) = toc;

  if abs(lambda(iter)-1)<epsi && (convrate(iter) < epsi)
    fprintf('The fixed point scheme converged at iter = %d \n', iter)
    fprintf('Iteration = %d/%d \nconvrate = %.5e \n|1-lambda| = %.5e \n k = %.10f \nelapsed time = %.2f s \n', ...
      iter,niter, convrate(max(iter,1)),abs(1-lambda(max(iter,1))),k, sum(Itertime) );
    break;
  end
  Psi0 = Psi1;
  

end
%normalize ttPsi1
% ttPsi1 = ttPsi1/norm(ttPsi1);
end