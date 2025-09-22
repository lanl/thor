function x=predict(x,u,A,F,B,Tfinal,dt,tint_order,isA,isF,eps_tt)
  %
  x = round(x,eps_tt);
  %
  tol = 1e-8;
  %
  rhs = @(tt,xx,ee)calculateRHS(tt,xx,eps_tt,u,A,F,B,isA,isF,x.n);
  %
  x = RK45(rhs, [0 Tfinal], x, dt, tol, eps_tt);
  %
%  x = x0;
%  %
%  tvals = 0:dt:Tfinal;
%  %
%  dtmin  = 1e-4;
%  dtmax  = 0.1;
%  dtinit = dt;
%  %
%  rtol = 1e-4;
%  atol = 1e-14*tt_ones([Nx Ny Nz Neq],4); 
%  %
%  Btable = butcher("Dormand-Prince-ERK");
%  numstg = numel(Btable(1,:))-1;
%  %
%  solverparameters = {"caller","rtol","atol","dtmin","dtmax","dtinit"};
%  par              = vars2struct(solverparameters{:});
%  par.ettinit      = eps_tt;
%  %
%  tic;
%  [~,x,~,~] = tt_solve_ERK(rhs,tvals,x,Btable,par);
%  toc
%  %
%  x = x{end};
  %
end