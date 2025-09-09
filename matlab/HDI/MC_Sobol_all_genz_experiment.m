close all; clear; clc;
% this function will test quasi-MC methods on all genz function at different
% tolerance
% it seems cubSobol_g gets the better error, but cubLattic_g is faster.
% So we will test both
run setup.m;
run setup_GAIL.m;
%%
Ds = [10, 20, 50, 100, 200, 500, 1000];
Tols = 10.^[-1 -3 -5 -7 -9 -11 -13];
R = cell(5,numel(Ds), numel(Tols)); % functions x Tols

for itest = 1:5
  for i = 1: numel(Ds)
    d = Ds(i);
    %set the problem and the domain
    [f,Itrue] = select_integrand_function(itest,d);
    hyperbox = [zeros(1,d); ones(1,d)];
    fprintf('MCcubSobol - Dimension = %d \n',d);
    for j = 1: numel(Tols)
      tol = Tols(j);
      in_param.abstol = tol;
      in_param.reltol= tol;
      
      %computing
      tic
      Q = cubSobol_g(f,hyperbox,in_param);
      R{itest,i,j}.time = toc;
      R{itest,i,j}.Err = norm(Itrue-Q);
      fprintf('tol = %.2e, err = %.2e, time = %.2f\n',tol, ...
        norm(Itrue-Q), R{itest,i,j}.time);
      
    end
  end
  save('MC_cubSobol_result.mat','Ds','Tols','R');
  fprintf('Result file is saved \n');
end