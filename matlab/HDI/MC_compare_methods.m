close all; clear; clc;
%this script will compare different MC methods

%% get function
d = 100;
tol = 1e-12;
itest = 1;
%%
[f,Itrue] = select_integrand_function(itest,d);
hyperbox = [zeros(1,d); ones(1,d)];

in_param.abstol = tol;
in_param.reltol= tol;
%% cubMC_g -- THIS IS ONLY MC method
if 1
  tic
  Q = cubMC_g(f,hyperbox,in_param);
  fprintf('cubMC err = %.2e \n',abs(Itrue-Q));
  toc;
end
%% cubSobol_g -- THIS IS quasi-MC with Sobol sequence
if 1
  tic
  Q = cubSobol_g(f,hyperbox,in_param);
  fprintf('cubSobol_g err = %.2e \n',abs(Itrue-Q));
  toc;
end

%% cubLattice_g
if 1
%   in_param.transform='';
  tic
  Q = cubLattice_g(f,hyperbox,in_param);
  fprintf('cubLattice_g err = %.2e \n',abs(Itrue-Q));
  toc;
end


