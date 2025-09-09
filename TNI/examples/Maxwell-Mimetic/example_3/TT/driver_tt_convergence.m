addpath(genpath('../../../../matlab/Maxwell-Mimetic/src/'))
addpath(genpath('../../../../matlab/utils/chebfun/'))
addpath(genpath('../../../../matlab/utils/tt-toolbox/'))
addpath(genpath('../../../../matlab/utils/ttfunc/'))

close all; clear; clc;
%% test the convergence of TT

Ns = [3:7];
Err = zeros(numel(Ns),4);
Time = zeros(numel(Ns,1));
MemEz = zeros(numel(Ns),2);
tol = 1e-3;
for i = 1: numel(Ns)
  fprintf('\n------------tt scheme-------------\n')
  nlev = Ns(i);
  [temperr,t, tm, fm] = tt_test3_fn(nlev, tol);
  %store results
  Err(i,:) = temperr;
  Time(i,1) = t;
  MemEz(i,1) = tm;
  MemEz(i,2) = fm;
  save(['../plot_data/tt_example_3.mat']);
  fprintf('results of nlev = %d is saved \n\n', nlev);
end
