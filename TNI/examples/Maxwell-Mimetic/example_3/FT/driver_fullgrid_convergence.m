addpath(genpath('../../../../matlab/Maxwell-Mimetic/src/'))
addpath(genpath('../../../../matlab/utils/chebfun/'))
addpath(genpath('../../../../matlab/utils/tt-toolbox/'))
addpath(genpath('../../../../matlab/utils/ttfunc/'))

close all; clear; clc;
%%
% test the convergence of fullgrid

Ns = 3:5;
Err = zeros(numel(Ns),4);
Time = zeros(numel(Ns,1));
fprintf('\n-------------------------\n')
for i = 1: numel(Ns)
  nlev = Ns(i);
  [temperr,t] = fullgrid_test3_fn(nlev);
  Err(i,:) = temperr;
  Time(i,1) = t;
  save('../plot_data/fullgrid_example_3.mat');
  fprintf('results of nlev = %d is saved \n\n', nlev);
end
