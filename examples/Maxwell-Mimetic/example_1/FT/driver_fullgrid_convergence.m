addpath(genpath('../../../../matlab/Maxwell-Mimetic/src/'))
addpath(genpath('../../../../matlab/utils/chebfun/'))
addpath(genpath('../../../../matlab/utils/tt-toolbox/'))
addpath(genpath('../../../../matlab/utils/ttfunc/'))

close all; clear; clc;

% test the convergence of fullgrid

Ns = 1:3;
Err = zeros(numel(Ns),4);
Time = zeros(numel(Ns,1));
final_sols = cell(numel(Ns,1));
fprintf('\n-------------------------\n')
for i = 1: numel(Ns)
  nlev = Ns(i);
  [temperr,t, Ex, Ey, Ez, Hx, Hy, Hz] = fullgrid_test1_k1_fn(nlev);
  %save variables
  Err(i,:) = temperr;
  Time(i,1) = t;
  final_sols{i,1}.Ex = Ex;
  final_sols{i,1}.Ey = Ey;
  final_sols{i,1}.Ez = Ez;
  final_sols{i,1}.Hx = Hx;
  final_sols{i,1}.Hy = Hy;
  final_sols{i,1}.Hz = Hz;

  save('../plot_data/fullgrid_test1_k1.mat');
  fprintf('results of nlev = %d is saved \n\n', nlev);
end
