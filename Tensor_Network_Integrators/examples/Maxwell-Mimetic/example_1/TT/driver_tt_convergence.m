addpath(genpath('../../../../matlab/Maxwell-Mimetic/src/'))
addpath(genpath('../../../../matlab/utils/chebfun/'))
addpath(genpath('../../../../matlab/utils/tt-toolbox/'))
addpath(genpath('../../../../matlab/utils/ttfunc/'))


close all; clear; clc;
%% test the convergence of fullgrid

Tols = [1e-5];
for jj = 1: numel(Tols)
  Ns = [1:5];
  Err = zeros(numel(Ns),4);
  Time = zeros(numel(Ns,1));
  MemEz = zeros(numel(Ns),2);
  final_sols = cell(numel(Ns),1);
  tol = Tols(jj);
  for i = 1: numel(Ns)
    fprintf('\n------------tt scheme-------------\n')
    fprintf('tol = %.2e\n', tol)
    nlev = Ns(i);
    [temperr,t, tm, fm, tEx, tEy, tEz, tHx, tHy, tHz] = tt_test1_k1_fn(nlev, tol);

    %store results
    Err(i,:) = temperr;
    Time(i,1) = t;
    MemEz(i,1) = tm;
    MemEz(i,2) = fm;

    final_sols{i,1}.tEx = tEx;
    final_sols{i,1}.tEy = tEy;
    final_sols{i,1}.tEz = tEz;
    final_sols{i,1}.tHx = tHx;
    final_sols{i,1}.tHy = tHy;
    final_sols{i,1}.tHz = tHz;

    filename = sprintf('../plot_data/tt_test1_k1_tol_%s.mat',num2str(tol));
    save(filename,'Ns','Err','Time','MemEz','final_sols','tol');
    fprintf('results of nlev = %d is saved \n\n', nlev);
  end
end