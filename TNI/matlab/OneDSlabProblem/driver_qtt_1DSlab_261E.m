close all; clear; clc;
if ~isdeployed
  run setup1DSlab.m;
end

%% setup parameters
param.a = 0;
param.b = 4.2;
% scatterting parameters
param.sigma_s = importdata('./OneDSlab_261E/sigma_s.txt');
% fission parameters
param.nusigma_f = load('./OneDSlab_261E/nusigma_f.txt');
param.chi = load('./OneDSlab_261E/chi.txt');

param.sigma = load('./OneDSlab_261E/sigma_t.txt');
param.nmu = 256; % index l
param.nx = 256; % index i - number of EDGES
param.niter = 200;
param.fixed_point_tol = 1e-4; %convergence
tol = 1e-4; %tt

savefilename = './Results/OneD_261E_running_qtt.mat';
for nx = [256, 512, 1024, 2048]
  param.nx = nx;
  filename = sprintf('./Results/OneD_261E_qtt_nx%d_result_1e4.mat',param.nx);

  %% TT fixed-point
  %tic
  [ttPsi1, ttk, ttconvrate, Itertime] = qtt_1DSlab_EVP_fn(param, tol, savefilename);
  %tt_time = toc

  %%
  save(filename,'ttPsi1','ttk','ttconvrate',...
    'Itertime','param','tol');
  fprintf('Elapsed time = %.5f s \n', sum(Itertime));
  fprintf('Result is saved !!!');
end