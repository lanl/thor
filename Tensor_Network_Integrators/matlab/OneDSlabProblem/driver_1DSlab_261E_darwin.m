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
param.nx = 2048; % index i - number of EDGES
param.niter = 200;
param.fixed_point_tol = 1e-6;
tol = 1e-6; %both convergence and tt

savefilename = './Results/OneD_261E_running_tt.mat';
for nx = 2048 %[256,512,1024,2048]
  param.nx = nx;
  %% TT fixed-point
  %tic
  [ttPsi1, ttk, ttconvrate, Itertime] = tt_1DSlab_EVP_fn(param, tol, savefilename);
  %tt_time = toc

  %%
  filename = sprintf('./Results/OneD_261E_tt_nx%d_result_1e6.mat',param.nx);
  save(filename,'ttPsi1','ttk','ttconvrate',...
    'Itertime','param','tol');
  fprintf('Elapsed time = %.5f s \n', sum(Itertime));
  fprintf('Result is saved !!!');
end