addpath(genpath('../../../matlab/Linear-STSC/src/'))
addpath(genpath('../../../matlab/utils/chebfun/'))
addpath(genpath('../../../matlab/utils/tt-toolbox/'))
addpath(genpath('../../../matlab/utils/ttfunc/'))
close all; clear;

%% load test cases
testname = 'ConstantCoeff'; %example 1 - constant coefficient
% testname = 'VariableCoeff'; %example 2 - variable coefficient
% testname = 'NonSmooth'; %example 3 - variable coefficient

fprintf('Working on testcase %s \n',testname);
load(sprintf('../testcases/%s.mat',testname));

%%
a = -1;
b = 1;
tol = 1e-13;
savefilename = sprintf('../plot_data/fullgrid_ST_FD_%s.mat', testname);

%% now the time step is the same as space step
Ns = [4,6,8];
R = cell(numel(Ns),1);

for i = 1:numel(Ns)
  tic
  n = Ns(i);
  TN = n;
  SN = n-1;
  N1 = SN+1;

  %%
  [fgsol, Ctt] = fg_ST_FD_FD_3D_Solver(testcase,TN,SN,tol,a,b);
  tt_time = toc;
  %% check the error
  exacttt = amen_cross(Ctt, @(x) cross_fun_nD(x,testcase.exactfn),tol);
  exacttt = tt_get_inner(exacttt,{2:TN,2:N1-1,2:N1-1,2:N1-1});
  exacttt = full(exacttt);
  Err = norm(fgsol-exacttt)/norm(exacttt);
  fprintf('Solution Error between TT and matrix= %.2e \n',...
    Err);

  fprintf('Elapsed time of qtt = %.2f seconds \n', tt_time)
  % fprintf('Compression Ratio of ttsol = %.2e \n', compress_ratio_tt(fgsol));


  %% collect results
  c.sol = fgsol;
  c.err = Err;
  c.time = tt_time;
  % c.amensolvedata = amensolvedata;
  R{i} = c;

  save(savefilename,'Ns','R','tol','testcase');
end




