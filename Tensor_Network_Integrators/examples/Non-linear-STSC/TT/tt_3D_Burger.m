addpath(genpath('../../../matlab/Linear-STSC/src/'))
addpath(genpath('../../../matlab/Non-linear-STSC/src/'))
addpath(genpath('../../../matlab/utils/chebfun/'))
addpath(genpath('../../../matlab/utils/tt-toolbox/'))
addpath(genpath('../../../matlab/utils/ttfunc/'))

close all; clear;

testname = "../testcases/Burger3D";

fprintf('Working on testcase %s \n',testname);
load(sprintf('./%s.mat',testname));
%%
a = 0;
b = 6;
tol = 1e-5;
Newtoneps = 1e-5;
tempname = split(testname,'/');
savefilename = sprintf('../plot_data/tt_%s.mat', tempname(end));
%% now the time step is the same as space step
Ns = [4:2:10];
niter=1;
R = cell(numel(Ns),niter);
for j = 1:niter
  for i = [1:numel(Ns)]
    tic
    n = Ns(i);
    SN = n-1;
    N1 = SN+1;
    %%
    fprintf('************************************\n')
    fprintf('n = %d \n',n);
    %%
    [fgsol, NewtonIter, Ctt] = ...
      tt_3D_Burger_Solver(testcase,SN,tol,Newtoneps,a,b);
    tt_time = toc;
    %% check the error
    exacttt = amen_cross(Ctt, @(x) cross_fun_nD(x,testcase.exactfn),tol);
    exacttt = tt_get_inner(exacttt,{2:N1,2:N1-1,2:N1-1,2:N1-1});
    
    Err = norm(fgsol-exacttt)/norm(exacttt);
    fprintf('Solution Error between TT and Exact= %.2e \n',...
      Err);
    fprintf('Newton Iter = %d \n', NewtonIter);
    fprintf('Elapsed time of tt = %.2f seconds \n', tt_time)
    % fprintf('Compression Ratio of ttsol = %.2e \n', compress_ratio_tt(fgsol));
    %% collect results
    c.sol = fgsol;
    c.err = Err;
    c.time = tt_time;
    c.NewtonIter = NewtonIter;
    % c.amensolvedata = amensolvedata;
    R{i,j} = c;

    save(savefilename,'Ns','R','tol','testcase');
  end
end


