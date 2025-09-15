%%
% This script will run the accuracy test
% methods: Trapezoidal and Simpson
close all; clear;clc;
run setup.m;

%% parameters
Ds = [2];
Ns = [10,20];
method = 'GL'; %method
norder = 2; % order
tol = 1e-16;
R = cell(5,numel(Ds), numel(Ns)); % cell array - dimensions x functions x orders

for itest = [2,4,5] %choose the test problem
  for i = 1:numel(Ds) %loop through number of dimensions
    d = Ds(i);
    for j = 1:numel(Ns) %loop through number of intervals
      ncell = Ns(j);
      fprintf('\n\n %s - Ds = %d, itest = %d, ncell = %d\n',...
        method, d, itest, ncell);
      r = tt_GL_fn(itest, d, ncell, method, norder, tol);
      R{itest,i,j} = r;
    end
  end
end
%%
if 0
  fprintf('Saving result ...')
  filename = sprintf('hconv_test_result_Genz245_%s.m',method);
  save(filename,'Ds','Ns','method','tol','R');
  fprintf('Done \n');
end
%% check the convergence for accuracy
if 1
  itest = 2
  err = squeeze(cellfun(@(c) c.Err, R(itest,1,:)));
  figure()
  plot(log(Ns),log(err));
  temp = polyfit(log(Ns),log(err),1);
  sgtitle(sprintf('convergence rate = %.5f \n',log(err(2)/err(1))/log(1/2)));
end