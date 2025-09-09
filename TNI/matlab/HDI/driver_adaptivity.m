%% 
% this script will run the exponential test of amencross mp 
% with all Genz functions
close all; clear;clc;
run setup.m;
% run setup_advanpix_drac.m;

%% parameters
Ds = [10 20 50 100];
Ns = [2,3,4];
norder = 1;
method = 'GL';
tol = 1e-16;
R = cell(1,numel(Ds), numel(Ns)); % cell array - dimensions x functions x orders
for itest = 7 %choose the test problem
  for i = 1:numel(Ds)
    d = Ds(i);
    for j = 1:numel(Ns)
      ncell = Ns(j);
      %       Errors(i,itest,j) = tt_GL_fn_poly(d,30,norder,ncell,1e-15);
      fprintf('\n\n Ds = %d, itest = %d, ncell = %d \n',d,itest,ncell);
      r = tt_GL_fn(itest, d, ncell, method, norder, tol);
      R{1,i,j} = r;
    end
  end
end

%%
if 1
  save('./Results/adapt_34_results.mat','Ds','norder','Ns','R');
end
%%
if 1
  figure()
%   cm ={'r-s','b--^','k:o'};
  for i = 1:numel(Ns)
    %   subplot(1,3,i)
    hold on;
    err = cellfun(@(c) c.Err, R(1,:,i));
    plot(Ds,err(:),'-*');
    xlabel('d');
    ylabel('Err');
  end
  legend(string(num2cell(Ns)))
  set(gca,'YScale','log');
  sgtitle(sprintf('tt Convergence - Genz 3 func'));
end