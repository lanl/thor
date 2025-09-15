%%
% this script will run the exponential test of amencross mp 
% with all Genz functions
close all; clear;clc;
run setup.m;
% run setup_advanpix_drac.m;

%% parameters
Ds = [10 20 50 100];
ncell = 1;
Nq = [1:10];
method = 'GL';
tol = 1e-16;
R = cell(2,numel(Ds), numel(Nq)); % cell array - dimensions x functions x orders
for itest = 8 %choose the test problem
  for i = 1:numel(Ds)
    d = Ds(i);
    for j = 1:numel(Nq)
      norder = Nq(j);
      %       Errors(i,itest,j) = tt_GL_fn_poly(d,30,norder,ncell,1e-15);
      fprintf('\n\n Ds = %d, itest = %d, nq = %d \n',d,itest,norder);
      
      r = tt_GL_fn(itest, d, 1, method, norder, tol);
      R{1,i,j} = r;
      r = tt_GL_fn(itest, d, 2, method, norder, tol);
      R{2,i,j} = r;
    end
  end
end

%%
if 0
  filename = sprintf('T10pconv_%s.mat',method);
  save('T10_pconv_result.mat','Ds','Nq','R');
end
%% plot for T10
load('T10_pconv_result.mat');
if 1
  load('T10_pconv_result.mat');
  figure()
  cm ={'r-s','b--^','k:o'};
  for i = 1:numel(Ds)
    subplot(2,2,i)
    hold on;
    err1 = cellfun(@(c) c.Err, R(1,i,:));
    err2 = cellfun(@(c) c.Err, R(2,i,:));
    plot(Nq,err1(:),'-^');
    plot(Nq,err2(:),'-*');
    xlabel('Number of Nodes');
    ylabel('Approximation Error');
    title(sprintf('d=%d',Ds(i)));
    set(gca,'YScale','log');
  end
  legend({'ncell=1','ncell=2'})
  
  sgtitle(sprintf(''));
end