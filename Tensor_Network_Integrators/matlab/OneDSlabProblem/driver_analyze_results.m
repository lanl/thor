close all; clear; clc;
run setup1DSlab.m;

%%
R{1} = load('./Results/OneD_261E_result_1e6.mat');
R{2} = load('./Results/OneD_261E_result_1e10.mat');
%%

%show the convrate
figure(1)
plot(R{1}.ttconvrate,'sr','DisplayName','1e-6')
hold on;
plot(R{2}.ttconvrate,'b*', 'DisplayName','1e-10')
set(gca,'YScale','log')
xlabel('Fixed point iteration')
ylabel('Convergence rate')
title('convergence rate')
legend()

%%
figure(2)
plot(R{1}.Itertime(R{1}.Itertime>0),'sr','DisplayName',...
  sprintf('1e-6 - Total time = %.2f minutes',sum(R{1}.Itertime(R{1}.Itertime>0))/60));
hold on;
plot(R{2}.Itertime(R{2}.Itertime>0),'b*', 'DisplayName', ...
  sprintf('1e-10 - Total time = %.2f minutes',sum(R{2}.Itertime(R{2}.Itertime>0))/60))
set(gca,'YScale','log')
xlabel('Fixed point iteration')
ylabel('Elapsed time (s)')
legend()

%% 
ttPsidiff = norm(R{1}.ttPsi1 - R{2}.ttPsi1)
kdiff = abs(R{1}.ttk - R{2}.ttk)
for i = 1:2
  compress_ratio_tt(R{i}.ttPsi1)
end