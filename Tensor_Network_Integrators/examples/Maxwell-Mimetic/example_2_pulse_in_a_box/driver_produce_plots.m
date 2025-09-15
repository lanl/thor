% produce figures
close all; clear; clc;

%% divergence figure
%load results
tdata = load('plot_data/tt_results.mat');
fdata = load('plot_data/fg_results.mat');
t = linspace(fdata.T0,fdata.T1,numel(fdata.divH));
figure()
plot(t,fdata.divH,'.r');
hold on;
plot(t,tdata.divH,'.--b');
% ylim([0,1]);
set(gca,'Yscale','log');
ylabel('Divergence of H');
xlabel('time (s)');
legend({'FG','TT'});
title('divergence of H');

%% energy figure
t = linspace(fdata.T0,fdata.T1,numel(fdata.Energy_E));
figure()
subplot(1,3,1)
plot(t,fdata.Energy_E,'-r');
hold on;
plot(t,fdata.Energy_H,'-b');
% set(gca,'Yscale','log');
ylabel('Energy');
xlabel('time (s)');
legend({'E','H'});
title('Full Grid');

subplot(1,3,2)
plot(t,tdata.Energy_E,'-r');
hold on;
plot(t,tdata.Energy_H,'-b');
% set(gca,'Yscale','log');
ylabel('Energy');
xlabel('time (s)');
legend({'E','H'});
title('Tensor Train');

subplot(1,3,3)
hold on;
plot(t,abs(fdata.Energy_E - tdata.Energy_E),'-r');
plot(t,abs(fdata.Energy_H - tdata.Energy_H),'-b');
xlabel('time (s)');
ylabel('Error')
set(gca,'YScale','log')
legend({'E','H'})
title('Difference in Energy')
sgtitle('Energy Conservation','FontSize',20)



