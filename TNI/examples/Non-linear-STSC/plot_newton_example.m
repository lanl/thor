addpath(genpath('../../matlab/Linear-STSC/src/'))
addpath(genpath('../../matlab/Non-linear-STSC/src/'))
addpath(genpath('../../matlab/utils/chebfun/'))
addpath(genpath('../../matlab/utils/tt-toolbox/'))
addpath(genpath('../../matlab/utils/ttfunc/'))

close all; clear; clc;

%% load data
fgdata = load('plot_data/fg_newton_example.mat');
ttstdata = load('plot_data/tt_newton_example.mat');
n = 12;
d = 4;
%% Compute compression ratio
ttstcomp = [];
for i = 1:numel(ttstdata.Itertime)
  r = ttstdata.R(:,i);
  temp = tt_rand(12,4,r);
  ttstcomp = [ttstcomp, compress_ratio_tt(temp)];
end

%%
fgmrk = '-*k';
ttstmrk = '--ok';

figure()

subplot(1,3,1)
grid on;
hold on;
plot(fgdata.Itertime,fgmrk,'DisplayName','Full Grid')
plot(ttstdata.Itertime,ttstmrk,'DisplayName','TT')
legend()
xlabel('Iteration')
ylabel('Time (seconds)')

subplot(1,3,2)
grid on;
hold on;
plot(fgdata.LocErr,fgmrk,'DisplayName','Full Grid')
plot(ttstdata.LocErr,ttstmrk,'DisplayName','TT')
legend()
set(gca,'YScale','log')
xlabel('Iteration')
ylabel('Norm of Residual')


subplot(1,3,3)
grid on;
hold on;
% plot(fgdata.Itertime,fgmrk,'DisplayName','fg')
plot(ttstcomp,ttstmrk,'DisplayName','TT')
legend()
ylim([5e-3,1]);
set(gca,'YScale','log')
xlabel('Iteration')
ylabel('Compression Ratio')



