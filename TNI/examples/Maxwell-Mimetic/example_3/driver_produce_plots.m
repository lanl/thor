% produce figures
close all; clear; clc;

%% convergence figures: error vs (1) storage, error vs (2) 1/h
ttdata = load('plot_data/tt_example_3.mat');
fdata = load('plot_data/fullgrid_example_3.mat');

%% compute error
Nf = fdata.Ns;
If = 1:numel(Nf);
Errf = fdata.Err(If,1) + fdata.Err(If,3);
MemEzf = ttdata.MemEz(If,2);
Timef = fdata.Time;

Nt = ttdata.Ns;
It = 1:numel(Nt);
Errt = ttdata.Err(It,1)+ ttdata.Err(It,3);
MemEzt(It) = ttdata.MemEz(It,1);
Timet = ttdata.Time(It);
%%
figure(1)
subplot(1,3,1)
plot(Nf, Errf, '-*r');
hold on;
plot(Nt,Errt,'o--b');
set(gca,'YScale','log')
xlabel('Grid size 2^x');
ylabel('|E-Etrue| + |H-Htrue| (log)')
legend('fg','tt')
title('Error vs Grid Size')

subplot(1,3,3)
plot(MemEzf, Errf, '-*r');
hold on;
plot(MemEzt,Errt,'o--b');
set(gca,'XScale','log','YScale','log')
xlabel('Number of Elements used in storage (log)');
ylabel('|E-Etrue| + |H-Htrue| (log)')
legend('fg','tt')
title('Error vs Storage')

subplot(1,3,2)
plot(Timef/3600, Errf, '-*r');
hold on;
plot(Timet/3600, Errt, 'o--b');
set(gca,'XScale','log','YScale','log')
ylabel('|E-Etrue| + |H-Htrue| (log)');
xlabel('Time in hours (log)')
legend('fg','tt')
title('Elapsed Time vs Error')
