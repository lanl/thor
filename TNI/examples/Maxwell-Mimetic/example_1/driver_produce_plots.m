% produce figures
close all; clear; clc;

%% convergence figures: error vs (1) storage, error vs (2) 1/h
ttdata = load('plot_data/tt_test1_k1_tol_1e-05.mat');
fdata = load('plot_data/fullgrid_test1_k1.mat');

%% compute error
Nf = fdata.Ns;
Errf = fdata.Err(Nf,1) + fdata.Err(Nf,3);
MemEzf = ttdata.MemEz(Nf,2);
Timef = fdata.Time;

Nt = ttdata.Ns;
Errt = ttdata.Err(Nt,1)+ ttdata.Err(Nt,3);
MemEzt(Nt) = ttdata.MemEz(Nt);
Timet = ttdata.Time(Nt);
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
