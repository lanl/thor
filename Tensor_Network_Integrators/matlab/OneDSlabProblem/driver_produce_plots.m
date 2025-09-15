close all; clear; clc;
run setup1DSlab.m;
load('./Results/1Dslab_EVP_tt.mat');

%% Elapsed Time
fgtime = cellfun(@(c) c.fg_time, R);
tttime = cellfun(@(c) c.tt_time, R);

figure(1)
hold on;
plot(Ns, fgtime, 'o-', 'DisplayName','fullgrid');
plot(Ns, tttime, 's-r', 'DisplayName','TT')
xlabel('Psi size - n x n x 8')
ylabel('Elapsed Time')
set(gca, 'YScale','log')
legend;

%% Difference in Solutions
kdiff = cellfun(@(c) norm(c.k -c.ttk)/norm(c.k), R);
Psidiff = cellfun(@(c) norm(c.Psi1 - full(c.ttPsi1))/norm(c.Psi1), R);

figure(2)
hold on;
plot(Ns, kdiff, 'o-', 'DisplayName','eigenvalue');
plot(Ns, Psidiff, 's-r', 'DisplayName','Psi')
xlabel('size of Psi - n x n x 8')
ylabel('Relative Norm')
title('Relative Difference between fg and tt solutions')
% set(gca, 'YScale','log')
legend;

%% Compressed ratio

Storage = cellfun(@(c) compress_ratio_tt(c.ttPsi1),R);
figure(3)
plot(Ns,Storage, 'o-')
xlabel('Psi size - n x n x 8')
ylabel('Compressed Ratio')
set(gca,'YScale','log')
title('Compressed Ratio of ttPsi')