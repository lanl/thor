close all; clear; clc;
run setup1DSlab.m;
%%

ktrue = 1.00334;

partisantime = [14.22,16.67,21.59,32.42];

Ri9 = cell(4,3);

Rdarw = cell(4,2);

Tols = {'1e4', '1e5', '1e6'};
Nxs = {'256','512','1024','2048'};
folder = './Results/I9';
folder2 = './Results/darwin';
%% load file
for i = 1:4
  filename = sprintf('%s/OneD_261E_tt_nx%s_result_1e6.mat',folder2,Nxs{i});
  Rdarw{i,1} = load(filename);
  filename = sprintf('%s/OneD_261E_qtt_nx%s_result_1e6.mat',folder2,Nxs{i});
  Rdarw{i,2} = load(filename);
  for j = 1:3
    filename = sprintf('%s/OneD_261E_qtt_nx%s_result_%s_i9.mat',folder,Nxs{i},Tols{j});
    Ri9{i,j} = load(filename);
  end
end

%% extract results
tttime6dw = cellfun(@(c) sum(c.Itertime),Rdarw(:,1));
tterr6 = cellfun(@(c) abs(c.ttk-ktrue),Rdarw(:,1))
qtttime6dw = cellfun(@(c) sum(c.Itertime),Rdarw(:,2));
qtterr6 = cellfun(@(c) abs(c.ttk-ktrue),Rdarw(:,2))
%% get qtt results
Time = cellfun(@(c) sum(c.Itertime),Ri9);
Err = cellfun(@(c) abs(c.ttk - ktrue),Ri9);
Compress = cellfun(@(c) compress_ratio_tt(c.ttPsi1),Ri9);
%%
Tols = {'1e-4', '1e-5', '1e-6'};
figure(1)
subplot(1,3,1)
hold on;
plot(Time,'-*');
plot(partisantime,'--*r');
ylabel('Time in seconds')
set(gca,'YScale','log')
legend([Tols,{'Partisan at 1e-6'}])
title('Elapsed Time')

% figure(2)
subplot(1,3,2)
hold on;
plot(Err,'-*');
ylabel('Error')
set(gca,'YScale','log')
legend([Tols])
title('Difference in eigenvalue vs Partisan')

% figure(3)
subplot(1,3,3)
hold on;
plot(Compress,'-*');
ylabel('Compress Ratio')
set(gca,'YScale','log')
legend([Tols])
title('Compress ratio of Psi')



