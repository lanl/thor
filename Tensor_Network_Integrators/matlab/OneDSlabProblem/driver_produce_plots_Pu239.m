close all; clear; clc;
run setup1DSlab.m;
load('./Results/1Dslab_Pu_239_results_darwin.mat');


%%
k = cellfun(@(c) c.k, R);
fgk = cellfun(@(c) c.fgk, R);
ttk = cellfun(@(c) c.ttk, R);

% psi = cellfun(@(c) c.psi, R);
% fgpsi = cellfun(@(c) c.fgpsi, R);
% ttpsi = cellfun(@(c) c.ttpsi, R);

% figure('Position',[314 818 1246 361])
meth = {'GES', 'ISFM', 'ISTT'};
if 0
  figure()
  %%
  % subplot(1,3,1)
  tiledlayout(1,3);

  nexttile;
  hold on;
  p1 = plot(Nl, k, '--k', 'DisplayName', meth{1});
  p2 = plot(Nl, fgk, '.k', 'DisplayName', meth{2});
  p3 = plot(Nl, ttk, 'dk', 'DisplayName', meth{3});
  % set(gcf, 'YScale', 'log')
  yline(1,'LineWidth',1.0);
  ylabel('Eigenvalues')
  xlabel( 'number of ordinate (angle) L')
  legend([p1(1),p2(1),p3(1)]);
  title('Convergence to ground truth k=1')
  %% Difference in Solutions
  kfgdiff = abs(k-fgk);
  kttdiff = abs(k-ttk);
  kfgttdiff = abs(fgk-ttk);

  % subplot(1,3,2);
  nexttile;
  hold on;
  plot(Nl, kfgdiff, '-k', 'linewidth',2.0,'DisplayName',sprintf('|k_{%s}-k_{%s}|',meth{1}, meth{2}));
  plot(Nl, kttdiff, '--k', 'linewidth',2.0,'DisplayName',sprintf('|k_{%s}-k_{%s}|',meth{1}, meth{3}));
  plot(Nl, kfgttdiff, 'ok', 'linewidth',2.0,'DisplayName',sprintf('|k_{%s}-k_{%s}|',meth{2}, meth{3}));
  set(gca, 'YScale','log')
  % plot(Ns, Psidiff, 's-r', 'DisplayName','Psi')
  xlabel( 'number of ordinate (angle) L')
  ylabel('Error')
  title('Differences in Eigenvalues')

  legend('location','best');

  %% differences in Psi

  psifgdiff = cellfun(@(c) norm(sort(c.psi(:))-sort(c.fgPsi(:)))...
    /norm(sort(c.psi(:))), R);
  psittdiff = cellfun(@(c) norm(sort(c.psi(:))-sort(full(c.ttPsi1)))...
    /norm(sort(c.psi(:))), R);
  psifgttdiff = cellfun(@(c) norm(sort(c.fgPsi(:))-sort(full(c.ttPsi1)))...
    /norm(sort(c.psi(:))), R);

  % subplot(1,3,3);
  nexttile;
  hold on;
  plot(Nl, psifgdiff, '-k', 'linewidth',2.0,'DisplayName', sprintf('|\\Psi_{%s}-\\Psi_{%s}|', meth{1}, meth{2}));
  plot(Nl, psittdiff, 'ok', 'linewidth',2.0,'DisplayName', sprintf('|\\Psi_{%s}-\\Psi_{%s}|', meth{1}, meth{3}));
  plot(Nl, psifgttdiff, '--k', 'linewidth',2.0,'DisplayName', sprintf('|\\exPsi_{%s}-\\Psi_{%s}|', meth{2}, meth{3}));
  set(gca, 'YScale','log')
  % plot(Ns, Psidiff, 's-r', 'DisplayName','Psi')
  xlabel( 'number of ordinate (angle) L')
  ylabel('Relative Error')
  title('Differences in Eigenvectors')

  legend('location','best');

end

%% Elapsed Time

if 1
  figure()

  time = cellfun(@(c) c.time, R);
  fgtime = cellfun(@(c) c.fg_time, R);
  tttime = cellfun(@(c) c.tt_time, R);

  % subplot(1,2,1);
  tiledlayout(1,2);

  nexttile;
  hold on;
  plot(Nl, time, '-k', 'linewidth',2.0,'DisplayName', meth{1});
  plot(Nl, fgtime, 'ok', 'linewidth',2.0,'DisplayName', meth{2});
  plot(Nl, tttime, '--k', 'linewidth',2.0,'DisplayName', meth{3});
  xlabel( 'number of ordinate (angle) L');
  ylabel('Elapsed Time');
  set(gca, 'YScale','log');
  legend('location','best');


  %% Compressed ratio
  Storage = cellfun(@(c) compress_ratio_tt(round(c.ttPsi1,1e-3)),R);
  % subplot(1,2,2)
  nexttile
  plot(Nl,Storage, '*--k','linewidth',2.0)
  xlabel( 'number of ordinate (angle) L')
  ylabel('Compression Ratio')
  set(gca,'YScale','log','box','off')
  % title('Compressed Ratio of eigenvector in QTT')
end