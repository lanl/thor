close all; clear; clc;
% this script will produce plots for experiment sections

%% setup plotting parameter here -- to be consistent accross machines
set(0,'defaultfigurecolor','w');
set(0,'DefaultAxesTitleFontWeight','normal');
set(0,'DefaultAxesFontSize',15)
set(0,'DefaultAxesFontName','Times New Roman');
% set(0,'DefaultLineMarkerSize',10)
set(0,'defaultLineLineWidth',3);
set(0,'defaultTextFontName', 'Time New Roman')
set(0, 'DefaultFigureRenderer', 'painters');

%% What to plot
itest = 5; % dataset to plot
plotacc = 1; % plot accuracy figure
plotMC = 0; % plot MC
plotadapt = 0; %plot adaptivity for ANOVA

%% Accuracy ~~~ only convergence of Genz 2,4,5 can be visualized
if plotacc
  %% h-convergence plots
  Trap = load('./Results/Acc_test_result_Trapezoidal.mat');
  %~~~~~~~~~~~~~~~~~~~ Trapezoidal
  TrapErr = squeeze(cellfun(@(c) c.Err, Trap.R(itest,:,:)));
  TrapDs = Trap.Ds;
  TrapNs = Trap.Ns;
  
  figure()
  Markers = {'+','>','*','x','v','d','^','s','o','<'};
  %%
  subplot(2,2,1)
  for i = 1:numel(TrapDs)
    hold on;
    %compute the slope
    temp = polyfit(log(TrapNs), log(TrapErr(i,:)), 1);
    label{1,i} = sprintf('d=%d',TrapDs(i));
    plot(TrapNs,TrapErr(i,:),strcat('-',Markers{i}));
    %   plot(log(TrapNs),log(TrapErr(i,:)));
  end
  set(gca,'XScale','log','YSCale','log');
  legend(label);
  xlabel('number of intervals');
  ylabel('Error');
%   title(sprintf('Trapezoidal - Genz %d function',itest));
  
  %~~~~~~~~~~~~~~~~~~~ Simpson
  Simp = load('./Results/Acc_test_result_Simpson.mat');
  SimpErr = squeeze(cellfun(@(c) c.Err, Simp.R(itest,:,:)));
  SimpDs = Simp.Ds;
  SimpNs = Simp.Ns;
  
  %%
  subplot(2,2,2)
  for i = 1:numel(SimpDs)
    hold on;
    %compute the slope
    temp = polyfit(log(SimpNs), log(SimpErr(i,:)), 1);
    label{1,i} = sprintf('d=%d', SimpDs(i));
    plot(SimpNs, SimpErr(i,:),strcat('-',Markers{i}));
    %   plot(log(TrapNs),log(TrapErr(i,:)));
  end
  set(gca,'XScale','log','YSCale','log');
  legend(label);
  xlabel('number of intervals');
  ylabel('Error');
%   title(sprintf('Simpson - Genz %d function',itest));
  
  %% p-convergence plots
  % ~~~~~~~~~~ Gauss Legendre
  G = load('Results/acc_pconv_GL.mat');
  GDs = G.Ds(1:5);
  GNq = G.Nq;
  GErr = squeeze(cellfun(@(c) c.Err, G.R(itest,:,:)));
  subplot(2,2,3)
  for i = 1:numel(GDs)
    hold on;
    plot(GNq(),GErr(i,:),strcat('-',Markers{i}));
    label{1,i} = sprintf('d=%d',GDs(i));
    xlabel('Nq');
    ylabel('Error');
  end
  legend(label)
  set(gca,'YScale','log');
%   sgtitle(sprintf('GL - Genz %d func',itest));
  
  % ~~~~~~~~~~ Gauss Legendre
  C = load('Results/acc_pconv_CC.mat');
  GDs = C.Ds;
  GNq = C.Nq;
  GErr = squeeze(cellfun(@(c) c.Err, C.R(itest,:,:)));
  subplot(2,2,4)
  for i = 1:numel(GDs)
    hold on;
    plot(GNq,GErr(i,:),strcat('-',Markers{i}));
    label{1,i} = sprintf('d=%d',GDs(i));
    xlabel('Nq');
    ylabel('Error');
  end
  legend(label)
  set(gca,'YScale','log');
%   sgtitle(sprintf('CC - Genz %d func', itest));
end

%% Efficiency - MC
% compare the performace with GL p-convergence result
% plot err vs time
if plotMC
  G = load('Results/acc_pconv_GL.mat');
  L = load('Results/MC_cubLattice_result.mat');
  S = load('Results/MC_cubSobol_result.mat');
  Ds = L.Ds;
  
  cM = jet(numel(Ds));
  figure()
  
  for i = 1:numel(Ds)
    plot(nan,nan,'-','Color',cM(i,:));
  end
  
  for i = 1:numel(Ds)
    %get data
    subplot(3,2,i)
    hold on;
    
    Gerr = squeeze(cellfun(@(c) c.Err, G.R(itest,i,:)));
    Gtime = squeeze(cellfun(@(c) c.time,G.R(itest,i,:)));
    Lerr = squeeze(cellfun(@(c) c.Err,L.R(itest,i,:)));
    Ltime = squeeze(cellfun(@(c) c.time,L.R(itest,i,:)));
    Serr = squeeze(cellfun(@(c) c.Err,S.R(itest,i,:)));
    Stime = squeeze(cellfun(@(c) c.time,S.R(itest,i,:)));
    
    %plot
    plot(Gtime(:), Gerr(:),'-*','MarkerSize',10);
    plot(Ltime(:), Lerr(:),'-*','MarkerSize',10);
    plot(Stime(:), Serr(:),'-*','MarkerSize',10);
    if i == numel(Ds)
      legend({'tt-GL','MCLattice','MCSobol'});
    end
    title(sprintf('d=%d',Ds(i)));
    set(gca,'YScale','log','XScale','log');
    xlabel('time (s)');
    ylabel('error')
  end
%   sgtitle(sprintf('Genz %d function',itest));
  
end

%% Adaptivity - Anova
if plotadapt
  %% Anova 1/2
  a1 = load('./Results/adapt_12_results.mat');
  Ds = a1.Ds;
  Ns = a1.Ns;
  figure()
  subplot(1,2,1)
  for i = 1:numel(Ns)
    hold on;
    err = cellfun(@(c) c.Err, a1.R(1,:,i));
    plot(Ds,err(:),'-*');
    xlabel('number of dimension');
    ylabel('Error');
  end
  legend(string(num2cell(Ns)))
  set(gca,'YScale','log');
  title(sprintf('ANOVA-1/2 Func'));
  %% Anova 3/4
  a1 = load('./Results/adapt_34_results.mat');
  Ds = a1.Ds;
  Ns = a1.Ns;
%   figure()
  subplot(1,2,2)
  for i = 1:numel(Ns)
    hold on;
    err = cellfun(@(c) c.Err, a1.R(1,:,i));
    plot(Ds,err(:),'-*');
    xlabel('number of dimension');
    ylabel('Error');
  end
  legend(string(num2cell(Ns)))
  set(gca,'YScale','log');
  title(sprintf('ANOVA-3/4 Func'));
end

%% Subroutines


