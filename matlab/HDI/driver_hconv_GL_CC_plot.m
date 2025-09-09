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

%% load data
CC2 = load('hconv_test_result_Genz2_CC.mat');
CC = load('hconv_test_result_Genz45_CC.mat');
GL = load('hconv_test_result_Genz245_GL.mat');
CC.R(2,:,:) = CC2.R(2,:,:);

%%
Itest = [2,4,5];
Ds = CC.Ds;
Ns = CC.Ns;

figure()
Markers = {'+','>','*','x','v','d','^','s','o','<'};

for ii = 1:3
  
  subplot(3,2,2*(ii-1)+1)
  GLErr = squeeze(cellfun(@(c) c.Err, GL.R(Itest(ii),:,:)));
  for i = 1:numel(Ds)
    hold on;
    %compute the slope
    if ii==1
      label{1,i} = sprintf('d=%d', Ds(i));
    end
    plot(Ns, GLErr(i,:),strcat('-',Markers{i}));
  end
  set(gca,'XScale','log','YSCale','log');
  if ii==1
    legend(label);
  end
  xlabel('Number of Intervals');
  ylabel('Approximation Error');
  
  %plot the triangle
  slope = log(GLErr(i,4)/GLErr(i,3))/log(1/2)
  
  subplot(3,2,2*ii)
  CCErr = squeeze(cellfun(@(c) c.Err, CC.R(Itest(ii),:,:)));
  for i = 1:numel(Ds)
    hold on;
    %compute the slope
    plot(Ns, CCErr(i,:),strcat('-',Markers{i}));
  end
  set(gca,'XScale','log','YSCale','log');
  xlabel('Number of Intervals');
  ylabel('Approximationg Error');
  
  %plot the triangle
  slope = log(CCErr(i,4)/CCErr(i,3))/log(1/2)
  
end
