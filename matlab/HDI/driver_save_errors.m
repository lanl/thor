close all; clear;clc;

Trap = load('./Results/Acc_test_result_Trapezoidal.mat');
Simp = load('./Results/Acc_test_result_Simpson.mat');
G = load('Results/acc_pconv_GL.mat');
C = load('Results/acc_pconv_CC.mat');

%%
itest=2;
GErr = squeeze(cellfun(@(c) c.Err, C.R(itest,:,:)));


