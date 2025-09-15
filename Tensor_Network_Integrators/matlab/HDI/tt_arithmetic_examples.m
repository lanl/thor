%%
close all; clear; clc;
run setup.m;
% A = tt_tensor(cell2core(rand(10,2),rand(2,10)));
B = tt_tensor(cos(rand(10,10)) + sin(rand(10,10)));
C = tt_tensor(cos(rand(10,10)) + sin(rand(10,10)));

S = round(times(B,C),1e-10,2)
