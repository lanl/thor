%%
close all; clear; clc;
run setup.m;
 A = [1 2 3; 4 5 6; 7 8 9; 10 11 12; 13 14 15];
t = tt_tensor(A);
% t{1} = t{1}(:,2:end,:); %this is how we remove first row
%%
G = core2cell(t);
% G{1}(:,1,:)=[0,0]';
G3 = zeros(1,6,2);
G3(:,2:end,:)= G{1};
G{1} = G3;
%%
t = cell2core(t,G);
%%
b = full(t,[6,3])