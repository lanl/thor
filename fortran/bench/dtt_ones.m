%> @file dtt_ones.m
%!!--------------------------------------------------------------------------~*
%!! Copyright (c) 2024 Triad National Security, LLC
%!! All rights reserved.
%!!--------------------------------------------------------------------------~*/

%!!
%!! @author Oleg Korobkin (korobkin@lanl.gov)
%!! @date   February 2025
%!! @brief  Benchmarking "tt_ones" tensor generator
%!!

fprintf('# MATLAB test: tt_ones, mode size = 100\n');
fprintf('# 1:dimension 2:time[s]\n');

batch_size = 1000;
mode_size  = 100;
d = 2;
for i = 1:11
    tic;
    for j = 1:batch_size
        v = tt_ones(mode_size, d);
    end
    tm = toc;
    fprintf('%10d  %12.5e  %16d\n', d, tm/batch_size, mem(v));
    d = d * 2;
end
