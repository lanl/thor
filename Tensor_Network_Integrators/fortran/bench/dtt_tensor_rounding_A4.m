%> @file dtt_tensori_rounding_A4.m
%!!--------------------------------------------------------------------------~*
%!! Copyright (c) 2024 Triad National Security, LLC
%!! All rights reserved.
%!!--------------------------------------------------------------------------~*/

%!!
%!! @author Oleg Korobkin (korobkin@lanl.gov)
%!! @date   February 2025
%!! @brief  Benchmarking dtt tensor rounding operation
%!!

batch_size = 10;
epsi = 1e-8;
low_rank = 5;

fprintf('# MATLAB test: round a random low-rank 4D tensor\n');
fprintf('# Modes = (M, 2M, 3M, 4M), low_rank = %d\n', low_rank);
fprintf('# Rounding tolerance: epsi = %e\n',epsi);
fprintf('# 1:mode_size 2:time[s] 3-7:ranks\n');

mode_size = 125;
for i = 1:15
    nn = [mode_size, 2*mode_size, 3*mode_size, 4*mode_size];
    tt = tt_rand(nn, size(nn(:),1), low_rank);
    tic;
    for j = 1:batch_size
        tt1 = round(tt, epsi);
    end
    tm = toc;
    fprintf('%10d  %12.5e  %4d %4d %4d %4d %4d\n', ...
            mode_size, tm/batch_size, ...
            tt1.r(1), tt1.r(2), tt1.r(3), tt1.r(4), tt1.r(5));
    mode_size = mode_size*2;
end
