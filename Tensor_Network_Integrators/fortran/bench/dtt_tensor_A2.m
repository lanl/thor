%> @file dtt_tensor_A2.m
%!!--------------------------------------------------------------------------~*
%!! Copyright (c) 2024 Triad National Security, LLC
%!! All rights reserved.
%!!--------------------------------------------------------------------------~*/

%!!
%!! @author Oleg Korobkin (korobkin@lanl.gov)
%!! @date   February 2025
%!! @brief  Benchmarking "tt_tensor" tensor generator
%!!

batch_size = 10;
epsi = 1e-8;
fprintf('# MATLAB test: tt_tensor(A2,epsi)\n');
fprintf('# Modes = (M, 2M)\n');
fprintf('# Rounding tolerance: epsi = %e\n', epsi);
fprintf('# 1:mode_size 2:time[s] 3-5:ranks\n');

for mode_size = 100:100:1000
    A2 = rand(mode_size, 2*mode_size);
    tic;
    for j = 1:batch_size
        v = tt_tensor(A2, epsi);
    end
    tm = toc;
    fprintf('%10d  %12.5e  %4d %4d %4d\n', ...
            mode_size, tm/batch_size, ...
            v.r(1), v.r(2), v.r(3));
end
