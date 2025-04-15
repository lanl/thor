%> @file dtt_tensor_A4.m
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
fprintf('# MATLAB test: tt_tensor(A4,epsi)\n');
fprintf('# Modes = (M, 2M, 3M, 4M)\n');
fprintf('# Rounding tolerance: epsi = %e\n',epsi);
fprintf('# 1:mode_size 2:time[s] 3-7:ranks\n');

for mode_size = 10:5:50
    A4 = rand(mode_size, 2*mode_size, 3*mode_size, 4*mode_size);
    tic;
    for j = 1:batch_size
        v = tt_tensor(A4, epsi);
    end
    tm = toc;
    fprintf('%10d  %12.5e  %4d %4d %4d %4d %4d\n', ...
            mode_size, tm/batch_size, ...
            v.r(1), v.r(2), v.r(3), v.r(4), v.r(5));
end
