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

fprintf('# MATLAB: matmul operation\n')
fprintf('# 1:mode_size 2:matmul[s] 3:matmul-transposed[s]\n')

batch_size = 100;

for mode_size = 1000:1000:10000
    A = rand(mode_size, mode_size);
    B = rand(mode_size, mode_size/100);
    tic;
    for test = 1:batch_size
        C = A*B;
    end
    tm_matmul = toc/batch_size;
    tic;
    for test = 1:batch_size
        C = A'*B;
    end
    tm_matmul_tp = toc/batch_size;
    fprintf('%10d  %12.5e  %12.5e\n', mode_size, tm_matmul, tm_matmul_tp);
end
