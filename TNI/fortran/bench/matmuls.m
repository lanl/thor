%> @file matmuls.m
%!!--------------------------------------------------------------------------~*
%!! Copyright (c) 2025 Triad National Security, LLC
%!! All rights reserved.
%!!--------------------------------------------------------------------------~*/

%!!
%!! @author Oleg Korobkin (korobkin@lanl.gov)
%!! @date   February 2025
%!! @brief  Testing matrix multiplication in Matlab
%!!

fprintf('# MATLAB: matrix multiplication timing\n')
fprintf('# 1:matrix size 2:matmul[s] 3:transpose[s]\n')

batch_size = 100;
mode_size  = 10;
for i = 1:9
    A = randn(mode_size, 2*mode_size);
    B = randn(2*mode_size, 3*mode_size);
    tic;
    for j = 1:batch_size
        C = A*B;
    end
    dt_matmul = toc/batch_size;
    tic;

    for j = 1:batch_size
        C = C';
    end
    dt_transpose = toc/batch_size;
    fprintf('%10d  %12.5e  %16d\n', mode_size, dt_matmul, dt_transpose);
    mode_size= 2*mode_size;
end
