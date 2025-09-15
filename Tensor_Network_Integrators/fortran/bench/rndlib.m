%> @file b_randn.m
%!!--------------------------------------------------------------------------~*
%!! Copyright (c) 2025 Triad National Security, LLC
%!! All rights reserved.
%!!--------------------------------------------------------------------------~*/

%!!
%!! @author Oleg Korobkin (korobkin@lanl.gov)
%!! @date   July 2025
%!! @brief  Benchmarking "randn" function in matlab
%!!

fprintf('# MATLAB test: rand & randn, array size = 10^2..10^8\n');
fprintf('# 1:size 2:time[s]\n');

batch_size = 100;
n  = 100;
for i = 2:8
    tic;
    for j = 1:batch_size
        v = rand(n, 1);
    end
    dt_rand = toc/batch_size;

    tic;
    for j = 1:batch_size
        v = randn(n, 1);
    end
    dt_randn = toc/batch_size;
    
    fprintf('%10d  %12.5e  %12.5e\n', n, dt_rand, dt_randn);
    n = n * 10;
end
