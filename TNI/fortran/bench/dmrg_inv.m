%> @file dtt_amen.m
%!!--------------------------------------------------------------------------~*
%!! Copyright (c) 2024 Triad National Security, LLC
%!! All rights reserved.
%!!--------------------------------------------------------------------------~*/

%!!
%!! @author Oleg Korobkin (korobkin@lanl.gov)
%!! @date   February 2025
%!! @brief  Benchmarking tensor inversion using DMRG
%!!

batch_size = 1;
epsi = 1e-8;
low_rank = 5;
d = 4;

fprintf('# MATLAB test: inversion with DMRG\n');
fprintf('# Tensor order = %d, low rank = %d\n', d, low_rank);
fprintf('# Modes = (M, M, ..., M)\n');
fprintf('# Rounding tolerance: epsi = %e\n',epsi);
fprintf('# 1:mode_size 2:dmrg_cross[s] 3:greedy2_cross[s]\n');

mode_size = 4;
for i = 1:14

    x = tt_rand(mode_size, d, low_rank);
    f = @(ind) 1.0/x(ind);
    g = @(ind) sum((ind(1:d)/mode_size).^2);

    tic;
    %for j = 1:batch_size
    %    y = dmrg_cross(d, mode_size, f, epsi, 'verb', 0);
    %end
    t_dmrg_inv = toc/batch_size;
    
    tic;
    %for j = 1:batch_size
    %    y = greedy2_cross(x.n, f, epsi, 'verb', 0);
    %end
    t_dmrgg_inv = toc/batch_size;

    tic;
    for j = 1:batch_size
        y = greedy2_cross(x.n, g, epsi, 'verb', 0);
    end
    t_dmrgg_fun = toc/batch_size;
    fprintf('%10d  %12.5e %12.5e %12.5e  %4d %4d %4d %4d\n', ...
            mode_size, t_dmrg_inv, t_dmrgg_inv, t_dmrgg_fun, ...
            y.r(1), y.r(2), y.r(3), y.r(4));
    mode_size = 2*mode_size;
end
