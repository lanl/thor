%> @file dtt_amen.m
%!!--------------------------------------------------------------------------~*
%!! Copyright (c) 2024 Triad National Security, LLC
%!! All rights reserved.
%!!--------------------------------------------------------------------------~*/

%!!
%!! @author Oleg Korobkin (korobkin@lanl.gov)
%!! @date   February 2025
%!! @brief  Benchmarking amen_mv and amen_solve tensor generator
%!!
global timers;

batch_size = 10;
epsi = 1e-8;
low_rank = 5;
d = 4;
nswp = 1;

nn = zeros(1, d);
qq = zeros(1, d);
ss = zeros(1, d);
rr = ones(1, d+1);
rr(2:d) = low_rank; % (1, 5, 5, 5, 1)
timers = zeros(1, 10);

fprintf('# MATLAB test: amen_mv and amen_solve\n');
fprintf('# Tensor order = %d, low rank = %d\n', d, low_rank);
fprintf('# Modes = (M, M, ..., M)\n');
fprintf('# Rounding tolerance: epsi = %e\n',epsi);
fprintf('# 1:mode_size 2:amen-mv[s] 3:amen-solve[s]\n');

mode_size = 2;
for i = 1:12
    qq(1:d) = mode_size;
    ss(1:d) = mode_size;
    nn = qq.*ss;

    x = tt_rand(qq, d, low_rank);
    tA = tt_rand(nn, d, low_rank);
    A = tt_matrix(tA, ss, qq);

    tic;
    timers(1:10)= 0.0;
    for j = 1:batch_size
        b = amen_mv(A, x, epsi, 'nswp', nswp, 'verb', 0);
    end
    t_amen_mv = toc/batch_size;

    tic;
    for j = 1:batch_size
        x = fort_amen_solve(A, b, epsi, 'nswp', nswp, 'verb', 0);
    end
    t_amen_solve = toc/batch_size;
    fprintf('%10d  %12.5e %12.5e\n',  mode_size, t_amen_mv, t_amen_solve);
    %timers(1:10)= timers(1:10)/batch_size;
    %fprintf('%10d  %12.5e  %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e\n', ...
    %        mode_size, sum(timers(1:10)), ...
    %        timers(1),timers(2),timers(3),timers(4),timers(5),...
    %        timers(6),timers(7),timers(8),timers(9),timers(10));

    mode_size = mode_size * 2;
end
