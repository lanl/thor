%> @file dtt_amen.m
%!!--------------------------------------------------------------------------~*
%!! Copyright (c) 2024 Triad National Security, LLC
%!! All rights reserved.
%!!--------------------------------------------------------------------------~*/

%!!
%!! @author Oleg Korobkin (korobkin@lanl.gov)
%!! @date   February 2025
%!! @brief  Benchmarking amen_mm: matrix x matrix multiplier
%!!
global timers;

batch_size = 1;%0;
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

fprintf('# MATLAB test: amen_mm\n');
fprintf('# Tensor order = %d, low rank = %d\n', d, low_rank);
fprintf('# Modes = (MxM, MxM, ..., MxM)\n');
fprintf('# Rounding tolerance: epsi = %e\n',epsi);
fprintf('# 1:mode_size 2:amen-mm[s]\n');

mode_size = 128;
for i = 1:1%2
    qq(1:d) = mode_size;
    ss(1:d) = mode_size;
    nn = qq.*ss;

    tA = tt_rand(nn, d, low_rank);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fid = fopen('/tmp/A.dat','w');              %
for j=1:d                                   %
    fprintf(fid, '%d,', qq(j));             %
end                                         %
fprintf(fid, '\n');                         %
for j=1:d                                   %
    fprintf(fid, '%d,', ss(j));             %
end                                         %
fprintf(fid, '\n');                         %
for j=1:d+1                                 %
    fprintf(fid, '%d,', rr(j));             %
end                                         %
fprintf(fid, '\n');                         %
for j=1:mem(tA)                             %
    fprintf(fid, '%14.7e\n', tA.core(j));   %
end                                         %
fclose(fid);                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    A = tt_matrix(tA, ss, qq);
    tB = tt_rand(nn, d, low_rank);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fid = fopen('/tmp/B.dat','w');              %
for j=1:d                                   %
    fprintf(fid, '%d,', qq(j));             %
end                                         %
fprintf(fid, '\n');                         %
for j=1:d                                   %
    fprintf(fid, '%d,', ss(j));             %
end                                         %
fprintf(fid, '\n');                         %
for j=1:d+1                                 %
    fprintf(fid, '%d,', rr(j));             %
end                                         %
fprintf(fid, '\n');                         %
for j=1:mem(tB)                             %
    fprintf(fid, '%14.7e\n', tB.core(j));   %
end                                         %
fclose(fid);                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    B = tt_matrix(tB, ss, qq);
    tX = tt_rand(nn, d, low_rank);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fid = fopen('/tmp/X.dat','w');              %
for j=1:d                                   %
    fprintf(fid, '%d,', qq(j));             %
end                                         %
fprintf(fid, '\n');                         %
for j=1:d                                   %
    fprintf(fid, '%d,', ss(j));             %
end                                         %
fprintf(fid, '\n');                         %
for j=1:d+1                                 %
    fprintf(fid, '%d,', rr(j));             %
end                                         %
fprintf(fid, '\n');                         %
for j=1:mem(tX)                             %
    fprintf(fid, '%14.7e\n', tX.core(j));   %
end                                         %
fclose(fid);                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    X = tt_matrix(tX, ss, qq);
%DEBUG% fid = fopen('bench/A.dat','r');
%DEBUG% qq = fscanf(fid, '%d,', d);
%DEBUG% ss = fscanf(fid, '%d,', d);
%DEBUG% rr = fscanf(fid, '%d,', d+1);
%DEBUG% ta = tt_rand(qq.*ss, d, rr);
%DEBUG% ta.core = fscanf(fid, '%e', mem(ta));
%DEBUG% fclose(fid);
%DEBUG% A = tt_matrix(ta, qq, ss)
%DEBUG% 
%DEBUG% fid = fopen('bench/B.dat','r');
%DEBUG% qq = fscanf(fid, '%d,', d);
%DEBUG% ss = fscanf(fid, '%d,', d);
%DEBUG% rr = fscanf(fid, '%d,', d+1);
%DEBUG% tb = tt_rand(qq.*ss, d, rr);
%DEBUG% tb.core = fscanf(fid, '%e', mem(tb));
%DEBUG% fclose(fid);
%DEBUG% B = tt_matrix(tb, qq, ss)
%DEBUG% 
%DEBUG% fid = fopen('bench/X.dat','r');
%DEBUG% qq = fscanf(fid, '%d,', d);
%DEBUG% ss = fscanf(fid, '%d,', d);
%DEBUG% rr = fscanf(fid, '%d,', d+1);
%DEBUG% tx = tt_rand(qq.*ss, d, rr);
%DEBUG% tx.core = fscanf(fid, '%e', mem(tx));
%DEBUG% fclose(fid);
%DEBUG% X = tt_matrix(tx, qq, ss)

    tic;
    for j = 1:batch_size
%DEBUG% C = amen_mm(A, B, epsi, 'nswp', nswp, 'verb', 4);
C = amen_mm(A, B, epsi, 'nswp', nswp, 'x0', X, 'verb', 4);
    end
    t_amen_mm = toc/batch_size;

    fprintf('%10d  %12.5e\n',  mode_size, t_amen_mm);
    mode_size = mode_size * 2;
end
