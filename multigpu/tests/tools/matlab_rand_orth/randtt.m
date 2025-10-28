format long g;
%close all; clear;
%WORKDIR = pwd();
%TTTBXDIR = "../../TT-Toolbox/";
%chdir(TTTBXDIR);
%setup;
%chdir(WORKDIR);

% S = tt_read_ascii('S.dat');
% norm(S)
% n = S.n(1); 
% d = size(S.n, 1);
% R = tt_rand(n, d, 1, 1);
% R.core(:) = 1.5;
% S1 = rounding_rand_orth(S, R, 1);
% S1
% S1.core
% exit

% testing the unfoldings
r1 = 2; r2 = 3; nk = 4;
A = zeros(r1,nk,r2);
for i = 1:r1
    A(i,:,:) = i*100;
end
for j = 1:nk
    A(:,j,:) = A(:,j,:) + j*10;
end
for k = 1:r2
    A(:,:,k) = A(:,:,k) + k;
end
Aflat = A(:);
Av = reshape(A, [r1*nk,r2]); % v-unfolding
Ah = reshape(permute(A, [1,3,2]), [r1,nk*r2]); % h-unfolding

%% % Create a random tt-tensor
%% n = 8; d = 4; r = 2; rnew = 2;
%% a = tt_rand(n, d, r, 1);
%% a256 = 256*a;
%% aaaa = a + a + a + a; % x4
%% aaaa = aaaa + aaaa + aaaa + aaaa; % x 16
%% aaaa = aaaa + aaaa + aaaa + aaaa; % x 64
%% aaaa = aaaa + aaaa + aaaa + aaaa; % x256
%% b = rounding_rand_orth(aaaa, r);
%% norm(a256 - b)

%% % Create a random tt-tensor
%% n = 4; d = 3; r = 2; rnew = 1;
%% a = tt_rand(n, d, r, 1);
%% b = rounding_rand_orth(a - a, rnew);
%% b.core

fprintf('\n- reading A2 from file A2.dat:\n')
A2 = tt_read_ascii('A2.dat')

fprintf('\n- reading A2r from file A2r.dat:\n')
A2r = tt_read_ascii('A2r.dat')

fprintf('\n- computing norm(S) = norm(A2r + A2):\n')
S = A2r + A2;
norm(S)

fprintf('\n- reading S from file S.dat:\n')
S = tt_read_ascii('S.dat')
norm(S)

exit

% read A from file, A.dat
d = 4;
fid = fopen('Athor.dat','r');
nn = fscanf(fid, '%d,', d);   n = nn(1);
rr = fscanf(fid, '%d,', d+1); r = rr(2);
A = tt_rand(nn, d, rr)
A.core = fscanf(fid, '%e', mem(A));
fclose(fid);

% test with A1: round A(ranl=1)
R = tt_rand(n, d, 1, 1);
R.core(:) = 1.5;
%A1 = rounding_rand_orth(A, R, 1);
%A1
%A1.core

% Test with A2 = A + A
A2 = A + A;
A2r = rounding_rand_orth(A2, A, r);
A2r.core
S = A2r - 2*A;
S.core
%% R = tt_rand(n, d, 1, 1);
%% R.core(:) = 1.5;
%% Sr = rounding_rand_orth(S, R, 1);
%% Sr.core
%% norm(Sr)
