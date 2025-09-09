%%
%% Matlab script to check that ty = mA*tx for tt-tensors tx and ty
%% and tt-matrix mA, written to tx.dat, ty.dat and mA.dat files resp.
%%


% determine the dimensions clumsily
pwd
fid = fopen('tx.dat','r');
d = size(strsplit(fgetl(fid),','),2);
fclose(fid);

% open the file tx.dat
fid = fopen('tx.dat','r');
nx = fscanf(fid, '%d,', d);
rx = fscanf(fid, '%d,', d+1);
tx = tt_rand(nx, d, rx)
tx.core = fscanf(fid, '%e', mem(tx));
fclose(fid);

% open the file ty.dat
fid = fopen('ty.dat','r');
ny = fscanf(fid, '%d,', d);
ry = fscanf(fid, '%d,', d+1);
ty = tt_rand(ny, d, ry)
ty.core = fscanf(fid, '%e', mem(ty));
fclose(fid);

% open the matrix file in tensor form, mA.dat
fid = fopen('mA.dat','r');
nn = fscanf(fid, '%d,', d);
rr = fscanf(fid, '%d,', d+1);
ta = tt_rand(nn, d, rr);
ta.core = fscanf(fid, '%e', mem(ta));
fclose(fid);

% create matrix and use amen_mv to do the matrix-vector multiplication
mA = tt_matrix(ta, ny, nx)
ty1 = amen_mv(mA, tx, 1e-12)

% check the answer
nrm_ty = norm(ty);
nrm_diff = norm(ty - ty1);
relative_diff = nrm_diff/nrm_ty

% check with fort_amen_solve
tx1 = fort_amen_solve(mA, ty, 1e-5)

% check the answer
nrm_tx = norm(tx);
nrm_diff = norm(tx - tx1);
relative_diff = nrm_diff/nrm_tx

