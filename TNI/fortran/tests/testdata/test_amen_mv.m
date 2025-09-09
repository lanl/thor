%%
%% Matlab script to check that mC = mA @ mB for tt-matrices
%%
fid = fopen('/tmp/mA.dat','r');
d = size(strsplit(fgetl(fid),','),2);
fclose(fid);

% open the matrix file in tensor form, mA.dat
fid = fopen('/tmp/mA.dat','r');
qq = fscanf(fid, '%d,', d);
ss = fscanf(fid, '%d,', d);
rr = fscanf(fid, '%d,', d+1);
ta = tt_rand(qq.*ss, d, rr);
ta.core = fscanf(fid, '%e', mem(ta));
fclose(fid);
mA = tt_matrix(ta, qq, ss)

% open the matrix file in tensor form, mC.dat
fid = fopen('/tmp/tx.dat','r');
xn = fscanf(fid, '%d,', d);
rr = fscanf(fid, '%d,', d+1);
tx = tt_rand(xn, d, rr)
tx.core = fscanf(fid, '%e', mem(tx));
fclose(fid);

% open the matrix file in tensor form, mB.dat
fid = fopen('/tmp/ty.dat','r');
yn = fscanf(fid, '%d,', d);
rr = fscanf(fid, '%d,', d+1);
ty = tt_rand(yn, d, rr);
ty.core = fscanf(fid, '%e', mem(ty));
fclose(fid);

% test mA * tx = ty1 == ty
ty1 = amen_mv(mA, tx, 1e-12)

% check the answer
nrm_y = norm(ty1);
nrm_diff = norm(ty - ty1);
relative_diff = nrm_diff/nrm_y
