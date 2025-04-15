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
fid = fopen('/tmp/mC.dat','r');
qq = fscanf(fid, '%d,', d);
ss = fscanf(fid, '%d,', d);
rr = fscanf(fid, '%d,', d+1);
ta = tt_rand(qq.*ss, d, rr);
ta.core = fscanf(fid, '%e', mem(ta));
fclose(fid);
mC = tt_matrix(ta, qq, ss)

% open the matrix file in tensor form, mB.dat
fid = fopen('/tmp/mB.dat','r');
qq = fscanf(fid, '%d,', d);
ss = fscanf(fid, '%d,', d);
rr = fscanf(fid, '%d,', d+1);
ta = tt_rand(qq.*ss, d, rr);
ta.core = fscanf(fid, '%e', mem(ta));
fclose(fid);
mB = tt_matrix(ta, qq, ss)


% test mA @ mB = mC1 == mC
mC1 = amen_mm(mA, mB, 1e-12)

% check the answer
nrm_C = norm(mC);
nrm_diff = norm(mC - mC1);
relative_diff = nrm_diff/nrm_C
