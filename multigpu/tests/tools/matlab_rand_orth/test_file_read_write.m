fprintf('\n- writing a random TT-tensor to /tmp/A.dat;\n');
d = 4; n = 3; r = 2;
A = tt_rand(n, d, r, 1)
tt_write_ascii(A, '/tmp/A.dat');

fprintf('- reading from /tmp/A.dat:\n');
B = tt_read_ascii('/tmp/A.dat')

fprintf('- tensor diff = %14.7e\n', norm(A - B));
fprintf('- flat core diff = %14.7e\n', norm(A.core - B.core));

