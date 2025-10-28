format long g;

fprintf('\n- reading A2 from file A2.dat:\n');
A2 = tt_read_ascii('/tmp/A2.dat')

fprintf('\n- reading A2r from file A2r.dat:\n');
A2r = tt_read_ascii('/tmp/A2r.dat')

fprintf('\n- computing S = A2r - A2;\n');
S = A2r - A2;

fprintf('\n- norm of S:\n');
norm(S)

fprintf('\n- reading S from file S.dat:\n');
S = tt_read_ascii('/tmp/S.dat')

fprintf('\n- norm of S from the file:\n');
norm(S)

