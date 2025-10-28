format long g;

fprintf('\n- reading A from file AA.dat:\n');
AA = tt_read_ascii('AA.dat')

fprintf('\n- reading B from file BB.dat:\n');
BB = tt_read_ascii('BB.dat')

fprintf('\n- computing S = AA + BB;\n');
S = AA + BB;

fprintf('\n- reading AB from file AB.dat:\n');
AB = tt_read_ascii('AB.dat')

fprintf('\n- norm of the differences between flattened cores:\n');
norm(AB.core - S.core)

