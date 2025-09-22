clear; close all; clc;

addpath(genpath('../../matlab/utils/tt-toolbox/'));
addpath(genpath('../../matlab/utils/ttfunc/'));
addpath(genpath('../../matlab/tt-iga/src/'));


%% Geometry & solver parameters
R        = 3;
r_out    = 1.0;
r_int    = 0.5;
p1 = 2; p2 = 2; p3 = 2;
bc_type   = "3";
g_type    = "1";
u_in      = 0.0;
u_out     = 0.0;
ffunction = "sinpxsinpysinpz";
tt_tol    = 5e-10;

%% List of refinements
nm_list = [2, 4, 8, 16];
nRun    = numel(nm_list);

% output file
outFile   = "Torus.csv";
if isfile(outFile)
  delete(outFile);
end

for idx = 1:nRun
    fprintf("===========================================================\n");
    nm = nm_list(idx);
    fprintf("Row %d/%d: \n", idx, nRun);
    % build geometry
    [wn, ctpxn, ctpyn, ctpzn, knot1n, knot2n, knot3n] = ...
        Torus222(R, r_int, r_out, nm, nm, nm);
    wn = ones(size(wn));  % B-spline weights

    % record grid size
    [n1, n2, n3] = size(ctpxn);
    fprintf("Torus %d/%d: ", idx, nRun);
    fprintf("Size %d x %d x %d: \n", n1, n2, n3);
    % TT solve
    tstart = tic;
    [~, ~, K_cr, f_cr, u_cr] = LinearSolveLaplaceTTfG0_3D_cross( ...
        knot1n, knot2n, knot3n, ctpxn, ctpyn, ctpzn, ...
        bc_type, g_type, u_in, u_out, ffunction, tt_tol);
    tt_time = toc(tstart);
    fprintf("Torus - TT solve time = %.2e s\n", tt_time);
    %---- save results ----
    T = table(n1, n2, n3, tt_time, K_cr, f_cr, u_cr, 'VariableNames',{ ...
            'n1','n2','n3','TT_solve_time', ...
            'K_compression','f_compression','u_compression'});
    if idx == 1
       % Write header + first row
       writetable(T, outFile);
    else
       % Append just the data row
       writetable(T, outFile, 'WriteMode','append', 'WriteVariableNames',false);
    end
    %---- clear large arrays to free memory ----);

    % clear large vars before next iteration
    clear wn ctpxn ctpyn ctpzn knot1n knot2n knot3n
    clear K_cr f_cr u_cr tt_time
end

fprintf("\nAll done!  Results written to %s\n", outFile);