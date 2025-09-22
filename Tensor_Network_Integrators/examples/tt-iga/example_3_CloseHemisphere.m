clear; close all; clc;

addpath(genpath('../../matlab/utils/tt-toolbox/'));
addpath(genpath('../../matlab/utils/ttfunc/'));
addpath(genpath('../../matlab/tt-iga/src/'));


%% Geometry & solver parameters
Rin       = 0.5;
Rout      = 1.0;
p1 = 2; p2 = 2; p3 = 2;
bc_type   = "2";             % as in your example
g_type    = "1";             % G⁰ continuity
u_in      = 0.0;
u_out     = 0.0;
ffunction = "sinpxsinpysinpz";
tt_tol    = 1e-10;

%% Refinement list
nm_list = [2, 4, 8, 16];
nRun    = numel(nm_list);

%% CSV setup
outFile = "Closehemisphere.csv";
if isfile(outFile)
    delete(outFile);
end

%% Loop over refinements
for idx = 1:nRun
    fprintf("===========================================================\n");
    nm = nm_list(idx);
    fprintf("Row %d/%d: ", idx, nRun);

    % build hemisphere geometry
    [wn, ctpxn, ctpyn, ctpzn, ...
        knot1n, knot2n, knot3n] = ...
      CloseHemiSphereFull221(Rin, Rout, nm, nm, nm);
    
    wn = ones(size(wn));  % use B-spline
    % grid size
    [n1, n2, n3] = size(ctpxn);
    fprintf("Close hemisphere tt - row %d/%d: \n", idx, nRun);
    fprintf("Size %d x %d x %d: \n", n1, n2, n3);

    % TT‐solve (includes BC & G⁰ internally)
    tstart = tic;
    [~, ~, K_cr, f_cr, u_cr] = LinearSolveLaplaceTTfG0_3D_cross( ...
        knot1n, knot2n, knot3n, ...
        ctpxn, ctpyn, ctpzn, ...
        bc_type, g_type, ...
        u_in, u_out, ffunction, ...
        tt_tol);
    tt_time = toc(tstart);
    fprintf("Close hemisphere - TT solve time = %.2e s\n", tt_time);
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