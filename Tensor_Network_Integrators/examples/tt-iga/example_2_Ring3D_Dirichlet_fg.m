%% ------------------------------------------------------------------------
% Batch‐solve Laplace in FG for a list of mesh‐sizes, save timing & errors
close all; clear; clc;

addpath(genpath('../../matlab/utils/tt-toolbox/'));
addpath(genpath('../../matlab/utils/ttfunc/'));
addpath(genpath('../../matlab/tt-iga/src/'));


%% 0) User‐set constants
r_out      = 1.0;
r_int      = 0.5;
nm_list    = [2, 4, 8, 16];   % mesh refinements
nRun    = numel(nm_list);
p1 = 2; p2 = 2; p3 = 2;                 % spline degrees
g_type     = "1";                       % G₀‐merge direction
bc_type    = "2";                       % Dirichlet faces direction
u_in  = 1; u_out = 2;                   % boundary values

% output file
outFile   = "Ring3D_fg.csv";
if isfile(outFile)
  delete(outFile);
end

%% 1) Loop over mesh sizes
for idx = 1:nRun
  fprintf("===========================================================\n");
  nm = nm_list(idx);
  
  %---- geometry and IGA info ----
  [wn, ctpxn, ctpyn, ctpzn, knot1n, knot2n, knot3n] = ...
      RingFull222(r_int, r_out, nm, nm, nm);
  wn = ones(size(wn));   % B‐spline weights
  [n1, n2, n3] = size(wn);
  fprintf("\n=== Running nm = %d × %d x %d mesh ===\n", n1, n2, n3);
  [ctpxv, ctpyv, ctpzv, IndexE, ctpxe, ctpye, ctpze, we, m3D] = ...
    IGAinfo(ctpxn, ctpyn, ctpzn, wn, knot1n, knot2n, knot3n);

  %---- enforce G₀ continuity in parametric domain ----
  [ctpxn_old, ctpyn_old, ctpzn_old] = deal(ctpxn, ctpyn, ctpzn);
  [ctpxn, ctpyn, ctpzn, IndexE] = GDomain_Index1D(...
    ctpxn, ctpyn, ctpzn, IndexE, g_type, false);

  %---- Dirichlet BC setup ----
  [uindex1, uindex2, ubc] = DirichletBC3D_Laplace(...
    ctpxn, ctpyn, ctpzn, IndexE, bc_type, u_in, u_out, false);

  %---- TT solve ----
  %[u1tt, ustt, Kcr] = LinearSolveLaplaceTTG0_3D_opt(...
  %  knot1n, knot2n, knot3n, ...
  %  ctpxn_old, ctpyn_old, ctpzn_old, ...
  %  ubc, bc_type, g_type, u_in, u_out, tt_tol);
  %tt_solve_time = toc;
  %fprintf("  TT solve time = %.2e s\n", tt_solve_time);

  %---- FG (sparse) solve ----
  tstart = tic;
  [ueL, ~, ~, ~] = LinearSolveLaplace_3D(ctpxe, ctpye, ctpze, ...
                                          knot1n, knot2n, knot3n, ...
                                          p1, p2, p3, ...
                                          we, m3D, IndexE, ...
                                          uindex1, uindex2, ubc);
  fg_solve_time = toc(tstart);
  fprintf('FG solve time = %.2e \n', fg_solve_time);

  %---- compute L2 error ----
  [L2fg, ~] = CalH1L2_Ring(...
    u_in, u_out, r_int, r_out, ...
    ueL, ctpxe, ctpye, ctpze, ...
    knot1n, knot2n, knot3n, ...
    p1, p2, p3, we, m3D);
  fprintf("  L2 norm       = %.2e\n", L2fg);


  %---- save results ----
  T = table(n1, n2, n3, fg_solve_time, L2fg, 'VariableNames',{ ...
            'n1', 'n2', 'n3', ...
            'FG_solve_time', ...
            'L2_fg'});
  if idx == 1
      % Write header + first row
      writetable(T, outFile);
  else
      % Append just the data row
      writetable(T, outFile, 'WriteMode','append', 'WriteVariableNames',false);
  end
  %---- clear large arrays to free memory ----
  clear ueL uindex1 uindex2...
        wn ctpxv ctpyv ctpzv ctpxe ctpye ctpze we m3D ...
        ctpxn ctpyn ctpzn ctpxn_old ctpyn_old ctpzn_old ...
        IndexE T;
  % keep only loop counters, params, outCSV, nm_list, etc.
end

fprintf("\nAll done!  Results written to %s\n", outFile);
