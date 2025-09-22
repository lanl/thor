%% ------------------------------------------------------------------------
% Batch‐solve Laplace in TT for a list of mesh‐sizes, save timing & errors
close all; clear; clc;

addpath(genpath('../../matlab/utils/tt-toolbox/'));
addpath(genpath('../../matlab/utils/ttfunc/'));
addpath(genpath('../../matlab/tt-iga/src/'));

%% 0) User‐set constants
r_out      = 1.0;
r_int      = 0.5;
nm_list    = [2, 4, 8, 16];   % mesh refinements in x/y
nRun       = numel(nm_list);
p1 = 2; p2 = 2; p3 = 2;                 % spline degrees
g_type     = "1";                       % G₀‐merge direction
bc_type    = "2";                       % Dirichlet faces direction
u_in  = 1; u_out = 2;                   % boundary values
tt_tol     = 1e-12;                     % TT truncation tol

% output file
outFile   = "Ring3D_tt.csv";
if isfile(outFile)
  delete(outFile);
end
%% 1) Loop over mesh sizes
for idx = 1:nRun
  fprintf("===========================================================\n");
  nm = nm_list(idx);
  fprintf("Row %d/%d: ", idx, nRun);
  %---- geometry and IGA info ----
  [wn, ctpxn, ctpyn, ctpzn, knot1n, knot2n, knot3n] = ...
      RingFull222(r_int, r_out, nm, nm, nm);
  wn = ones(size(wn));   % B‐spline weights
  [n1, n2, n3] = size(wn);
  fprintf("Ring 3D tt - row %d/%d: \n", idx, nRun);
  fprintf("Size %d x %d x %d: \n", n1, n2, n3);
  [ctpxv, ctpyv, ctpzv, IndexE, ctpxe, ctpye, ctpze, we, m3D] = ...
    IGAinfo(ctpxn, ctpyn, ctpzn, wn, knot1n, knot2n, knot3n);

  %---- enforce G₀ continuity in parametric domain ----
  [ctpxn_old, ctpyn_old, ctpzn_old] = deal(ctpxn, ctpyn, ctpzn);
  [ctpxn, ctpyn, ctpzn, IndexE] = GDomain_Index1D(...
    ctpxn, ctpyn, ctpzn, IndexE, g_type, false);

  %---- Dirichlet BC setup ----
  [~,~,ubc] = DirichletBC3D_Laplace(...
    ctpxn, ctpyn, ctpzn, IndexE, bc_type, u_in, u_out, false);

  %---- TT solve ----
  tstart = tic;
  [~, ustt, Kcr, fcr, ucr] = LinearSolveLaplaceTTG0_3D_cross(...
    knot1n, knot2n, knot3n, ...
    ctpxn_old, ctpyn_old, ctpzn_old, ...
    ubc, bc_type, g_type, u_in, u_out, tt_tol);
  tt_time = toc(tstart);
  fprintf("Ring3D - TT solve time = %.2e s\n", tt_time);

  %---- compute L2 error ----
  uf   = full(ustt); 
  uett = uf(IndexE);
  [L2tt, ~] = CalH1L2_Ring(u_in, u_out, r_int, r_out, ...
                              uett, ctpxe, ctpye, ctpze, ...
                              knot1n, knot2n, knot3n, ...
                              p1, p2, p3, we, m3D);
  %if n1*n2*n3 < 1e6
  %   [L2tt, ~] = CalH1L2_Ring(u_in, u_out, r_int, r_out, ...
  %                            uett, ctpxe, ctpye, ctpze, ...
  %                            knot1n, knot2n, knot3n, ...
  %                            p1, p2, p3, we, m3D);
  %else
  %   [L2tt, ~] = CalH1L2_Ring_par(u_in, u_out, r_int, r_out, ...
  %                            uett, ctpxe, ctpye, ctpze, ...
  %                            knot1n, knot2n, knot3n, ...
  %                            p1, p2, p3, we, m3D);
  %end
  fprintf("Ring3Dtt - L2 norm       = %.2e\n", L2tt);

  %---- save results ----
  T = table(n1, n2, n3, tt_time, L2tt, Kcr, fcr, ucr, 'VariableNames',{ ...
            'n1', 'n2', 'n3', ...
            'TT_solve_time', ...
            'L2_TT', ...
            'K_compression', ...
            'f_compression', ...
            'u_compression'});
  if idx == 1
      % Write header + first row
      writetable(T, outFile);
  else
      % Append just the data row
      writetable(T, outFile, 'WriteMode','append', 'WriteVariableNames',false);
  end
  fprintf('Ring tt completed size (%d x %d x %d) in %d seconds \n', n1, n2, n3, tt_time);
  fprintf("=================================================================\n");
  %---- clear large arrays to free memory ----
  clear ustt Kcr fcr ucr uf uett ...
        wn ctpxv ctpyv ctpzv ctpxe ctpye ctpze we m3D ...
        ctpxn ctpyn ctpzn ctpxn_old ctpyn_old ctpzn_old ...
        IndexE T;
  % keep only loop counters, params, outCSV, nm_list, etc.
end
fprintf("\nAll done!  Results written to %s\n", outFile);
