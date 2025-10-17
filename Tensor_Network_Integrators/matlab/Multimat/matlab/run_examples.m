close all; clear;clc;
run setup_hydro_env.m;
global mesh_scaling;
%% specify problem here
probname = '2shocks'; %'sod','2shocks'
mesh_scaling = 0; %0;

%%
fprintf('\n--- Running TT solver ---\n');
tt_solver(probname, mesh_scaling);

fprintf('\n--- Running FG sparse solver ---\n');
fg_sparse_solver(probname, mesh_scaling);

fprintf('\n--- Loading results ---\n');
tt_file = sprintf('./%s_mesh_%d_TT.mat', probname, mesh_scaling);
fg_file = sprintf('./%s_mesh_%d_FG.mat', probname, mesh_scaling);

if ~isfile(tt_file)
  error('Missing TT result file: %s', tt_file);
end
if ~isfile(fg_file)
  error('Missing FG result file: %s', fg_file);
end

tt = load(tt_file);
fg = load(fg_file);

fprintf('\n--- Comparing results ---\n');
tt_assert_error(fg.rho,  tt.rho_tt,  -1, 'Density (rho)');
tt_assert_error(fg.ei,   tt.ei_tt,   -1, 'Internal Energy (ei)');
tt_assert_error(fg.pres, tt.pres_tt, -1, 'Pressure (pres)');

fprintf('\n--- Summary ---\n');
fprintf('TT elapsed time: %.2f seconds\n', tt.elapsed_time);
fprintf('FG elapsed time: %.2f seconds\n', fg.elapsed_time);
fprintf('TT cycles      : %d\n', tt.ncycle);
fprintf('FG cycles      : %d\n', fg.ncycle);
