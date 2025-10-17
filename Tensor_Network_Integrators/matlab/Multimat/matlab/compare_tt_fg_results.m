function compare_tt_fg_results(probname, mesh_scaling)
    % compare_tt_fg_results - Compares TT and FG solver outputs.
    %
    % Usage:
    %   compare_tt_fg_results('2shocks', 0)
    %
    % Assumes you have TT and FG result files saved as:
    %   ./<probname>_mesh_<mesh_scaling>_TT.mat
    %   ./<probname>_mesh_<mesh_scaling>_FG.mat

    %% Construct file names
    tt_file = sprintf('./%s_mesh_%d_TT.mat', probname, mesh_scaling);
    fg_file = sprintf('./%s_mesh_%d_FG.mat', probname, mesh_scaling);

    %% Load data
    if ~isfile(tt_file)
        error('TT result file not found: %s', tt_file);
    end
    if ~isfile(fg_file)
        error('FG result file not found: %s', fg_file);
    end

    tt = load(tt_file);  % Contains rho_tt, ei_tt, pres_tt, etc.
    fg = load(fg_file);  % Contains rho, ei, pres, etc.

    %% Compare fields
    fprintf('Comparing TT and FG results for "%s" with mesh scaling %d:\n', probname, mesh_scaling);
    fprintf('--------------------------------------------------------------\n');
    tt_assert_error(fg.rho,     tt.rho_tt,   -1, 'Density (rho)');
    tt_assert_error(fg.ei,      tt.ei_tt,    -1, 'Internal Energy (ei)');
    tt_assert_error(fg.pres,    tt.pres_tt,  -1, 'Pressure (pres)');

    %% Compare runtime and cycles
    fprintf('\nRuntime comparison:\n');
    fprintf('  TT elapsed time: %.2f seconds\n', tt.elapsed_time);
    fprintf('  FG elapsed time: %.2f seconds\n', fg.elapsed_time);

    fprintf('\nCycle comparison:\n');
    fprintf('  TT cycles: %d\n', tt.ncycle);
    fprintf('  FG cycles: %d\n', fg.ncycle);

    fprintf('--------------------------------------------------------------\n');
end

%% Helper function
function compare_field(name, tt_field, fg_field)

fprintf('\n%s:\n', name);
fprintf('  Max Abs Diff : %e\n', max_abs);

end