close all; clear; clc;
%%
load("Ring_Full_benchmark.mat");

%%
% Keep only valid (non-empty struct) entries
valid_idx = cellfun(@(c) isstruct(c) && ~isempty(c), R);
R_valid = R(valid_idx);  % Filtered cell array

% Extract fields using cellfun
rel_err        = cellfun(@(s) s.rel_err, R_valid);
tt_solve_time  = cellfun(@(s) s.tt_solve_time, R_valid);
tt_build_time  = cellfun(@(s) s.tt_build_time, R_valid);
fg_solve_time  = cellfun(@(s) s.fg_solve_time, R_valid);
fg_build_time  = cellfun(@(s) s.fg_build_time, R_valid);
cmpr           = cellfun(@(s) s.cmpr, R_valid);
% Compute DOFs from num_nodes
dof = cellfun(@(s) prod(s.num_nodes), R_valid);

% Derived quantities
tt_time = tt_build_time + tt_solve_time;
fg_time = fg_build_time + fg_solve_time;

%%
% === Plotting ===
figure;

% --- 1st subplot: DOF vs Relative Error
subplot(1,3,1);
loglog(dof, rel_err, 'o-k', 'LineWidth', 1.5);
xlabel('DOF'); ylabel('Relative Error');
title('DOF vs Relative Diff FG vs TT');
set(gca,'XScale','log','YScale','log');
grid on;

% --- 2nd subplot: DOF vs Timings
subplot(1,3,2);
hold on;
plot(dof, tt_solve_time, 'o-r', 'DisplayName', 'TT Solve');
plot(dof, tt_build_time, 'o-b', 'DisplayName', 'TT Build');
plot(dof, fg_solve_time, '*--r', 'DisplayName', 'FG Solve');
plot(dof, fg_build_time, '*--b', 'DisplayName', 'FG Build');
plot(dof, tt_time, 's-k', 'DisplayName', 'TT Total');
plot(dof, fg_time, '^-k', 'DisplayName', 'FG Total');
xlabel('DOF'); ylabel('Time (s)');
title('DOF vs Timing Components');
legend('Location', 'northwest');
set(gca,'XScale','log','YScale','log');
grid on;

% --- 3rd subplot: DOF vs Compression
subplot(1,3,3);
plot(dof, cmpr, 'd-k', 'LineWidth', 1.5);
xlabel('DOF'); ylabel('Compression Ratio');
title('DOF vs Compression');
set(gca,'XScale','log','YScale','log');
grid on;