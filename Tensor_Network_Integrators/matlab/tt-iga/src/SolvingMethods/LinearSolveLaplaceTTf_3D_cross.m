function [u1tt, ustt, K_cr, f_cr, u_cr] = LinearSolveLaplaceTTf_3D_cross(knot1n, knot2n, knot3n,...
                                                         ctpxn, ctpyn, ctpzn,...
                                                         u_in, u_out, bc_type,...
                                                         ffunction,...
                                                         tt_tol)

[n1, n2, n3] = size(ctpxn);
%% Tensor decompose mass matrix (for dynamic problem)
%MT = GetMtt_3D(knot1n, knot2n, knot3n, ctpxn, ctpyn, ctpzn, tt_tol);
%fprintf('|MT-Mt| = %.2e \n',check_tt_error(Mt, MT));

%% Tensor decompose stiffness matrix
%tstart = tic;
fprintf("Computing tangent matrix K\n");
KT = GetKtt_3D_cross(knot1n, knot2n, knot3n, ctpxn, ctpyn, ctpzn, tt_tol);
%fprintf('Computing K cross time = %.2e \n', toc(tstart));
%fprintf('|KT-Kt| = %.2e \n',check_tt_error(Kt, KT));

%% Tensor decompose force vector
fprintf("Computing force vector f\n");
fT = Getftt_3D_cross(knot1n, knot2n, knot3n, ctpxn, ctpyn, ctpzn, ffunction, tt_tol);

%% compute K11tt and f1tt
G = core2cell(KT);
g = core2cell(fT);
if bc_type == "1"
    G11 = G;
    k1 = 1;
    k2 = n1;
    G11{1} = G11{1}(:, k1+1:k2-1, k1+1:k2-1, :);
    K11tt = cell2core(tt_matrix, G11);
    %%%
    g11 = g;
    g11{1} = g11{1}(:, k1+1:k2-1, :, :);
    f1tt = cell2core(tt_matrix, g11);
    f1tt = tt_tensor(f1tt);
elseif bc_type == "2"
    G11 = G;
    k1 = 1;
    k2 = n2;
    G11{2} = G11{2}(:, k1+1:k2-1, k1+1:k2-1, :);
    K11tt = cell2core(tt_matrix, G11);
    %%%
    g11 = g;
    g11{2} = g11{2}(:, k1+1:k2-1, :, :);
    f1tt = cell2core(tt_matrix, g11);
    f1tt = tt_tensor(f1tt);
elseif bc_type == "3"
    G11 = G;
    k1 = 1;
    k2 = n3;
    G11{3} = G11{3}(:, k1+1:k2-1, k1+1:k2-1, :);
    K11tt = cell2core(tt_matrix, G11);
    %%%
    g11 = g;
    g11{3} = g11{3}(:, k1+1:k2-1, :, :);
    f1tt = cell2core(tt_matrix, g11);
    f1tt = tt_tensor(f1tt);
elseif bc_type == "12"
    G11 = G;
    k11 = 1;
    k12 = n1;    
    G11{1} = G11{1}(:, k11+1:k12-1, k11+1:k12-1, :);
    k21 = 1;
    k22 = n2;
    G11{2} = G11{2}(:, k21+1:k22-1, k21+1:k22-1, :);
    K11tt = cell2core(tt_matrix, G11);
    %%%
    g11 = g;
    g11{1} = g11{1}(:, k11+1:k12-1, :, :);
    g11{2} = g11{2}(:, k21+1:k22-1, :, :);
    f1tt = cell2core(tt_matrix, g11);
    f1tt = tt_tensor(f1tt);
elseif bc_type == "23"
    G11 = G;
    k11 = 1;
    k12 = n2;    
    G11{2} = G11{2}(:, k11+1:k12-1, k11+1:k12-1, :);
    k21 = 1;
    k22 = n3;
    G11{3} = G11{3}(:, k21+1:k22-1, k21+1:k22-1, :);
    K11tt = cell2core(tt_matrix, G11);
    %%%
    g11 = g;
    g11{2} = g11{2}(:, k11+1:k12-1, :, :);
    g11{3} = g11{3}(:, k21+1:k22-1, :, :);
    f1tt = cell2core(tt_matrix, g11);
    f1tt = tt_tensor(f1tt);
elseif bc_type == "13"
    G11 = G;
    k11 = 1;
    k12 = n1;    
    G11{1} = G11{1}(:, k11+1:k12-1, k11+1:k12-1, :);
    k21 = 1;
    k22 = n3;
    G11{3} = G11{3}(:, k21+1:k22-1, k21+1:k22-1, :);
    K11tt = cell2core(tt_matrix, G11);
    %%%
    g11 = g;
    g11{1} = g11{1}(:, k11+1:k12-1, :, :);
    g11{3} = g11{3}(:, k21+1:k22-1, :, :);
    f1tt = cell2core(tt_matrix, g11);
    f1tt = tt_tensor(f1tt);
end

%fprintf('|f1 - f1tt| = %.2e \n',check_tt_error(f1,f1tt));
%% Solve
fprintf("AMEN solve\n");
%u1tt = amen_solve2(K11tt, f1tt, tt_tol);
u1tt = amen_solve2(K11tt, f1tt, tt_tol, ...
                  'resid_damp', 0.75, ...
                  'trunc_norm', 'fro', ...
                  'nswp', 300, ...
                  'rmax', 300, ...
                  'kickrank', 2);

u1tt = round(u1tt, tt_tol); 

%u1tt = tt_gmres_block(K11tt, f1tt, tt_tol);
%u1tt = dmrg_solve2(K11tt, f1tt, tt_tol);
%u1tt = dmrg_solve3(K11tt, f1tt, tt_tol);

ustt = get_utt2(bc_type, u_in, u_out, u1tt);

K_cr = 1/compress_ratio_tt(K11tt);
fprintf('compression ratio K = %.2e\n', K_cr);

f_cr = 1/compress_ratio_tt(f1tt);
fprintf('compression ratio f = %.2e\n', f_cr);

u_cr = 1/compress_ratio_tt(u1tt);
fprintf('compression ratio u = %.2e\n', u_cr);
