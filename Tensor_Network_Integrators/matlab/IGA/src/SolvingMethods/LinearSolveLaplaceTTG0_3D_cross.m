function [u1tt, ustt, K_cr, f_cr, u_cr] = LinearSolveLaplaceTTG0_3D_cross(...
    knot1n, knot2n, knot3n, ...
    ctpxn_old, ctpyn_old, ctpzn_old, ...
    ubc, bc_type, g_type, u_in, u_out, tt_tol)

  %% 1) build the TT‐stiffness
  KT = GetKtt_3D_cross(...
    knot1n, knot2n, knot3n, ...
    ctpxn_old, ctpyn_old, ctpzn_old, ...
    tt_tol);
  %% Tensor decompose mass matrix (for dynamic problem)
  %MT = GetMtt_3D(knot1n, knot2n, knot3n, ctpxn_old, ctpxn_old, ctpxn_old, tt_tol);

  %% 2) unpack cores into a cell array
  Tcell = core2cell(KT);

  %% 3) figure out sizes and indices
  nVec  = size(ctpxn_old);        % [n1,n2,n3]
  idxG  = sscanf(g_type, '%d');   % 1,2 or 3  (G₀‐merge direction)
  idxB  = sscanf(bc_type, '%d');  % 1,2 or 3  (Dirichlet‐trim direction)

  %% 4) apply the G₀ “merge last→first” and drop last slice in dim idxG
  N = nVec(idxG);
  A = Tcell{idxG};
  A(:,1,:,:)    = A(:,1,:,:)    + A(:,N,:,:);  % merge face at N into face 1
  A(:,:,1,:)    = A(:,:,1,:)    + A(:,:,N,:,:);
  Tcell{idxG}  = A(:,1:N-1, 1:N-1, :);
  nVec(idxG)   = N - 1;

  %% 5) reassemble the modified stiffness
  KTm = cell2core(tt_matrix, Tcell);

  %% 6) build interior‐interior (K11) and interior‐boundary (K12) blocks
  Gcell = core2cell(KTm);
  k1 = 1;  k2 = nVec(idxB);

  % K11: trim out the two Dirichlet layers in dim idxB
  G11 = Gcell;
  G11{idxB} = G11{idxB}(:, k1+1:k2-1, k1+1:k2-1, :);
  K11tt    = cell2core(tt_matrix, G11);

  % K12: same interior rows, but pick only the two boundary slices
  G12 = Gcell;
  G12{idxB} = G12{idxB}(:, k1+1:k2-1, [k1,k2], :);
  K12tt     = cell2core(tt_matrix, G12);
  %K12tt     = round(K12tt, tt_tol);

  %% 7) pack the known‐bc vector into a TT‐tensor
  % after G₀‐merge the grid is (n1,n2,n3)=nVec; along idxB we have two faces
  u2sz = nVec;
  u2sz(idxB) = 2;
  UbcArr = reshape(ubc, u2sz);
  u2tt   = tt_tensor(UbcArr, tt_tol);
  u2tt   = round(u2tt, tt_tol);

  %% 8) form the reduced right‐hand‐side and solve
  f1tt = - K12tt * u2tt;
  f1tt   = round(f1tt, tt_tol);

  u1tt = amen_solve2(K11tt, f1tt, 1e-14, ...
                  'resid_damp', 0.75, ...
                  'trunc_norm', 'fro', ...
                  'nswp', 250, ...
                  'rmax', 250, ...
                  'kickrank', 2);

  u1tt = round(u1tt, tt_tol);

  %% 9) rebuild the full solution (including Dirichlet faces)
  ustt = get_utt(bc_type, u_in, u_out, u1tt);
  %% 10) Compress ratio of K
  K_cr = 1/compress_ratio_tt(K11tt);
  fprintf('compression ratio K = %.2e\n', K_cr);
  f_cr = 1/compress_ratio_tt(f1tt);
  fprintf('compression ratio f = %.2e\n', f_cr);
  u_cr = 1/compress_ratio_tt(u1tt);
  fprintf('compression ratio u = %.2e\n', u_cr);
end
