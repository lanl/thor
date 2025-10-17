function [u1tt, ustt, K_cr, f_cr, u_cr] = LinearSolveLaplaceTTfG0_3D_cross(knot1n, knot2n, knot3n, ...
                                                  ctpxn_old, ctpyn_old, ctpzn_old, ...
                                                  bc_type, ...
                                                  g_type, u_in, u_out, ...
                                                  ffunction, tt_tol)

  %% 1) Build TT‐stiffness and TT‐force
  fprintf("Computing tangent matrix K\n");
  KT = GetKtt_3D_cross(knot1n, knot2n, knot3n, ...
                 ctpxn_old, ctpyn_old, ctpzn_old, tt_tol);
  fprintf("Computing force vector f\n");
  fT = Getftt_3D_cross(knot1n, knot2n, knot3n, ...
                 ctpxn_old, ctpyn_old, ctpzn_old, ...
                 ffunction, tt_tol);
  %% 2) Unpack into cell arrays
  Tcell = core2cell(KT);
  tcell = core2cell(fT);

  %% 3) Prepare indexing
  sz    = size(ctpxn_old);      % [n1,n2,n3]
  nVec  = double(sz(:))';       % convert to row [n1,n2,n3]
  idxG  = sscanf(g_type,'%d');  % G₀‐direction index 1/2/3
  idxB  = sscanf(bc_type,'%d'); % BC‐trim index 1/2/3

  %% 4) Apply G₀‐merge + trim in direction idxG
  N = nVec(idxG);
  % 4a) stiffness core
  A = Tcell{idxG};
  A(:,1,:,:)    = A(:,1,:,:)    + A(:,N,:,:);  % merge last → first slice
  A(:,:,1,:)    = A(:,:,1,:)    + A(:,:,N,:,:);
  Tcell{idxG}  = A(:,1:N-1,1:N-1,:);
  % 4b) force core
  B = tcell{idxG};
  B(:,1,:,:)    = B(:,1,:,:)    + B(:,N,:,:);
  tcell{idxG}  = B(:,1:N-1,:,:);
  % shrink grid size
  nVec(idxG) = N - 1;

  %% 5) Trim out Dirichlet rows/cols in direction idxB
  k1 = 1;  k2 = nVec(idxB);
  % 5a) stiffness
  C = Tcell{idxB};
  Tcell{idxB} = C(:, k1+1:k2-1, k1+1:k2-1, :);
  % 5b) force
  D = tcell{idxB};
  tcell{idxB} = D(:, k1+1:k2-1, :, :);

  %% 6) Reassemble TT‐matrix and TT‐vector
  K11tt   = cell2core(tt_matrix, Tcell);
  f1core  = cell2core(tt_matrix, tcell);
  f1tt    = tt_tensor(f1core);

  %% 7) Solve and round
  fprintf("AMEN solve\n");
  u1tt = amen_solve2(K11tt, f1tt, tt_tol, ...
                     'resid_damp',0.75, ...
                     'trunc_norm','fro', ...
                     'nswp', 800, ...
                     'rmax', 600, ...
                     'kickrank',2);
  u1tt = round(u1tt, tt_tol);

  %% 8) Recover full solution with Dirichlet faces
  ustt = get_utt(bc_type, u_in, u_out, u1tt);

  %Compression ratio
  K_cr = 1/compress_ratio_tt(K11tt);
  fprintf('compression ratio K = %.2e \n', K_cr);
  f_cr = 1/compress_ratio_tt(f1tt);
  fprintf('compression ratio f = %.2e \n', f_cr);
  u_cr = 1/compress_ratio_tt(u1tt);
  fprintf('compression ratio u = %.2e \n', u_cr);
end
