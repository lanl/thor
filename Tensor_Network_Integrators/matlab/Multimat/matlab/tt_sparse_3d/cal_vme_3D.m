function X = cal_vme_3D(k, j, i, mesh, mat, ttadv, dir, fracten, ...
                        nbdry, kmax, jmax, imax, nmat, vol_cell, dx, xl_prob)
% CAL_VME_3D  Compute signed per-material fluxes (dvol, dmass, dener) for one face.
% Matches ten_quantities_crossing_edge_3d’s logic and sign conventions.
% Returns: X = [dvol(1:nmat), dmass(1:nmat), dener(1:nmat)]

  % Preallocate per-material outputs (signed)
  dvol_for_edge  = zeros(1, nmat);
  dmass_for_edge = zeros(1, nmat);
  dener_for_edge = zeros(1, nmat);

  % --- Skip if in active cell region (as requested) ---
if (k >= nbdry+1 && k <= kmax) && ...
    (j >= nbdry+1 && j <= jmax) && ...
    (i >= nbdry+1 && i <= imax)
  %   X = [dvol_for_edge, dmass_for_edge, dener_for_edge];
  %   return;
  % end

  % Fractional Courant at this face (assumed: vel*dt/dx(dir))
  frac = fracten(k, j, i);
  if abs(frac) <= ttadv.small
    X = [dvol_for_edge, dmass_for_edge, dener_for_edge];
    return;
  end

  % Distance traveled along 'dir'
  dist = frac * dx(dir);

  % Upwind/source cell indices (same as tensor routine)
  if frac > 0
    kl = k; jl = j; il = i;
    if dir == 1, il = i - 1; end
    if dir == 2, jl = j - 1; end
    if dir == 3, kl = k - 1; end
  else
    kl = k; jl = j; il = i;
  end

  % Active materials in upwind cell
  mats = [];
  for m = 1:nmat
    if mesh.vf_3dmat(kl, jl, il, m) >= mat.vfmin
      mats(end+1) = m; %#ok<AGROW>
    end
  end
  if isempty(mats)
    X = [dvol_for_edge, dmass_for_edge, dener_for_edge];
    return;
  end

  sign_frac = 2*(frac > 0) - 1;  % +1 for +dir, -1 for -dir

  if numel(mats) == 1
    % -------- Single-material fast path --------
    m = mats(1);
    dvol  = sign_frac * vol_cell * abs(frac);
    dmass = dvol * mesh.rho_3dmat(kl, jl, il, m);
    dener = dvol * mesh.ei_3dmat (kl, jl, il, m);

    dvol_for_edge(m)  = dvol;
    dmass_for_edge(m) = dmass;
    dener_for_edge(m) = dener;

  else
    % -------- Multi-material geometric advection --------
    ijk_up  = [il, jl, kl];
    xl_cell = xl_prob + (ijk_up - nbdry - 1) .* dx;

    inward_norm = [0.0, 0.0, 0.0];
    xl_slab = xl_cell; xr_slab = xl_cell;

    if frac > 0
      inward_norm(dir) = 1.0;         % into slab from edge plane
      slab_faceid = 2*dir - 1;        % 1,3,5
      xl_slab(dir) = xl_cell(dir) + dx(dir) - abs(dist);
      xr_slab(dir) = xl_cell(dir) + dx(dir);
    else
      inward_norm(dir) = -1.0;        % into slab from edge plane
      slab_faceid = 2*dir;            % 2,4,6
      xl_slab(dir) = xl_cell(dir);
      xr_slab(dir) = xl_cell(dir) + abs(dist);
    end

    % Build polyhedron for upwind cell
    [nnode_tot, coords_tot, nface_for_mpoly, ...
     nnode_for_face_ea_mpoly, nodelist_for_face_ea_mpoly] = ...
       ttadv.build_mpoly3d(mesh, ijk_up);

    % IMPORTANT: pass xl_cell (cell lower corner), not xl_slab
    [nmat_adv, matid_adv, vol_adv, mass_adv, ener_adv] = ...
      mat.advect3d(mesh, ijk_up, xl_cell, dx, numel(mats), mats, ...
                   xl_slab, xr_slab, inward_norm, slab_faceid, ...
                   nnode_tot, coords_tot, nface_for_mpoly, ...
                   nnode_for_face_ea_mpoly, nodelist_for_face_ea_mpoly);

    % Write signed per-material results
    for p = 1:nmat_adv
      m = matid_adv(p);
      dvol_for_edge(m)  = sign_frac * vol_adv(p);
      dmass_for_edge(m) = sign_frac * mass_adv(p);
      dener_for_edge(m) = sign_frac * ener_adv(p);
    end
  end
end
  % Output vector
  X = [dvol_for_edge, dmass_for_edge, dener_for_edge];
end





%% %%%
% function X = cal_vme_3D(k, j, i, mesh, mat, ttadv, dir, fracten, ...
%                         nbdry, kmax, jmax, imax, nmat, vol_cell, dx, xl_prob)
% 
% % Preallocate assuming at most 2 materials
% dvol_for_edge  = [0, 0];
% dmass_for_edge = [0, 0];
% dener_for_edge = [0, 0];
% 
% % Skip if in active cell region
% if (k >= nbdry+1 && k <= kmax) && ...
%    (j >= nbdry+1 && j <= jmax) && ...
%    (i >= nbdry+1 && i <= imax)
% %   X = [dvol_for_edge, dmass_for_edge, dener_for_edge];
% %   return;
% % end
% 
% frac = fracten(k,j,i);
% 
% if abs(frac) > ttadv.small
%   dist = frac * dx(dir);  % define before using
% 
%   if frac > 0
%     kl = k; jl = j; il = i;
%     if dir == 1, il = i-1; end
%     if dir == 2, jl = j-1; end
%     if dir == 3, kl = k-1; end
%   else
%     kl = k; jl = j; il = i;
%   end
% 
%   % Identify materials in source cell
%   mats = [];
%   for m = 1:nmat
%     if mesh.vf_3dmat(kl, jl, il, m) >= mat.vfmin
%       mats(end+1) = m;
%     end
%   end
% 
%   if isempty(mats)
%     X = [dvol_for_edge, dmass_for_edge, dener_for_edge];
%     return;
%   end
% 
%   if length(mats) == 1
%     m = mats(1);
%     dvol = vol_cell * abs(frac);
%     dmass = dvol * mesh.rho_3dmat(kl,jl,il,m);
%     dener = dvol * mesh.ei_3dmat(kl,jl,il,m);
%     dvol_for_edge(m)  = dvol;
%     dmass_for_edge(m) = dmass;
%     dener_for_edge(m) = dener;
%   else
%     ijk = [il, jl, kl];
% 
%     if frac > 0
%       xl_slab = xl_prob + (ijk - nbdry - 1).*dx;
%       xr_slab = xl_slab;
%       xl_slab(dir) = xl_slab(dir) + dx(dir) - abs(dist);
%       inward_norm = zeros(1,3); inward_norm(dir) = 1.0;
%       slab_faceid = 2*dir - 1;
%     else
%       xl_slab = xl_prob + (ijk - nbdry - 1).*dx;
%       xr_slab = xl_slab;
%       xr_slab(dir) = xr_slab(dir) + abs(dist);
%       inward_norm = zeros(1,3); inward_norm(dir) = -1.0;
%       slab_faceid = 2*dir;
%     end
% 
%     [nnode_tot, coords_tot, nface_for_mpoly, ...
%      nnode_for_face_ea_mpoly, nodelist_for_face_ea_mpoly] = ...
%       ttadv.build_mpoly3d(mesh, ijk);
% 
%     [nmat_adv, matid_adv, vol_adv, mass_adv, ener_adv] = ...
%       mat.advect3d(mesh, ijk, xl_slab, dx, length(mats), mats, ...
%       xl_slab, xr_slab, inward_norm, slab_faceid, ...
%       nnode_tot, coords_tot, nface_for_mpoly, ...
%       nnode_for_face_ea_mpoly, nodelist_for_face_ea_mpoly);
% 
%     for p = 1:nmat_adv
%       m = matid_adv(p);
%       dvol_for_edge(m)  = vol_adv(p) * sign(frac);
%       dmass_for_edge(m) = mass_adv(p) * sign(frac);
%       dener_for_edge(m) = ener_adv(p) * sign(frac);
%     end
%   end
% end
% end % active region
% 
% % Return final vector
% X = [dvol_for_edge, dmass_for_edge, dener_for_edge];
% 
% end
% 
% % function X = cal_vme_3D(k, j, i, mesh,mat, ttadv, dir, fracten, nbdry, ...
% %   kmax, jmax, imax,nmat,vol_cell)
% % 
% % %% preallocate
% % dvol_for_edge  = [0, 0];
% % dmass_for_edge= [0, 0];
% % dener_for_edge = [0, 0];
% % 
% % if (k >= nbdry+1 && k <= kmax) && ...
% %     (j >= nbdry+1 && j <= jmax) && ...
% %     (i >= nbdry+1 && i <= imax)
% %   X = [dvol_for_edge,dmass_for_edge,dener_for_edge];
% %   return;
% % end
% % frac = fracten(k,j,i);
% % %%
% % if abs(frac) > ttadv.small
% % 
% %   if frac > 0
% %     kl = k; jl = j; il = i;
% %     if dir == 1, il = i-1; end
% %     if dir == 2, jl = j-1; end
% %     if dir == 3, kl = k-1; end
% %   else
% %     kl = k; jl = j; il = i;
% %   end
% % 
% %   % Identify materials at source cell
% %   mats = [];
% %   for m = 1:nmat
% %     if mesh.vf_3dmat(kl, jl, il, m) >= mat.vfmin
% %       mats(end+1) = m;
% %     end
% %   end
% % 
% %   if isempty(mats)
% %    X = [dvol_for_edge,dmass_for_edge,dener_for_edge];
% %    return;
% %   end
% % 
% %   if length(mats) == 1
% %     % Simple case: only one material
% %     m = mats(1);
% %     % obj.matid_for_edge(k,j,i,m) = 1;
% %     dvol = vol_cell * abs(frac);
% %     dmass = dvol * mesh.rho_3dmat(kl,jl,il,m);
% %     dener = dvol * mesh.ei_3dmat(kl,jl,il,m);
% %     dvol_for_edge(m)  = dvol;
% %     dmass_for_edge(m) = dmass;
% %     dener_for_edge(m) = dener;
% %   else
% %     % Mixed cell, need polygon advection
% %     ijk = [il, jl, kl];
% % 
% %     if frac > 0
% %       xl_slab = (xl_prob + (ijk - nbdry - 1).*dx);
% %       xr_slab = xl_slab;
% %       xl_slab(dir) = xl_slab(dir) + dx(dir) - abs(dist);
% %       inward_norm = zeros(1,3);
% %       inward_norm(dir) = 1.0;
% %       slab_faceid = 2*dir-1;
% %     else
% %       xl_slab = (xl_prob + (ijk - nbdry - 1).*dx);
% %       xr_slab = xl_slab;
% %       xr_slab(dir) = xr_slab(dir) + abs(dist);
% %       inward_norm = zeros(1,3);
% %       inward_norm(dir) = -1.0;
% %       slab_faceid = 2*dir;
% %     end
% % 
% %     [nnode_tot, coords_tot, nface_for_mpoly, ...
% %       nnode_for_face_ea_mpoly, nodelist_for_face_ea_mpoly] = ...
% %       ttadv.build_mpoly3d(mesh, ijk);
% % 
% %     [nmat_adv, matid_adv, vol_adv, mass_adv, ener_adv] = ...
% %       mat.advect3d(mesh, ijk, xl_slab, dx, length(mats), mats, ...
% %       xl_slab, xr_slab, inward_norm, slab_faceid, ...
% %       nnode_tot, coords_tot, nface_for_mpoly, ...
% %       nnode_for_face_ea_mpoly, nodelist_for_face_ea_mpoly);
% % 
% %     % Store results
% %     for p = 1:nmat_adv
% %       m = matid_adv(p);
% %       dvol_for_edge(m)  = vol_adv(p) * sign(frac);
% %       dmass_for_edge(m) = mass_adv(p) * sign(frac);
% %       dener_for_edge(m) = ener_adv(p) * sign(frac);
% %     end
% %   end
% % end % frac
% % 
% % X = [dvol_for_edge,dmass_for_edge,dener_for_edge];
% % 
% % end % function
