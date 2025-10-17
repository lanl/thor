function obj = apply_face_fluxes_to_cell(obj, K, J, I, dir, nmat, sizes_edge)
% APPLY_FACE_FLUXES_TO_CELL
% Update ONLY cell (K,J,I) using its two adjacent faces along 'dir'.
% - sizes_edge: [nx, ny, nz] for the faces array (as in your code)
% - obj must have dvol_for_edge/dmass_for_edge/dener_for_edge of size
%   (k=sz(3), j=sz(2), i=sz(1), m=nmat) and per-cell arrays vol/mass/ener.

% Map cell -> minus/plus face indices (kE,jE,iE order)
switch dir
  case 1  % x-sweep: faces at i and i+1
    kM = K; jM = J; iM = I;      % minus face
    kP = K; jP = J; iP = I + 1;  % plus  face
  case 2  % y-sweep: faces at j and j+1
    kM = K; jM = J;      iM = I;
    kP = K; jP = J + 1;  iP = I;
  case 3  % z-sweep: faces at k and k+1
    kM = K; jM = J; iM = I;
    kP = K + 1; jP = J; iP = I;
end

% ----- MINUS FACE -----
[sgnM, vM, mM, eM] = face_flux_at(obj, kM, jM, iM, nmat, sizes_edge);
if sgnM ~= 0
  gainM = (sgnM > 0); % +: (cell-1)->(cell) => current cell gains
  for mm = 1:nmat
    if vM(mm)==0 && mM(mm)==0 && eM(mm)==0, continue; end
    dv = abs(vM(mm)); dm = abs(mM(mm)); de = abs(eM(mm));
    if gainM
      obj.vol_for_cell (K,J,I,mm) = obj.vol_for_cell (K,J,I,mm) + dv;
      obj.mass_for_cell(K,J,I,mm) = obj.mass_for_cell(K,J,I,mm) + dm;
      obj.ener_for_cell(K,J,I,mm) = obj.ener_for_cell(K,J,I,mm) + de;
    else
      obj.vol_for_cell (K,J,I,mm) = obj.vol_for_cell (K,J,I,mm) - dv;
      obj.mass_for_cell(K,J,I,mm) = obj.mass_for_cell(K,J,I,mm) - dm;
      obj.ener_for_cell(K,J,I,mm) = obj.ener_for_cell(K,J,I,mm) - de;
    end
  end
end

% ----- PLUS FACE -----
[sgnP, vP, mP, eP] = face_flux_at(obj, kP, jP, iP, nmat, sizes_edge);
if sgnP ~= 0
  gainP = (sgnP < 0); % -: (cell+1)->(cell) => current cell gains
  for mm = 1:nmat
    if vP(mm)==0 && mP(mm)==0 && eP(mm)==0, continue; end
    dv = abs(vP(mm)); dm = abs(mP(mm)); de = abs(eP(mm));
    if gainP
      obj.vol_for_cell (K,J,I,mm) = obj.vol_for_cell (K,J,I,mm) + dv;
      obj.mass_for_cell(K,J,I,mm) = obj.mass_for_cell(K,J,I,mm) + dm;
      obj.ener_for_cell(K,J,I,mm) = obj.ener_for_cell(K,J,I,mm) + de;
    else
      obj.vol_for_cell (K,J,I,mm) = obj.vol_for_cell (K,J,I,mm) - dv;
      obj.mass_for_cell(K,J,I,mm) = obj.mass_for_cell(K,J,I,mm) - dm;
      obj.ener_for_cell(K,J,I,mm) = obj.ener_for_cell(K,J,I,mm) - de;
    end
  end
end
end

% ---------- subfunction ----------
function [sgn, v, m, e] = face_flux_at(obj, kE, jE, iE, nmat, sizes_edge)
% Bounds (note: faces are indexed as k=sz(3), j=sz(2), i=sz(1))
if kE < 1 || jE < 1 || iE < 1 || ...
    kE > sizes_edge(3) || jE > sizes_edge(2) || iE > sizes_edge(1)
  sgn = 0; v = []; m = []; e = []; return;
end
v = squeeze(obj.dvol_for_edge (kE,jE,iE,1:nmat));
m = squeeze(obj.dmass_for_edge(kE,jE,iE,1:nmat));
e = squeeze(obj.dener_for_edge(kE,jE,iE,1:nmat));
if ~any(v(:) | m(:) | e(:)), sgn = 0; return; end

% "first-nonzero" sign rule to match flat code
nz = find(v ~= 0, 1, 'first');
if isempty(nz)
  nz = find(m ~= 0, 1, 'first');
  if isempty(nz)
    nz = find(e ~= 0, 1, 'first');
    if isempty(nz), sgn = 0; return; end
    sgn = sign(e(nz));
  else
    sgn = sign(m(nz));
  end
else
  sgn = sign(v(nz));
end
end