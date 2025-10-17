function obj_flat = tensor_to_flat_edge_flux(obj_tensor, dir, tol)
%TENSOR_TO_FLAT_EDGE_FLUX Convert per-face, per-material tensors to flat arrays.
%
% obj_flat = tensor_to_flat_edge_flux(obj_tensor, dir, tol)
%
% Inputs
%   obj_tensor : object (from ten_quantities_crossing_edge_3d) with fields:
%                dvol_for_edge(k,j,i,m), dmass_for_edge(k,j,i,m), dener_for_edge(k,j,i,m)
%                (optionally) nmat_for_edge(k,j,i) -- not required
%   dir        : sweep direction (1:x, 2:y, 3:z) used to define face set
%   tol        : threshold to treat a tensor entry as zero (default 0)
%
% Output
%   obj_flat   : struct with fields matching quantities_crossing_edge_3d layout:
%                nmat_for_edge(k,j,i)   int32
%                matid_for_edge(:)      int32 (linear list)
%                dvol_for_edge(:)       double (linear list)
%                dmass_for_edge(:)      double (linear list)
%                dener_for_edge(:)      double (linear list)
%
% The face traversal order and (k,j,i) indexing match both production routines.

  if nargin < 3, tol = 0; end
  global mesh;

  % Mesh sizes & face extents (match producers)
  ncell   = mesh.ncell_prob;
  nbdry   = mesh.nbdry_prob;
  ncell_ext  = ncell + 2*nbdry;
  ncell_bdry = ncell + nbdry;
  nnode_ext  = ncell_ext + 1;
  nnode_bdry = nnode_ext - nbdry;
  nmat = mesh.nmat_mesh;

  sizes_edge       = ncell_ext;     % faces array size in cell dims
  sizes_edge(dir)  = nnode_ext(dir);% but nodes along sweep dir

  imax = ncell_bdry(1); jmax = ncell_bdry(2); kmax = ncell_bdry(3);
  if dir == 1, imax = nnode_bdry(1); end
  if dir == 2, jmax = nnode_bdry(2); end
  if dir == 3, kmax = nnode_bdry(3); end

  % Validate tensor sizes
  tens_sz = size(obj_tensor.dvol_for_edge);
  expect_sz = [sizes_edge(3), sizes_edge(2), sizes_edge(1), nmat];
  if numel(tens_sz) ~= 4 || any(tens_sz ~= expect_sz)
    error('Tensor size mismatch. Got [%s], expected [%s].', ...
          num2str(tens_sz), num2str(expect_sz));
  end

  % Preallocate worst-case capacity: every face has all materials
  nfaces = sizes_edge(1) * sizes_edge(2) * sizes_edge(3);
  cap = nfaces * nmat;

  matid_for_edge = zeros(cap,1,'int32');
  dvol_for_edge  = zeros(cap,1);
  dmass_for_edge = zeros(cap,1);
  dener_for_edge = zeros(cap,1);
  nmat_for_edge  = zeros(sizes_edge(3), sizes_edge(2), sizes_edge(1), 'int32');

  loc_edge = 1; % next free slot in flat arrays

  % Walk faces in producer order
  for k = nbdry+1:kmax
    for j = nbdry+1:jmax
      for i = nbdry+1:imax
        % Pull the per-material tensors for this face
        v = squeeze(obj_tensor.dvol_for_edge (k,j,i,:));
        m = squeeze(obj_tensor.dmass_for_edge(k,j,i,:));
        e = squeeze(obj_tensor.dener_for_edge(k,j,i,:));

        if tol > 0
          active = find( (abs(v) > tol) | (abs(m) > tol) | (abs(e) > tol) );
        else
          active = find( v | m | e ); % any nonzero
        end

        if isempty(active)
          % no transport recorded at this face
          nmat_for_edge(k,j,i) = int32(0);
          continue;
        end

        % Deterministic ordering: ascending material id
        active = sort(active(:).'); %#ok<NASGU>

        % Write to flat arrays
        n_here = numel(active);
        nmat_for_edge(k,j,i) = int32(n_here);

        idxs = loc_edge : (loc_edge + n_here - 1);
        if idxs(end) > cap
          % Grow (should be rare; initial cap is conservative)
          grow = max(cap, nfaces); % double capacity
          matid_for_edge = [matid_for_edge; zeros(grow,1,'int32')];
          dvol_for_edge  = [dvol_for_edge;  zeros(grow,1)];
          dmass_for_edge = [dmass_for_edge; zeros(grow,1)];
          dener_for_edge = [dener_for_edge; zeros(grow,1)];
          cap = cap + grow;
        end

        mats = int32(active);
        matid_for_edge(idxs) = mats;
        % Use the same signed values as stored in the tensors
        dvol_for_edge(idxs)  = double(v(active));
        dmass_for_edge(idxs) = double(m(active));
        dener_for_edge(idxs) = double(e(active));

        loc_edge = loc_edge + n_here;
      end
    end
  end

  % Trim to actual length used
  last = loc_edge - 1;
  matid_for_edge = matid_for_edge(1:last);
  dvol_for_edge  = dvol_for_edge (1:last);
  dmass_for_edge = dmass_for_edge(1:last);
  dener_for_edge = dener_for_edge(1:last);

  % Pack output struct (mirrors quantities_crossing_edge_3d fields)
  obj_flat = struct();
  obj_flat.nmat_for_edge  = nmat_for_edge;
  obj_flat.matid_for_edge = matid_for_edge;
  obj_flat.dvol_for_edge  = dvol_for_edge;
  obj_flat.dmass_for_edge = dmass_for_edge;
  obj_flat.dener_for_edge = dener_for_edge;
end