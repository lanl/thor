function [IndexE, m3D] = N3DIndex_optimized(knot1, knot2, knot3)
  % Compute exact IGA 3D connectivity for any open knot vectors.

  % macro‐element boxes (not strictly needed for indexing, but returned anyway)
  m3D = Getmacro3D(knot1, knot2, knot3);

  % degrees
  [~, p1] = EvaluateKnot(knot1);
  [~, p2] = EvaluateKnot(knot2);
  [~, p3] = EvaluateKnot(knot3);

  % number of control points per direction
  n1 = length(knot1) - p1 - 1;
  n2 = length(knot2) - p2 - 1;
  n3 = length(knot3) - p3 - 1;

  % find knot‐span starts (element indices)
  el_x = find(diff(knot1) > 0);
  el_y = find(diff(knot2) > 0);
  el_z = find(diff(knot3) > 0);

  nx = numel(el_x);
  ny = numel(el_y);
  nz = numel(el_z);
  NOE = nx*ny*nz;
  nConn = (p1+1)*(p2+1)*(p3+1);

  IndexE = zeros(NOE, nConn);
  e = 0;
  %disp('Compute connectivity');
  %tic

  % loop in the same z→y→x order your reference used:
  for k = 1:nz
    for j = 1:ny
      for i = 1:nx
        e = e + 1;
        % pick the p+1 basis *ending* at each span
        ix = (el_x(i)-p1) : el_x(i);
        iy = (el_y(j)-p2) : el_y(j);
        iz = (el_z(k)-p3) : el_z(k);
        [II,JJ,KK] = ndgrid(ix, iy, iz);
        % global index, x fastest, then y, then z
        idx = II(:) + n1*(JJ(:)-1) + n1*n2*(KK(:)-1);
        IndexE(e,:) = idx.';
      end
    end
  end
%toc
end
